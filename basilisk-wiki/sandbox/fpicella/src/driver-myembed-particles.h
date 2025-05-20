/*
# Rigid particles in Stokes flow
using myembed (Arthur)
and particle (Antoon)
*/

#include "fpicella/src/periodic-shift-treatment.h"
#include "fpicella/src/compute_embed_color_force_torque_RBM.h"
/*
Required definitions
I set them like this if undefined elsewhere */

#ifndef NPARTICLES
#define NPARTICLES 1
#endif

static FILE *singleParticleFile[NPARTICLES] = {NULL}; // for the moment only one set of particles...?

/*
Buffer fields, required for definition of _global_ cs and fs fractions...
*/
scalar csLOCAL[];
face vector fsLOCAL[];

/*
Some practical shortcuts
*/
#define circle(px,py,pr) (sq (PS(x,px)) + sq (PS(y,py)) - sq (pr))
#define PARTICLE circle(p().x,p().y,p().r)
#define COLOR (sq (PS(COORD.x,p().x)) + sq (PS(COORD.y,p().y)) - sq (p().r)) // more complex, hate foreach_dimension() :D

#define THETA atan2(PS(COORD.y,p().y),PS(COORD.x,p().x))
/*
Simple way to account for variable particle orientation...
*/
#define TS (THETA - p().theta.z)

#define sinTheta     sin(p().theta.z)
#define cosTheta     cos(p().theta.z)

/**
Modify the structure of the particle pointer*/
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; \
double r; \
coord theta; \
coord omega; \
coord F; \
coord T; \
coord B; \
coord M; \
double thrust; double alpha; double beta; 
/**
F, T, translational-rotational forces MEASURED on particle. (i.e. HYDRO).
B, M, translational-rotational BODY forces apllied on particle.*/

#include "fpicella/src/myembed-particles.h"

int alternating_series(int n) {
    if (n == 0) return 0; // Base case
    int value = (n + 1) / 2;
    return (n % 2 == 1) ? value : -value;
}


void locate_particle_in_zero(){
  int _l_particle = 0;
    while (pn[_l_particle] != terminate_int) {
      for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
/*
First, barebones version, all particles have the same radius...
...but provided the numerical method improved, I could play around with 
different shapes, size...
For the moment, I stick with cylinders and spheres
*/
        p().r = RADIUS;

        p().x = +3.*p().r*alternating_series(_j_particle);
        p().y = +0.*p().r*alternating_series(_j_particle);
        p().z = 0.;
        p().theta.z=+_j_particle*M_PI*1.;
				p().u.x = 0.0;
				p().u.y = 0.0;
				p().omega.x = 0.;
				p().omega.y = 0.;
				p().omega.z = 0.;
				p().theta.x = 0.;
				p().theta.y = 0.;
				p().theta.z = 0.;
				// Body force (i.e. sedimentation?)
				p().B.x = 0.;
				p().B.y = 0.;
				p().B.z = 0.;
      }
      _l_particle++;
    }
}

/*Set list of particles to work with.
For the moment, I've fot a single set.*/
Particles ParticleList; // // Call the particles!
/*Initialize particles */
void initialize_particles ()
{
  ParticleList = init_tp_circle(NPARTICLES);
  locate_particle_in_zero();
}


///*
//Default particle initialization
//*/
//event init (i = 0){
//	initialize_particles();
//}


// Particle output, super-compact (and working in serial!) way :)
// FP, 20250212 11h37
event output_particle_initialize(i=0){ // first iteration, define name and open files...
  foreach_particle(){
    char filename[100];
    sprintf(filename, "particle_%03d.dat", _j_particle);
    singleParticleFile[_j_particle] = fopen(filename,"w");
  }
}
#define FORMAT_SPEC_7 ("%+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n") // to avoid writing every time it...
event output_particle(i++){
  foreach_particle(){
    fprintf(singleParticleFile[_j_particle],FORMAT_SPEC_7,t,p().x,p().y,p().theta.z,p().u.x,p().u.y,p().omega.z);
    fflush(singleParticleFile[_j_particle]);
  }
}

/*
### Compute volume fractions associated to the presence of the particle
Idea here, is to compute the LOCAL fractions, and to add them to the
full problem.
*/
void compute_particle_fractions(){
/**
	Updated version: each time reset ALL cs and fs fields to 1.
	This is done so to avoid, in the case container is treated using mask
	to _decapitate_ particle's fractions by successive multiplication of
	non 1 or 0 values on the boundaries. */
	foreach()
		cs[] = 1.;
	foreach_face()
		fs.x[] = 1.;
/**
	Get back to the original implementation*/
	foreach_particle(){
		solid(csLOCAL,fsLOCAL,PARTICLE);
		foreach()
			cs[] *= csLOCAL[];
		foreach_face()
			fs.x[] *= fsLOCAL.x[];
	}
}

/**
## Fluid-solid coupling

The following no-slip Dirichlet boundary condition for velocity and
hommogeneous Neumann boundary condition for pressure on the discrete
rigid boundary $\delta \Gamma_{i,\Delta}$ allow us to couple the
motion of the fluid and the discrete rigid body $\Gamma_{i,\Delta}$:
$$
\left\{
\begin{aligned}
&
\mathbf{u} = \mathbf{u}_{\Gamma,i} = \mathbf{u_{\Gamma,i}} +
\omega_{\Gamma,i} \times \left(\mathbf{x} - \mathbf{x}_{\Gamma,i}\right)
\\
&
{\nabla}_{\Gamma,i} p = 0.
\end{aligned}
\right.
$$
The homogeneous Neumann boundary condition is suitable for a fixed
rigid body and we have found that using it with a moving rigid body
does not significantly affect the computed solution.

#### No-slip Dirichlet boundary condition for velocity

The function *velocity_noslip_x()* computes the previously defined
no-slip boundary condition for the $x$-component of the velocity
$\mathbf{u}$. */

foreach_dimension()
static inline double velocity_noslip_x (Point point,
					double xc, double yc, double zc)
					// xc, yc and zc are coordinates of the cell's center
{
  assert (cs[] > 0. && cs[] < 1.);

	foreach_particle(){
		coord COORD = {xc,yc,zc}; // required to properly compute COLOR function
															// owing to a logic in foreach_dimension()
		if (COLOR>-0.5 && COLOR<0.5){
      /**
      We first compute the relative position $\mathbf{r_{i}} =
      \mathbf{x} - \mathbf{x}_{\Gamma,i}$. */ 
      // The coordinate x,y,z are not permuted with foreach_dimension()
      coord r = {xc, yc, zc};
      foreach_dimension() {
	r.x -= p().x;
	if (Period.x) {
	  if (fabs (r.x) > fabs (r.x + (L0)))
	    r.x += (L0);
	  if (fabs (r.x) > fabs (r.x - (L0)))
	    r.x -= (L0);
	}
      }
  
		/*
		We then compute the surface velocity as no-slip RIGID BODY MOTION
		*/
#if dimension == 2
      coord sgn = {-1, 1};
      return (p().u.x) + sgn.x*p().omega.x*(r.y); // in 2D, omega.x=omega.y=omega.z
#else // dimension == 3
    return (p().u.x) + p().omega.y*(r.z) - p().omega.z*(r.y);
#endif // dimension
		}
	}
  // Other fixed embedded boundaries (u=0)
#ifndef imposedFlow
  return 0.;
#else
	return imposedU.x;
#endif
/*
	Some other BC on _fixed_ embedded boundaries? Just change them here :)
*/
}


/*
Non-zero neumann bc for pressure.
Taken from Ghigo's [implementation]()
and adapted for Stokes configuration.
*/
/**
#### Neumann boundary condition for pressure

The function *pressure_gradient()* returns in *da* the contribution of
all the components of the particle acceleration and viscous stresses
to the pressure gradient. */

/*
For the case of Steady Stokes flow, the part owing to acceleration will be
accounted to zero. Viscous stresses instead will be set on.
*/

static inline void pressure_acceleration (Point point,
					  double xc, double yc, double zc,
					  coord * da)
{
  assert (cs[] > 0. && cs[] < 1.);
// On all embedded cells that do not belong to a particle
// i.e. to a container, channel...
	foreach_dimension()
		da->x = 0.;
// if not...
	foreach_particle(){
		coord COORD = {xc,yc,zc}; // required to properly compute COLOR function
															// owing to a logic in foreach_dimension()
		if (COLOR>-0.5 && COLOR<0.5){
	    /**
	    We first compute the acceleration of the particle *gu*. */
	      
	    coord gu;
	    foreach_dimension()
	      gu.x = 0.;
	    /**
	    We then compute the viscous contribution. */
	      
	    coord gmu;
	    foreach_dimension()
	      gmu.x = 0.;
	    
	     foreach_dimension() { 
	       scalar s = u.x; 
	       double a = 0.; 
	       foreach_dimension() 
	     	a += (mu.x[1]*face_gradient_x (s, 1) - mu.x[0]*face_gradient_x (s, 0))/Delta; 
	       double b, c = embed_flux (point, u.x, mu, &b); 
	       a += (b + c*u.x[]); 
	       gmu.x = -a; 
	     } 
	          
	    /**
	    Finally, we combine both contributions. */
	  
	    foreach_dimension()
	      da->x = (gu.x + gmu.x); 
		}
	}
}

foreach_dimension()
static inline double pressure_acceleration_x (Point point,
					      double xc, double yc, double zc)
{
  coord da;
  pressure_acceleration (point, xc, yc, zc, &da);
  return da.x;
}

/**
The following function finally computes the Neumann boundary condition
for pressure, where the inward unit normal $\mathbf{n}_{\Gamma}$ is
pointing from fluid to solid. To switch to a homogenous Neumann
boundary condition, the user can set the scalar attribute
*neumann_zero* to true. The default is false. */

attribute {
  bool neumann_zero;
}

static inline double pressure_neumann (Point point,
//				       coord dup, coord dwp,
//				       coord wp, coord pp,
				       double xc, double yc, double zc)
{
  coord da = {0., 0., 0.};
  //pressure_acceleration (point, dup, dwp, wp, pp, xc, yc, zc, &da);
  pressure_acceleration (point, xc, yc, zc, &da);

  coord b, n;
  embed_geometry (point, &b, &n);
  
  double dpdn = 0.;
  foreach_dimension()
    dpdn += (-da.x)*n.x;
  dpdn *= rho[]/(cs[] + SEPS);
  return dpdn;
}






void hydro_forces_torques()
{
	foreach_particle(){
		fraction(csLOCAL,PARTICLE);
		// allocate some auxiliary variables
		coord Fp, Fmu;
  	coord Tp, Tmu;
  	coord center = {p().x,p().y,p().z};
  	coord Omega = {p().omega.x,p().omega.y,p().omega.z}; // for the moment, it works ONLY in 2D!
  	embed_color_force_RBM  (p, u, mu, csLOCAL, center, Omega, p().r, &Fp, &Fmu);
  	embed_color_torque_RBM (p, u, mu, csLOCAL, center, Omega, p().r, &Tp, &Tmu);
		foreach_dimension()
			p().F.x = Fp.x + Fmu.x;
		/*
		The rest works only in 2D for the moment
		*/
		p().T.x = Tp.x + Tmu.x;
		p().T.y = p().T.x;	
		p().T.z = p().T.x;	
	}
}

void velocity_for_force_free()
{
	foreach_particle(){
		/**
		Translational velocity. */
		coord propulsionAngle = {cosTheta,sinTheta};
		foreach_dimension()
			p().u.x += (p().F.x+p().B.x+p().thrust*propulsionAngle.x)*dt;
		/**
		Angular velocity. */
		foreach_dimension()
			p().omega.x += (p().T.x+p().M.x)*dt*100.; 
		/**
		Here dt and the arbitrary quantity (here 100) are just parameters,
		in a steady case, they do not have a real physical meaning.
		Consider them as relaxation parameters for the timestepper.*/
		/**
		In 2D...must get back to -z variables.*/
		#if dimension == 2
			p().omega.z = p().omega.x;
	}
}
void particle_location_update()
{
	foreach_particle()
		foreach_dimension()
			p().x += p().u.x*dt;
	/**
	Treat periodicity...*/
	// // // TO BE DONE
	/**
	Simple update of angular position.*/
	foreach_particle()
		p().theta.z += p().omega.z*dt;
		
}
