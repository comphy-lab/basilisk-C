/**
# Cylinder moving in the x-direction in a fully periodic domain with DLMFD
*/

/** A cylinder twice as heavy as the surrounding fluid is advected in
  the x-direction by a pressure gradient */

# define LEVEL 8 
# include "grid/quadtree.h"
# define DLM_Moving_particle 1
# define adaptive 1
# define NPARTICLES 1
/* # define debugInterior 1 */

# if adaptive
# define MAXLEVEL (LEVEL + 1)
# endif

/** 
# Physical parameters 
*/
# define Uc 1. //caracteristic velocity
# define rhoval 1. // fluid density
# define diam (1.) // particle diameter
# define ReD  200.0 // Reynolds number based on the particle's diameter
# define fs_density_ratio  2. // fluid solid density ratio
# define Ld_ratio 5. // box size-particle diameter ratio
# define Ldomain (Ld_ratio*diam)
# define rhosolid (fs_density_ratio*rhoval) //particle density 
# define tval (rhoval*Uc*diam/ReD) // fluid dynamical viscosity 
 
/** 
# Output and numerical parameters 
*/
# define Tc (diam/Uc) // caracteristic time scale
# define mydt (Tc/400) // maximum time-step (the time step is also adaptive in time but it won't exceed this value)
# define maxtime (4.)
# define tsave (Tc/200.)

/** 
We include the ficitious-domain implementation with a toy-model
granular solver
*/
# include "dlmfd-toygs.h"
# include "view.h"

double deltau;
scalar un[];

int main (int argc, char *argv[]) {
  L0 = Ldomain;

  /* set time step */
  DT = mydt;
      
  /* initialize grid */
  init_grid (1 << (LEVEL));
  
  /* boundary conditions */
  periodic(right); 
  periodic(top);

  /* Convergence criteria */
  TOLERANCE = 1e-3;
  
  run(); 
}

/**
We initialize here the fluid and particle variables. */

event init (i = 0) {
  /* set origin */
  origin (0, 0);

  /* Initialize acceleration (face) vectors for pressure gradient to drive the flow */
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }

  /* If new simulation: set fluid/particles initial conditions */
  if (!restore (file = "dump")) {

    /* Initial condition for the fluid */
    foreach() {
      foreach_dimension() {
	u.x[] = 0;
      }
    }
  
    /* initial condition: particles position */
    particle * p = particles;
    for (int k = 0; k < NPARTICLES; k++) {
	
      GeomParameter gp = {0.};
      gp.center.x = L0 - diam/2. - 3.*L0/((double) pow(2,LEVEL));
      gp.center.y = L0/2;
      gp.radius = diam/2;
      p[k].iscube = 0;
      p[k].iswall = 0;
      p[k].g = gp;

      /* initial condition: particle's velocity */
      coord c = {0., 0., 0.};

      /* Translational velocity U */
      p[k].U = c;
      /* Rotational velocity w */
      p[k].w = c;
    }
  }
  else // Restart of a simulation
    {
      /* the default init  event will take care of it */
    }
}


/* We log the number of iterations of the 
    multigrid solver for pressure and viscosity */
event logfile (i++) {
  deltau = change(u.x,un);
  fprintf (stderr, "%d %g %d %d %g\n", i, t, mgp.i, mgu.i, deltau);
  if (t > maxtime) return 1;
}

event acceleration(i++) {
  face vector av = a;
  coord dp = {20./L0/rhoval, 0, 0};
  foreach_face(){
    av.x[] = dp.x;
  }
}

event output_data (t += tsave; t < maxtime) {
  stats statsvelo;
  view (fov = 22.3366, quat = {0,0,0,1}, tx = -0.465283, ty = -0.439056, bg = {1,1,1}, width = 890, height = 862, samples = 1);
  
  statsvelo = statsf (u.x);
  clear();
  squares ("u.x", map = cool_warm, min = statsvelo.min, max = statsvelo.max); 
  cells();
  save ("movie.mp4");
}

event lastdump (t = maxtime) {
  dump(file = "dump");
}


/** 
# Results
~~~gnuplot particle's x-position
plot "particle-data-0" u 1:2 w l
~~~

~~~gnuplot particle's x-component velocity
plot "particle-data-0" u 1:4 w l
~~~

~~~gnuplot particle's y-position
plot "particle-data-0" u 1:3 w l
~~~

~~~gnuplot particle's y-component velocity
plot "particle-data-0" u 1:5 w l
~~~

~~~gnuplot particle's angular velocity
plot "particle-data-0" u 1:6 w l
~~~
*/

/**
# Animation
fluid's velocity (x-component)
<video width="890" height="862" controls>
<source src="cylinder-periodic-x/movie.mp4" type="video/mp4">
</video>
*/
