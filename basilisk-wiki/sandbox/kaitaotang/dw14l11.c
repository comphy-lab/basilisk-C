/**
# Droplet breakup in crossflow

We wish to study the behaviour of a single droplet breaking up in a
wind, far from any boundaries.

We use the centered Navier--Stokes solver and log performance
statistics. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

/**
We have two phases e.g. air and water. For large viscosity and density
ratios, the harmonic mean for the viscosity tends to work better than
the default arithmetic mean. We "overload" the default by defining the
*mu()* macro before including the code for two phases. */

// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"

/**
We also need surface tension, and in 3D only we will use the
$\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) to display the vortices using
Basilisk View. */

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "navier-stokes/conserving.h"
#include "view.h"

/**
The tag.h library is included to allow for retrieval of fragment statistics. */

#include "tag.h"
#include <sys/stat.h>
#include <sys/types.h>

/**
We can control the maximum runtime. */

#include "maxruntime.h"


/**
We utilise the Manifold Death algorithm of [Chirco et al.](https://doi.org/10.1016/j.jcp.2022.111468) 
to artificially perforate thin bag films with square holes at late times. This allows for establishing grid
convergence of fragment statistics, see discussions by [Tang et al.](https://doi.org/10.1017/jfm.2023.605)*/
#define SQUARES
#include "signature.h"

/**
The density ratio is 833 and the dynamic viscosity ratio is 55, corresponding to air-water systems. */

#define RHOR 833.
#define MUR 55.


/**
The Weber number is 15 and the Ohnesorge number is 0.001, which fall within the inviscid bag breakup regime. */
# define WE 15.
# define OH 0.001
# define MAXTIME 150.


/**
The domain width is set 15 times larger than the droplet diameter *2R0*, 
and the initial centre-of-mass position *X0* is set close to the left boundary. 
This allows for enough room for the subsequent downstream transport of the
drop as it flattens and breaks up. U is the speed of the ambient flow which
is set as unity.*/
#define WIDTH 30.
#define X0 5.
#define R0 1
#define U 1.

// Set the eccentricity e0 within (0,1) for an oblate drop
#define e0 0.5 


/**
Here we specify the maximum grid resolution level $L_{max}$ (MAXLEV) and the 
minumum level $L_{min}$ for the Basilisk Adaptive Mesh Refinement (AMR) Scheme.
The signature level $L_{sig}$ is set as 14 for the MD algorithm.*/
int MAXLEV = 15, MINLEV = 5, SIGN_LEV = 14;


/**
File pointers for recording the kinematic and energetic properties of the droplet.*/
FILE * fk = NULL;
FILE * fe = NULL;

/** Here sb is the droplet volume. xb, yb, zb, vbx, vby, vbz are the components
of the droplet centre-of-mass location and velocity.*/
double xb = 0., yb = 0., zb = 0., sb = 0.;
double vbx = 0., vby = 0., vbz = 0.;

// Liquid- and gas-phase dissipation rates
double diss_w = 0., diss_a = 0;


/** We set Dirichlet boundary conditions for the incoming velocity at the left 
boundary and pressure at the right boundary. The latter corresponds to a pressure
outflow condition. */
u.n[left] = dirichlet(U);
u.t[left] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = neumann(0.);

p[left] = neumann(0.);
p[right] = dirichlet(0.);

// Scalar fields required by the MD algorithm. 
scalar phii[], sign[], M[];

// Scalar fields for output.
scalar l2[], omegay[], pp[];

/** MD parameters. *max_change* is the maximum number of holes allowed to form per
call of the algorithm, *1/prob* is the probability an eligible thin film cell at $L_{sig}$ 
is perforated, and *large* enables the detection of thin films on a local 5*5*5 stencil 
(reduced to 3*3*3 stencils if set as false). */
const int max_change = 2000; 
const int prob = 17500;
const double T_MD = 0.1;
bool large = true;



/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    MAXLEV = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement. */
  size (WIDTH);
  origin (-X0, -L0/2, -L0/2.);
  init_grid (1 << MINLEV);

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */
  
  rho1 = 1.;
  rho2 = 1./RHOR;
  f.sigma = rho2*sq(U)*(2.*R0)/WE;
  mu1 = OH*sqrt(rho1*f.sigma*2.*R0);
  mu2 = mu1/MUR;

  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  TOLERANCE = 1e-4;

  char name[80];
  sprintf (name, "3Denergy-%d.dat", MAXLEV);
  fe = fopen (name, "w");

  sprintf (name, "3Dkinetics-%d.dat", MAXLEV);
  fk = fopen (name, "w");

  /**
  We create folders for storing MD debug information, drop interface 
  profiles and fragment statistics. */
  mkdir("./chirco_debugs", 0777);
  mkdir("./interface_profiles", 0777);
  mkdir("./frag_stats", 0777);

  run();

  fclose (fe); fclose (fk);
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a spherical shell encompassing the drop. We
then initialise the volume fraction for a bubble initially at (0,0,0)
of diameter unity. In the fraction command, the signs are chosen to
ensure that a droplet, not a bubble is formed (i.e. f=1 inside, f=0
outside).*/

event init (t = 0) {
  if (!restore (file = "restart")) {
    double a0 = R0 * pow(1.-sq(e0), -1./6.);
    double b0 = a0 * pow(1.-sq(e0), 1./2.);
    
    refine (sq(x) + sq(y) + sq(z) - sq(1.1*a0) < 0 && sq(x) + sq(y) + sq(z) - sq(0.9*b0) > 0 && level < MAXLEV);
    fraction (f, -pow(sq(x/b0) + sq(y/a0) + sq(z/a0), 0.5) + 1.);
  }
  else fprintf (stderr, "Dumpfile restored!\n");
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event adapt (i++) {
  double femax = 1e-2, uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){femax, uemax, uemax, uemax}, MAXLEV, MINLEV);
}


/**
This function calculates the liquid- and gas-phase instantaneous 
dissipation rates based on the velocity vector $u$. */

int dissipation_rate (vector u, double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1] - u.x[-1])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    
    rateWater += mu1*f[]*sqterm; //water
    rateAir   += mu2*(1. - f[])*sqterm; //air
  }

  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}


/**
## Outputs

Every 0.05 simulation time unit, we output the time, volume, position, and
velocity of the bubble. */

event logfile (t = 0; t <= MAXTIME; t += 0.05) {
  xb = 0., yb = 0., zb = 0., sb = 0.;
  vbx = 0., vby = 0., vbz = 0.;
  double ke_d = 0., ke_a = 0., a = 0.;

  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb) 
      reduction(+:vbx) reduction(+:vby) reduction(+:vbz) 
      reduction(+:sb) reduction(+:ke_d) reduction(+:ke_a) reduction(+:a)) {
    double dv_l = f[]*dv();
    xb += x*dv_l;
    yb += y*dv_l;
    zb += z*dv_l;
    vbx += u.x[]*dv_l;
    vby += u.y[]*dv_l;
    vbz += u.z[]*dv_l;
    sb += dv_l;

    // Kinetic energy of the droplet and ambient airflow
    ke_d += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho1*f[];
    ke_a += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho2*(1-f[]);
  }

  // Droplet surface area
  a = interface_area(f);


  // Viscous dissipation rate
  double rates[2];
  dissipation_rate(u, rates);
  diss_w = rates[0], diss_a = rates[1];
  
  // Droplet CM's location and velocity
  xb = xb/sb; yb = yb/sb; zb = zb/sb; vbx = vbx/sb; vby = vby/sb; vbz = vbz/sb;

  // Droplet's streamwise length and spanwise width
  double h_max = -100.0, h_min = 100.0, t_max = -100.0, t_min = 100.0;
  foreach(reduction(max:h_max) reduction(min:h_min) reduction(max:t_max) reduction(min:t_min))
    if (f[] > 1e-6 && f[] < 1.-1e-6) {
      h_max = fmax(h_max, y); h_min = fmin(h_min, y); 
      t_max = fmax(t_max, x); t_min = fmin(t_min, x);
    }

  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb,
	   vbx, vby, vbz, h_max-h_min, t_max-t_min);
  fflush (stderr);

  fprintf (fk,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb, vbx, vby, vbz, h_max-h_min, t_max-t_min);
  fflush (fk);

  fprintf (fe, 
          "%.8f %.8f %.8f %.8f %.8f %.8f\n", 
          t, ke_d, ke_a, diss_w, diss_a, (a-(4*M_PI*sq(R0)))*f.sigma);
  fflush (fe);
}

/**
Every 0.02 time unit, we output a full snapshot of the simulation, to enable
restart and for visualisation. In three dimensions, we compute the value of 
the $\lambda_2$ field which will be used for visualisation of vortices, as 
well as the streamwise vorticity $\omega_y = \partial_x u_z - \partial_z u_x$ 
and pressure field pp[]. Note that we use a separate scalar field pp[] to copy 
information from the original pressure field p[] due to a bug that prevents 
proper restart of the program if this is not implemented.*/

event snapshot (t = 0; t <= MAXTIME; t += 0.02) {  
#if dimension == 3
  lambda2 (u, l2);
  foreach() {
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
    pp[] = p[];
  }
  boundary ({omegay, pp});
#endif
  
  char name[80];
  sprintf (name, "restart");
  dump (file = name);
}


/**
Every 0.05 time unit, we output ASCII files recording the droplet interface
profile. As each parallel process generates a portion of the interface, these
have to be combined using the cat command in post-processing. */

event outputInterface(t = 0; t <= MAXTIME; t += 0.05) {
  char resultname[100];
  sprintf(resultname, "./interface_profiles/results%4.2f_%d.dat", t, pid());
  FILE * fp = fopen(resultname, "w");

  scalar xpos[], ypos[], zpos[];
  position (f, xpos, {1, 0, 0});
  position (f, ypos, {0, 1, 0});
  position (f, zpos, {0, 0, 1});
  foreach()
      if (xpos[] != nodata && fabs(zpos[]) <= 0.5*WIDTH/(1 << MAXLEV))
        fprintf (fp, "%g %g %g\n", xpos[], ypos[], zpos[]);
  
  fclose (fp);
  
  /*
  char command[80];
  sprintf(command, "LC_ALL=C  cat results* > interface_%g.dat", t);
  system(command);
  sprintf(command, "LC_ALL=C  rm results*");
  system(command);
  */
}


// This event computes the size and velocity of all fragments using tag.h.
event count_droplets (t = 90; t <= MAXTIME; t += 0.02) {
  scalar m_cd[];
  foreach() m_cd[] = f[] > 0.01;
  int n = tag (m_cd);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach (serial)* to avoid doing a parallel traversal when
  using OpenMP. This is because we don't have reduction operations for
  the *v* and *b* arrays (yet). */

  double v[n], ux_drop[n], uy_drop[n], uz_drop[n];
  coord b[n];

  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = ux_drop[j] = uy_drop[j] = uz_drop[j] = 0.;

  foreach_leaf()
    if (m_cd[] > 0) {
      int j = m_cd[] - 1;
      v[j] += dv()*f[];

      coord p = {x,y,z};
      foreach_dimension()
	    b[j].x += dv()*f[]*p.x;

      ux_drop[j] += dv()*f[]*u.x[];
	  uy_drop[j] += dv()*f[]*u.y[];
	  uz_drop[j] += dv()*f[]*u.z[];
    }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ux_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

  static FILE * fdrop = fopen("droplets.dat", "w");

  if (n > 1 && pid() == 0) 
    for (int j=0; j<n; j++) {
      fprintf (fdrop, "%d %g %d %g %g %g %g %g %g %g\n", i, t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], ux_drop[j]/v[j], uy_drop[j]/v[j], uz_drop[j]/v[j]);
    }
    fflush(fdrop);
    
    char resultname[100];
    sprintf(resultname, "./frag_stats/frag_stat_%4.2f.dat", t);
    FILE * ffrag = fopen(resultname, "w");
  
    for (int j=0; j<n; j++) {
      fprintf (ffrag, "%g %d %g %g %g %g %g %g %g\n", t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], ux_drop[j]/v[j], uy_drop[j]/v[j], uz_drop[j]/v[j]);
    }
    
    fclose(ffrag);
}


/**
This event calls the MD algorithm every T_MD simulation time unit to
detect and artificially perforate bag films that are sufficiently thin. */
event neck_detect(t = 0; t <= MAXTIME; t += T_MD){
  foreach(){
    phii[] = 2*f[] - 1;
    sign[] = 7;
  }  
  
  for (int ilev = depth() - 1; ilev >= SIGN_LEV; ilev--)  
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  
  // Computing the local geometric signature field *sign*.
  compute_signature_neigh_level (f, phii, sign, SIGN_LEV);

  if (pid()==0)  
    printf("time %g level used for moments %d and depth is %d \n", t, SIGN_LEV, depth()); 
  
  for (int ilev = SIGN_LEV; ilev < depth(); ilev++)  
    foreach_level(ilev){
      sign.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        sign.prolongation (point, sign);
        phii.prolongation (point, phii);
      }
    }

  /** This routine scans all eligible cells at $L_{sig}$ and randomly perforate
  these cells with a probability of 1/prob.*/
  change_topology (f, sign, M, SIGN_LEV, max_change, large, prob); 
}



/**
This event uses bview to produce a movie containing the front view of the drop,
allowing for detailed observation of the film perforation process. */
#define POPEN(name, mode) fopen (name ".ppm", mode)

event movie (t += 0.02)
{
  view (quat = {-0.500, -0.500, -0.500, -0.500}, fov = 15, near = 0.01, far = 1000, tx = 0.000, ty = -0.001, tz = -0.888, width = 720, height = 720, bg = {1,1,1}, samples = 4);
  clear();
  draw_vof ("f");
  box();

  static FILE * fp = POPEN ("movie", "w");
  save (fp = fp);
}