/**
# Rain Droplet blow by stream
This test case is meant to compare centered Navier-Stokes formulations
with the momentum conserving schemes. The same setup is available
in PARIS-SIMULATOR suite.

A water droplet initially at rest is free falling affected by an upstream
air flow,
*/

#define OLDMOMCONS 0
#define NEWMOMCONS 1
#define ADAPTMESH 0
#define THETA 2. // Blending factor for the MOMCONS scheme

//Default refinements
int minlevel = 5;
int MAXLEVEL = 9;

#define GFSOUT 0
#define BVIEWPOST 1


//Liquid
#define RHOL 1000.
#define MUL 8.9e-4
#define SIGMA 0.0728

//Air properties
#define RHOG 1.2
#define MUG 1.98e-5

// Drop and init parameters
#define DIAMETER 3e-3
#define ACCEL 9.81      //On x axis
#define USTREAM 8.046   //On opossite axis to gravity
// We = rho_g*U^2*D/sigma = 1.2*64*0.003/0.0728 = 3.2 --> low Weber, no BU!
// Re = rho_g*U*D/MU = (1.2*8*0.003/1.98)*1e5 = 1460


#define PI 3.14157
#define MAXTIME 1.0e-2  //Simulation time

#include "grid/octree.h"
#if !OLDMOMCONS
#include "navier-stokes/centered.h"
#include "two-phase.h"
#if NEWMOMCONS
#include "navier-stokes/conserving.h"
#endif
#else
#include "momentum.h"
vector u[];
#endif
#include "tension.h"
#include "view.h"

static FILE * fp = NULL;

double maxruntime = HUGE;

/**
Boundary conditions*/

#if !OLDMOMCONS
u.n[top] = neumann(0);
u.n[bottom]  = dirichlet(USTREAM);
u.t[bottom]  = dirichlet(0);
#else
q.n[top] = neumann(0);
q.n[bottom]  = dirichlet(USTREAM*RHOG);
q.t[bottom]  = dirichlet(0);
#endif


p[top]    = dirichlet(0);
p[bottom]   = neumann(0);


/**
The main function can take two optional parameters: the max and min level
for refinement. */

int main (int argc, char * argv[]) {
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);
  if (argc > 2)
    minlevel = atoi (argv[2]);
    
  /**
  We set the domain geometry and initial refinement. */
  L0 = 0.012;
  size (L0);
  origin(-L0/2, -L0/2, -L0/2);
  init_grid (pow(2.0,minlevel));

  DT = 1.25e-5;

  // CFL number
  CFL = 0.4;
  /**
  Physical parameters are set by the macros: 
  densities, viscosities and surface tension. */
  
  f.sigma = SIGMA;

  rho1 = RHOL;
  rho2 = RHOG;
  
  mu1 = MUL;
  mu2 = MUG;

  theta = THETA;
  /**
  Tolerance of the Poisson problem should never be reduced,
  this could lead to non divergence free velocity fields
  and mass conservation errors. */
  
  TOLERANCE = 1e-5; //default is 1e-3.

  char name[80];
  sprintf (name, "E-K2.dat");
  fp = fopen (name, "w");
  run();
  fclose (fp);
}

/**
Initial condition is given by drop position only, there is no flow at the
begining. First step should give a velocity field close to potential flow
solution around the drop, plus small velocities inside the droplet
to avoid infite shear in the interface.*/


void careful_refinement (){

  for (int index = minlevel; index <= MAXLEVEL; index++)
    refine(sq(x)+sq(y)+sq(z) - sq(0.5*DIAMETER + L0/pow(2.0,index)) < 0
            &&
           sq(x)+sq(y)+sq(z) - sq(0.5*DIAMETER - L0/pow(2.0,index)) > 0
            && level < index+1);
}


event init (t = 0) {
  /**
  Domain was reduced in previous versions. This affects dump, so is not done 
  anymore. Check mask-function for this functionality.*/

  if (!restore (file = "dump"))
  #if ADAPTMESH
      careful_refinement();
  #endif
    fraction (f, sq(0.5*DIAMETER) - (sq(x) + sq(y) + sq(z)));
  
  event("properties");
  
}
 
event end (t = MAXTIME) {
  printf ("i = %d t = %g\n", i, t);
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= ACCEL;
}

/**
AMR is applied for velocity tolerance considering using 1% of its RMS value
as cmax. Fraction field is also used for adaptation. */

#if ADAPTMESH
event adapt (i++) {
 
#if OLDMOMCONS
  foreach()
    foreach_dimension()
      u.x[] += q.x[]/rho[];

  boundary({all});
#endif

  double uemax = 1e-1;
  double fmax = 1e-3;
  adapt_wavelet ({f,u}, (double[]){fmax,uemax,uemax,uemax}, MAXLEVEL, minlevel);

}
#endif

/**
Logfile is reporting drop position and velocity every ten-timesteps.
Multigrid iterations, mesh size, cpu time and speed (cells/s) are also reported.
*/

event logfile (i += 10,first) {
  double xd = 0., yd = 0., zd = 0., sd = 0., keLiq = 0.,keGas = 0.;
  double vdx = 0., vdy = 0., vdz = 0.;
  
  if (i == 0){
    fprintf (ferr,
	     "# 1:t 2:dt 3:sd 4:xd/sd 5:yd/sd 6:zd/sd 7:vdx/sd 8:vdy/sd 9:vdz/sd "
	     "10:mgp.i 11:mgu.i 12:grid->tn 13:perf.t 14:perf.speed\n");
  }

  foreach(reduction(+:xd) reduction(+:yd) reduction(+:zd) reduction(+:sd)
	  reduction(+:vdx) reduction(+:vdy) reduction(+:vdz) reduction(+:keGas)) {
    double dv = f[]*dv();
    xd += x*dv;
    yd += y*dv;
    zd += z*dv;
    sd += dv;
    
#if OLDMOMCONS
    vdx += q.x[]/rho[]*dv;    vdy += q.y[]/rho[]*dv;    vdz += q.z[]/rho[]*dv;
    keGas += (sq(q.x[]) + sq(q.y[])+ sq(q.z[]))/rho[]*(1-f[])*dv()/2;
#else
    vdx += u.x[]*dv;    vdy += u.y[]*dv;    vdz += u.z[]*dv;
    keGas += (sq(u.x[]) + sq(u.y[])+ sq(u.z[]))*rho[]*(1-f[])*dv()/2;
#endif
  }

  keLiq = rho1*(sq(vdx)+sq(vdy)+sq(vdz))/2/sd;

  double RefKE = (0.5*USTREAM*USTREAM*RHOL*DIAMETER/2*DIAMETER/2*DIAMETER/2*4.0/3.0*PI);
  keGas = keGas/RefKE;
  keLiq = keLiq/RefKE;

  fprintf (ferr,
	   "%.4e %.4e %.4e %.2e %.2e %.2e %.2e %.2e %.2e"
	   " %d %d %ld %.2e %.2e\n", 
	   t, dt, sd, xd/sd, yd/sd, zd/sd,
	   vdx/sd, vdy/sd, vdz/sd, mgp.i, mgu.i,
	   grid->tn, perf.t, perf.speed);

  fflush (ferr);

  if (t==0)
    fprintf(fp,"# time, KE_gas KE_liq\n");
  
  fprintf (fp, "%g %g %g \n", t, keGas,keLiq) ;
  fflush (fp);
}


/**
Gfsview snapshots are taken at regular time intervals.
Velocity is computed in momentum conserving formulation.
Vorticity field is also saved.*/

#if GFSOUT
event snapshot (t = 0.0; t+=MAXTIME/20; t <= MAXTIME)
{
  if(i > 0)
    dump (file = "dump");
    
  vector omega[];
  #if OLDMOMCONS
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary ((scalar *){u});
  #endif

  foreach(){
    omega.x[] = (u.z[0,1] + u.z[0,-1] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
    omega.y[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1] + u.z[-1])/(2.*Delta);
    omega.z[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
  }

  boundary((scalar *){omega});
  
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  output_gfs (file = name, t = t, list = {f,u,p,omega});
}
#endif

/**
Gfsview snapshots are taken at regular time intervals.
Velocity is computed in momentum conserving formulation.
Vorticity field is also saved.*/

#if BVIEWPOST

#if 1
event snapshot (t = 0.0; t += MAXTIME/10; t <= MAXTIME) {

  /**Set name of output file: */
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[], ff[];
  foreach() {
    pid[] = fmod(pid()*(npe() + 37), npe());
    ff[] = f[] < 1e-4 ? 0 : f[] > 1. - 1e-4 ? 1. : f[];
  }
  boundary ({pid,ff});
  
  dump (list = (scalar *){f,u}, file = name);
}
#endif

event photo (t = 0.0; t += MAXTIME/10; t <= MAXTIME)
{
  vector omega[];
  #if OLDMOMCONS
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary ((scalar *){u});
  #endif

  foreach(){
    omega.x[] = (u.z[0,1] + u.z[0,-1] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
    omega.y[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1] + u.z[-1])/(2.*Delta);
    omega.z[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
  }

  boundary((scalar *){omega});
  
  char fname[80];
  clear();
  view (fov = 20, camera = "front", width = 700, height = 700, samples = 1);
  draw_vof ("f");
  squares ("omega.z", linear = false, spread = 2000, n = {0,0,1});
  sprintf (fname, "vort-%4f.ppm", t);
  save (fname); 
}

event animation (i = 1; i += 20; t <= MAXTIME)
{
  static FILE * fp = popen ("ppm2mp4 movie.mp4", "w");

  clear();
  view (fov = 20, camera = "iso", width = 700, height = 700, samples = 1);

  draw_vof ("f");
  squares ("u.y", linear = false, min = 0, max = 15, n = {0,0,1});
  box();
  save (fp = fp);

}
#endif


/**
## Results

~~~gnuplot Evolution of the drop kinetic energy
set xlabel "time"
set ylabel "dimensionless kinetic energy"
set yrange [*:*]
set xrange [0:0.02]
plot 'reference-E-k2.txt' w l title "Liquid reference",\
'E-K2.dat' u 1:2 title "Gas KE",\
'E-K2.dat' u 1:3 title "Liquid KE"
set term pngcairo
set out "E-k2.png"
replot
~~~

We can check an intermediate state

![Front view at time t = 0.005 s, colored by vorticity [-2000;2000] *1/s*](raindrop/vort-0.005000.ppm)

And the full simulation evolution

![Planed colored by velocity [0;15]*m/s*](raindrop/movie.mp4)

*/