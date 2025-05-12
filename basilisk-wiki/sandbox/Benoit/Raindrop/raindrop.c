/**
# Falling Droplet
This test case is meant to compare centered Navier-Stokes formulations
with the momentum conserving schemes. The same setup is available
in PARIS-SIMULATOR suite.

A water droplet initially at rest is free falling affected by an upstream
air flow,
*/

#define OLDMOMCONS 0
#define NEWMOMCONS 1 // Conserving
#define ADAPTMESH 0
#define THETA 2. // Blending factor for the MOMCONS scheme

//Default refinements
int minlevel = 4;
int MAXLEVEL = 6;

#define GFSOUT 0
#define BVIEWPOST 1


//Liquid
#define RHOL 998.2
#define MUL 8.9e-4
#define SIGMA 0.0728

#define RHOG 1.2
//Air properties
#define MUG 1.98e-5

// Drop and init parameters
#define DIAMETER 3e-3
#define ACCEL 9.81      //On x axis
#define U0 5.0   //On opossite axis to gravity
// We = rho_g*U^2*D/sigma = 1.2*64*0.003/0.0728 = 3.2 --> low Weber, no BU!
// Re = rho_g*U*D/MU = (1.2*8*0.003/1.98)*1e5 = 1460


#define PI 3.141592653589793238
#define MAXTIME 0.1e-1  //Simulation time


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
#include "navier-stokes/perfs.h"

static FILE * fp = NULL;

double maxruntime = HUGE;

/**
Boundary conditions*/

#if !OLDMOMCONS
u.n[top] = neumann(0);
u.n[bottom]  = dirichlet(U0);
u.t[bottom]  = dirichlet(0);
#else
q.n[top] = neumann(0);
q.n[bottom]  = dirichlet(U0*RHOG);
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

#if ADAPTMESH
  init_grid (pow(2.0,minlevel));
#else
  init_grid (pow(2.0,MAXLEVEL));
#endif



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
  sprintf (name, "prop.dat");
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



/*(sqrt(sq(.x[]) + sq(u.y[]) + sq(u.z[])))*/




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

  boundary(all);
#endif

  
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-3,1e-3,1e-3}, MAXLEVEL, minlevel);

}
#endif

/**
Logfile is reporting drop position and velocity every (ten-)timesteps.
Multigrid iterations, mesh size, cpu time and speed (cells/s) are also reported.
*/

event logfile (i += 1,first) {
  double xd = 0., yd = 0., zd = 0., sd = 0., keLiq = 0.,keGas = 0.;
  double vdx = 0., vdy = 0., vdz = 0.;
  

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

  double RefKE = (0.5*U0*U0*RHOL*DIAMETER/2*DIAMETER/2*DIAMETER/2*4.0/3.0*PI);
  keGas = keGas/RefKE;
  keLiq = keLiq/RefKE;



  fprintf(fp,"# 1:t 2:dt 3:sd 4:xd/sd 5:yd/sd 6:zd/sd 7:vdx/sd 8:vdy/sd 9:vdz/sd "
     "10:mgp.i 11:mgu.i 12:grid->tn 13:perf.t 14:perf.speed 15:keGas 16:keLiq\n");
  fprintf (fp,
	   "%.4e %.4e %.4e %.2e %.2e %.2e %.2e %.2e %.2e"
	   " %d %d %ld %.2e %.2e %g %g\n", 
	   t, dt, sd, xd/sd, yd/sd, zd/sd,
	   vdx/sd, vdy/sd, vdz/sd, mgp.i, mgu.i,
	   grid->tn, perf.t, perf.speed,keGas,keLiq);

  fflush (fp);
}


event measurements (t += 0.001)
{

    norm func3 = normf(f);
    double rmsi  = func3.rms;
    double avgi  = func3.avg;
    double volumei  = func3.volume;
    double maxnormi  = func3.max;
   
    fprintf (stderr," %g %g %g %g %g %d %g %g\n", t,rmsi,avgi,volumei,maxnormi, mgp.i,mgp.resa, mgp.resb);
}


event photo (t = 0.0; t += MAXTIME/20; t <= MAXTIME)
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
  sprintf (fname, "raindrop-%4f.ppm", t);
  save (fname); 
}

event animation (i = 1; i += 10; t <= MAXTIME)
{
  static FILE * fp = popen ("ppm2mp4 movie.mp4", "w");

  clear();
  view (fov = 20, camera = "iso", width = 700, height = 700, samples = 1);

  draw_vof ("f");
  squares ("u.y", linear = false, min = 0, max = 15, n = {0,0,1});
  box();
  save (fp = fp);

}




/**
## Results

~~~gnuplot
set xlabel 'Time'
unset logscale y
set key left top
set ylabel 'Y velocity'
set title 'Velocity field as a function of the time'
set autoscale
set xrange [0:0.01]


f(x) = m * x + q       

fit f(x) 'prop.dat' using 1:8 via m, q


mq_value = sprintf("Paramètres :\nm = %f (pente)\nq = %f (valeur à x = 0)", m, q)
set label 1 at 0.006,0.04 mq_value

plot 'prop.dat' u 1:8 w point pointtype 6 ps 0.5 lc rgb "red" title 'vy pour U0 = 5',\
     f(x) lw 3 lc rgb 'black'
     
     
set out


~~~
### Kinetic energy
~~~gnuplot

set key font ",13"
set title font ",13"
set xlabel "Time" font ",13"
set ylabel "Kinetic energy" font ",13"
  
set title 'Kinetic energy as a function of the time'
set logscale y
set key right top

set yrange [0.0000001:0.01]
set xrange [0:0.003]
set logscale y 10
set format y '10^{%L}'
unset format x
plot 'prop.dat' u 1:16 w point pointtype 7 ps 0.5 lc rgb "green" title "Ec"

~~~


*/
