/**
# Diffusion from a spherical shell
A diffusion equation in spherical coordinates is solved using
a 1D grid and the ([spherisym.h](/sandbox/ecipriano/src/spherisym.h)) metrics.
*/

#include "grid/multigrid1D.h"
#include "spherisym.h"
#include "run.h"
#include "timestep.h"
#include "diffusion.h"

/**
We set a large total simulation time value in order to ensure reaching
steady state conditions. */

#define tEnd 1000
int maxlevel;
double rho, cp, lambda;

/**
We allocate fields for the diffusive tracer *T* and
the analytical solution *Ta*. */

scalar T[], Ta[];

/**
We set the boundary conditions. */

double TL = 0.8;
double TR = 0.0;

T[left] = dirichlet (TL);
T[right] = dirichlet (TR);

Ta[left] = dirichlet (TL);
Ta[right] = dirichlet (TR);

int main (void) {

  /**
  We set the material properties for the problem. */
  rho = 1.;
  cp = 1.;
  lambda = 1.e-2;

  /**
  The radius of the spherical shell is *X0* while the
  total radius of the environment is *L0*. */

  X0 = 0.2e-3;
  L0 = 1000*X0;
  
  for (maxlevel = 8; maxlevel <= 12; maxlevel++) {
    origin (X0);
    size (L0);
    init_grid (1 << maxlevel);
    DT = 10.;
    run();
  }
}

double dtmax;
double BM, mQ, Q;

event init (i = 0) {

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");

  /**
  We set the initial conditions for the temperature field. */

  foreach()
    T[] = TR;
  boundary({T});

  /**
  We compute the analytical solution. */

  foreach() {
    double T1 = TL; double T2 = TR;
    double r1 = X0; double r2 = L0;
    Ta[] = -r1*r2/(x*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2);
  }
  boundary({Ta});
}

event stability (i++,last) {
  dt = dtnext (DT);
}

event tracer_diffusion (i++) {
  /**
  We compute the diffusivity coefficient. */

  face vector D[];
  foreach_face()
    D.x[] = lambda/rho/cp*fm.x[];
  boundary((scalar *){D});

  scalar theta[];
  foreach()
    theta[] = max(cm[], 1.e-20);
  boundary({theta});

  /**
  We solve the diffusion equations using the 
  Poisson–Helmholtz solver. */

  diffusion (T, dt, D, theta=theta);
}

/**
We write the profile of *T* along the radius every
10 iterations, for any maxlevel of refinement. */

event logprofile (i += 10) {
  char name[80];
  sprintf (name, "outfile-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  double step = L0/(1 << maxlevel);
  for (double x = X0; x < L0; x += step) {
    double error = fabs (interpolate(Ta,x) - interpolate(T,x))/interpolate(Ta,x);
    fprintf (fp, "%g %g %g %.16f\n", x, interpolate(T,x), interpolate(Ta,x), error);
  }
  fprintf (fp, "\n\n");
}

event stop (t = tEnd) {
  return 0;
}

/**
~~~gnuplot  Profile of *T* along the radius for different levels of refinement.
stats 'outfile-8' nooutput

set xr[0:0.02]
set yr[0:0.1]

p "outfile-12" index (STATS_blocks-2) u 1:3 w l t "Analytic", \
  "outfile-8"  index (STATS_blocks-2) w l t "maxlevel 8", \
  "outfile-9"  index (STATS_blocks-2) w l t "maxlevel 9", \
  "outfile-10" index (STATS_blocks-2) w l t "maxlevel 10", \
  "outfile-11" index (STATS_blocks-2) w l t "maxlevel 11", \
  "outfile-12" index (STATS_blocks-2) w l t "maxlevel 12"
~~~
*/
