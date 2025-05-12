/**
# Implicit Saint-Venant solutions for waves

We solve the Saint-Venant equations either explicitly or implicitly,
in one dimension. */

#include "grid/multigrid1D.h"
#if EXPLICIT
# if ML
#   include "layered/hydro.h"
# else
#   include "saint-venant.h"
# endif
#else
# if ML
#   include "layered/hydro.h"
#   include "layered/implicit.h"
# else
#   include "saint-venant-implicit.h"
# endif
#endif

/**
We vary the wave amplitude to test both linear and nonlinear regimes. */

double amp;

int main()
{
  init_grid (512);
  periodic (right);

  /**
  We set the acceleration of gravity and the maximum timestep (for the
  implicit solver). */
  
  G = 10.;
  DT = 5e-3;
  theta = 1.3;
  
  /**
  We do several runs with increasing wave amplitude. */

  amp = 0.1;
  run();
  amp = 0.4;
  run();
  amp = 1.;
  run();
  system ("cat out-* > log");
}

/**
The initial condition is a simple sinusoidal wave. The water depth is one. */

event init (i = 0) {
  foreach()
    h[] = 1. + amp*sin(2.[-1]*pi*x);
}

/**
We log the evolution of the timestep. */

event logfile (i++) {
  char name[80];
  sprintf (name, "log-%g", amp);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g\n", t, dt, statsf(h).sum);
}

/**
We output wave profiles at a few time intervals. */

event output (t <= 0.2; t += 0.1) {
  char name[80];
  sprintf (name, "out-%g", amp);
  static FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g\n", x, h[]);
  fprintf (fp, "\n");
  fflush (fp);
}

/**
For small wave amplitudes, non-linear effects are small and we get the
following wave evolution for the explicit and implicit solvers.

~~~gnuplot Explicit and implicit solutions for a wave amplitude of 0.1
set xlabel 'x'
set ylabel 'z'
plot 'out-0.1' w l t 'implicit', '../explicit/out-0.1' w l t 'explicit', \
     '../implicit-ml/out-0.1' w l t 'implicit (ML)', \
     '../explicit-ml/out-0.1' w l t 'explicit (ML)'
~~~

The solutions are not too far given that the explicit scheme requires
698 time steps to verify the CFL condition, against only 51 steps for
the implicit scheme.

For an amplitude of 0.4, the explicit solution is close to a weak
shock and requires 840 timesteps against 186 timesteps for the
implicit scheme.

~~~gnuplot Explicit and implicit solutions for a wave amplitude of 0.4
plot 'out-0.4' w l t 'implicit', '../explicit/out-0.4' w l t 'explicit', \
     '../implicit-ml/out-0.4' w l t 'implicit (ML)', \
     '../explicit-ml/out-0.4' w l t 'explicit (ML)'
~~~

For an amplitude of 1., strong shocks develop for both solutions. The
implicit scheme is a bit more oscillatory near the shocks but this
remains acceptable. The explicit solution requires 1092 steps against
471 steps for the implicit scheme.

~~~gnuplot Explicit and implicit solutions for a wave amplitude of 1
plot 'out-1' w l t 'implicit', '../explicit/out-1' w l t 'explicit', \
     '../implicit-ml/out-1' w l t 'implicit (ML)', \
     '../explicit-ml/out-1' w l t 'explicit (ML)'     
~~~
*/
