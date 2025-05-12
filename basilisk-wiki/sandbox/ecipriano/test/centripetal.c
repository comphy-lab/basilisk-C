/**
# Suspending force - Volumetric Formulation

Test for the articifial suspending force applied to the
whole liquid volume. The test is perfomed in zero-gravity
conditions and just the centripetal force contribution
is considered. The goal is to observe the influence of
this force on the liquid internal recirculationa and the
oscillations in the velocity field.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "centripetal.h"
#include "view.h"

/**
We define the maximum level of refinement, the total
simulation time, the initial radius of the droplet and its
center coordinates. Other variables are for post-processing
purposes.
*/

int maxlevel = 6;
double tend = 0.3;
double R0 = 0.5e-3;
double XC, YC;
double trmin, trmax;

/**
We initialize a tracer field, which is used to observe
the liquid recirculation.
*/

scalar tr[];

int main (void) {
  /**
  Material properties used for this simulations. */

  rho1 = 626.7, rho2 = 17.51;
  mu1 = 1.e-4, mu2 = 1.e-5;

  /**
  The tracer is advected using the vof-fluxes. */

  f.tracers = {tr};

  /**
  We set a maximum time-step and the problem geometry. */

  DT = 1.e-3;
  L0 = 5.*R0;
  XC = 0.5*L0, YC = 0.;
  N = 1 << maxlevel;
  run ();
}

#define circle(x,y,R)(sq(R) - sq(x - XC) - sq(y - YC))

event init (i = 0) {
  /**
  A droplet is initialized on the bottom wall of the domain,
  and the tracer is inizialized in order to be zero in the
  gas phase and equal to the *y* coordinate inside the liquid.
  */

  fraction (f, circle (x,y,R0));
  foreach()
    tr[] = y*f[];

  /**
  We find the max and min value of tr in order to impose
  *max* and *min* in the *squares()* function of bview.
  */

  trmin = 1000., trmax = 0.;
  foreach (reduction(min:trmin), reduction(max:trmax)) {
    if (f[] > 1.e-10) {
      trmin = min (trmin, tr[]);
      trmax = max (trmax, tr[]);
    }
  }

  /**
  The default parameters of the centripetal force are
  gathered in the structure *sfm*. These parameters are
  overwritten as follows. The parameter *sigma* is a fake
  surface tension which is used just to compute a realistic
  stability condition.
  */

#ifdef CENTRIPETAL
  sfm.p = (coord){XC,YC};
  sfm.eps = 1.e-5;
  sfm.sigma = 0.03;
#endif
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, (double[]){0,1.e-3,1.e-3}, maxlevel);
}
#endif

/**
We write the maximum value of the velocity field. */

event logfile (t += 0.01) {
  fprintf (stderr, "%f %.10e\n", t, max (statsf(u.x).max, statsf(u.y).max));
}

/**
We write a video with the evolution of the tracer. */

event movie (t += 0.001; t <= tend) {
  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("tr", min = trmin, max = trmax);
  draw_vof ("f", lw = 1.5);
  save ("movie.mp4");
}

/**
## Results

The centripetal force applied to the whole liquid volume
promotes liquid internal recurculation causing the motion
of the tracer.

![Evolution of the tracer](centripetal/movie.mp4)(width="800" height="600")

While this effect can be negligible in some conditions,
when the external convective fluxes lead to higher
liquid recirculation effects, it can be resolved using
the interfacial version of the suspending force: [icentripetal.c](icentripetal.c).
That method allows a lower maximum velocity to be obtained.

~~~gnuplot Evolution of the maximum value of the velocity field
set xlabel 'time [s]'
set ylabel 'max velocity'
plot "log" u 1:2 w l
~~~

*/

