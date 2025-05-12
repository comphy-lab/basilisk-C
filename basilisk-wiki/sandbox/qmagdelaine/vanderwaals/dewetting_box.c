/**
# Dewetting of film in a box

This program illustrates the use of the header file
[vanderwaals.h](/sandbox/qmagdelaine/vanderwaals/vanderwaals.h) to add Van der
Waals-like forces from the subtrate and notably:

* that it seems compatible with *adapt_wavelet()*,
* how you can easily add these forces on the different walls of the numeric
domain.

[vanderwaals.h](/sandbox/qmagdelaine/vanderwaals/vanderwaals.h) and this program
are directly inspired from [dewetting.c](/sandbox/popinet/dewetting.c) and are
an implementation from [Mahady, Afkhami, Kondic, Phys. Fluids 28, 062002
(2016)](http://doi.org/10.1063/1.4949522). 

![Dewetting on the four walls](dewetting_box/video_vof_tracer.mp4)
![Dewetting on the four walls: pressure](dewetting_box/video_pressure.mp4)
*/

/**
We define the geometrical, temporal and resolution parameters: */

#define MIN_LEVEL 4
#define LEVEL 7
#define MAX_LEVEL 7

#define F_ERR 1e-3
#define U_ERR 1e-2

#define T_END 4.
#define DELTA_T 0.04 // for measurements and videos

/**
We solve the Navier-Stokes equation for a two phase flow with surface tension
and Van der Waals forces. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "vanderwaals.h"
#include "view.h"
#define BG 0.6 // background gray level for the videos
#define DG 0. // dark gray for the videos


/**
## Physical parameters */

#define Oh (1./sqrt(10.))
#define hstar 0.01
#define theta_eq (pi/2.)
#define R 0.4

int main() {
  N = 1 << LEVEL;
  f.sigma = 1.;
  mu1 = mu2 = sq(Oh)*f.sigma*2.*R;
  u.t[bottom] = dirichlet(0);
  run();
}

/**
## Initialization

In the initialization of the simulation, we set all the parameters associated to
the Van der Waals forces. The choice of the walls has to be here and not in the
main function, because the list of the walls is set to none in the *defaults*
event, in [vanderwaals.h](/sandbox/qmagdelaine/vanderwaals/vanderwaals.h). */

event init (t = 0) {
  f.theta = theta_eq, f.hs = hstar, f.m = 3., f.k = 2.;
  sprintf (f.walls[0], "bottom");
  sprintf (f.walls[1], "top");
  sprintf (f.walls[2], "left");
  sprintf (f.walls[3], "right");
  fraction (f, max(max(5.*hstar*(1. + 0.01*noise()) - y,
                       5.*hstar*(1. + 0.01*noise()) - L0 + y),
                   max(5.*hstar*(1. + 0.01*noise()) - x,
                       5.*hstar*(1. + 0.01*noise()) - L0 + x)));
}

#if TREE
event adapt (i++) {
  boundary({f, u});
  adapt_wavelet({f, u}, (double[]){F_ERR, U_ERR, U_ERR}, MAX_LEVEL, MIN_LEVEL);
}
#endif

/**
## Post-processings and videos

We juste save the interface and videos. */

event outputs (t = 0.; t += DELTA_T; t <= T_END) {
  static FILE * fp = fopen ("current_time", "w");
  fprintf (fp, "%d %g %g\n", i, t, statsf (u.x).max);
  fflush(fp);
  
  output_facets (f, stderr);
  fflush(stderr);
  
  view (fov = 19., width = 640, height = 640, samples = 1, bg = {BG, BG, BG},
        tx = -0.5, ty = -0.5);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("p", min = -15., max = 5., linear = false, map = cool_warm);
  save ("video_pressure.mp4");

  view (fov = 19., width = 640, height = 640, samples = 1, bg = {BG, BG, BG},
        tx = -0.5, ty = -0.5);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {0.282, 0.459, 0.776}, filled = 1,
             fc = {0.751, 0.785, 0.846});
  draw_vof("f", edges = false, filled = -1, fc = {1., 1., 1.});
  save ("video_vof_tracer.mp4");
}

/**
# Results

~~~gnuplot Evolution of the interface
set size square
raspberry="#FA0F50"
plot 'log' w l lc rgb raspberry
~~~

# References

~~~bib
@article{mahady2016,
  title={A numerical approach for the direct computation of flows including
  fluid-solid interaction: modeling contact angle, film rupture, and dewetting},
  author={Mahady, Afkhami, Kondic},
  journal={Physics of Fluids},
  volume={28},
  issue={6},
  year={2016}
}

*/