/**
This is a three-phase flow. Two test cases are embedded in this
sample, the case 1 is not spreading and the case 2 is spreading. */

#define R_VOFLIMIT 1e-3

#include "axi.h"
#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
// #include "conserving-3p.h"
#include "view.h"

/**
Setting the surface tensions for the two cases. */

// case1
#if 1
# define SIGMA_1 0.02
# define SIGMA_2 0.04
# define SIGMA_3 0.03
// case2
#else
# define SIGMA_1 0.04
# define SIGMA_2 0.04
# define SIGMA_3 -0.01
#endif

/**
Minimum and maximum level of refinement are set here, */

const int MAXLEVEL = 8, MINLEVEL = 4;
double volref[3];

int main()
{
  f1.sigma = SIGMA_1;
  f2.sigma = SIGMA_2;
  f3.sigma = SIGMA_3;
  rho1 = 1.0;
  rho2 = 1000.0;
  rho3 = 1000.0;
  mu1 = 1.0e-2;
  mu2 = 1.0e-0;
  mu3 = 1.0e-0;
  run();
}

event init(t = 0)
{
  double r0 = 0.10;
  double x0 = 0.50;

  refine (sq(x - x0) + sq(y) + sq(z) < sq(1.1 * r0) &&
	  sq(x - x0) + sq(y) + sq(z) > sq(0.9 * r0) && level < MAXLEVEL);
  refine ((x < 1.1 * x0 && x > 0.9 * x0 && y > r0) && level < MAXLEVEL);

  foreach () {
    f1[] = f2[] = f3[] = 0.;
    if (sq(x - x0) + sq(y) + sq(z) < sq(r0))
      f3[] = 1.0;
    else if (x < x0)
      f2[] = 1.0;
    f1[] = 1.0 - f2[] - f3[];
  }

  /**
  We calculate the volume integration from the initial condition
  as the reference values. */

  volref[0] = statsf(f1).sum;
  volref[1] = statsf(f2).sum;
  volref[2] = statsf(f3).sum;
}

/**
Refine the grid by the fraction values and velocities */

event adapt(i++)
{
  adapt_wavelet ({f1, f2, f3, u.x, u.y},
		 (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-3}, 
		 maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event logfile (i++)
{
  fprintf (stderr, "%g %g %g %g %g\n", t,
	   statsf(f1).sum - volref[0],
	   statsf(f2).sum - volref[1],
	   statsf(f3).sum - volref[2],
	   normf (u.x).max);
}

event writedata (t += 1.0; t <= 35.)
{
  char name[80];
  sprintf (name, "dump-%g", t);
  dump (name);
}

event movies (i += 10)
{
  clear();
  view (fov = 6.989, quat = {0,0,-0.707,0.707},
	tx = 1e-6, ty = -0.5, width = 780, height = 382);
  draw_vof ("f1", lw = 2);
  draw_vof ("f2", lw = 2);
  draw_vof ("f3", lw = 2);
  mirror ({1}) {
    draw_vof ("f1", lw = 2);
    draw_vof ("f2", lw = 2);
    draw_vof ("f3", lw = 2);
  }
  //  squares ("u.x", spread = -1);
  //  squares ("p");
  cells();
  save ("movie.mp4");
  
  clear();
  view (fov = 0.593863, tx = 0.122977, ty = -0.489743);
  draw_vof ("f1", lw = 2, lc = {1,0,0});
  draw_vof ("f2", lw = 2, lc = {0,1,0});
  draw_vof ("f3", lw = 2, lc = {0,1,1});
  cells();
  save ("zoom.mp4");
}

/**
![Animation of the interfaces and adaptive mesh.](testcase-3p/movie.mp4)

![Zoom on the triple point.](testcase-3p/zoom.mp4)

~~~gnuplot Errors in volume conservation
set xlabel 't'
set ylabel 'Volume error'
plot 'log' u 1:2 w l t 'f1', '' u 1:3 w l t 'f2', '' u 1:4 w l t 'f3'
~~~

~~~gnuplot Maximum horizontal velocity
set xlabel 't'
set ylabel 'u max'
set logscale y
set grid
plot 'log' u 1:5 w l t ''
~~~
*/