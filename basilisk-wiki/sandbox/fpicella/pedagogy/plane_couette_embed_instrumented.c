/**
# Plane Couette flow using embedded boundaries
It's a pedagogical case for the resolution of 
2D steady planar [Couette flow](https://en.wikipedia.org/wiki/Couette_flow)

![Here we work with a *plane* Couette flow, but a cylindrical one is quite more common for [experimental demonstrations](https://www.youtube.com/watch?v=57IMufyoCnQ)...](https://i.ytimg.com/vi/57IMufyoCnQ/maxresdefault.jpg){width="400px"}

![Sketch of a 2D plane Couette Flow, from [wikipedia](https://en.wikipedia.org/wiki/Couette_flow)](https://img.favpng.com/6/14/15/taylor-couette-flow-fluid-viscosity-laminar-flow-png-favpng-A10dQVKprLLTxviA66Tmqqqt1.jpg){width="400px"}

### Notes:
Case is shifted so that the walls are at y = 0 and y = H


EPS is used to guarantee the correct functioning of [embed](https://basilisk.fr/src/embed.h) routines.
*/
//#include "grid/multigrid.h"
#include "grid/quadtree.h"
#include "embed.h"
#include "navier-stokes/centered.h"

#define H     1.0
#define UTOP  1.0
#define UBOT  0.0
#define RE    1.0
#define EPS   1e-14
/*
An unique definition of the embedded geometry, that will not change during the simulation.*/
#define CHANNEL_PHI(y,H,EPS) intersection ((y) - (EPS), (H) + (EPS) - (y))


int maxlevel = 7;
int minlevel = 5;

face vector muv[];
scalar un[];

int main()
{
  // Box larger than the channel; the actual fluid is the embedded strip 0 < y < H.
  L0 = 4.0;
  X0 = -L0/2.;
  Y0 = -L0/2.;

  //stokes = true;
  periodic (right);

	N = 1 << minlevel;
  DT = 1e-2;
  TOLERANCE = 1e-8;

  init_grid (N);
  run();
}

event properties (i++)
{
  double nu = UTOP*H/RE;
  foreach_face()
    muv.x[] = nu*fm.x[];
  boundary ((scalar *) {muv});
}

event init (t = 0)
{
  mu = muv;

  // Fluid is the strip 0 < y < H.
  // Outside that strip is solid.
	solid (cs, fs, CHANNEL_PHI(y, H, EPS));

  // No-slip moving walls imposed directly on Cartesian components.
  // Lower wall at y = 0 -> UBOT
  // Upper wall at y = H -> UTOP
  u.x[embed] = dirichlet(y > H/2. ? UTOP : UBOT);
  u.y[embed] = dirichlet(0.);

  // Good practice for embedded Stokes tests
  for (scalar s in {u, p, pf})
    s.third = true;

  foreach() {
    if (cs[] > 0.) {
      double yy = clamp(y, 0., H);
      u.x[] = UBOT + (UTOP - UBOT)*yy/H;
      u.y[] = 0.;
      un[]  = u.x[];
    }
    else
      un[] = 0.;
  }
  boundary ((scalar *) {u, un});

	#if TREE
  /**
  Optional but useful: pre-refine around the embedded boundaries before
  time-stepping starts.
  */
  astats ss;
  int it = 20;
  do {
    ss = adapt_wavelet ((scalar *){cs}, (double[]){1e-12}, maxlevel, minlevel);
		refine (level < maxlevel && cs[] >= 1.); // so to ensure maximum refinement within fluid zone.
		solid (cs, fs, CHANNEL_PHI(y, H, EPS));
  } while ((ss.nf || ss.nc) && --it);
#endif
}


event adapt (i++)
{
  adapt_wavelet ((scalar *){cs,u},
                 (double[]){1e-3, 1e-3, 1e-3},
                 maxlevel, minlevel);

  refine (level < maxlevel && cs[] >= 1.);
	solid (cs, fs, CHANNEL_PHI(y, H, EPS));
}

/*
A complete logfile event + gnuplot session, used to extract data from the computation.*/

#include <stdio.h>
#include <math.h>

/**
## Results

~~~gnuplot Couette convergence. Getting as low as 1e-10, not bad!
set grid
set xlabel "iteration"
set ylabel "du"
set title "Convergence history"
set logscale y

plot 'du.dat' using 1:3 \
     with lines lw 4 lc rgb "blue" title "du"
~~~

~~~gnuplot Couette velocity profile. Beautifully linear.
set grid
set xlabel "u_x"
set ylabel "y"
set title "Velocity profile at convergence"
unset logscale y

plot 'profile_analytic.dat' using 3:2 \
     with lines lw 5 lc rgb "red" title "analytic",\
     'profile_interp.dat' using 3:2 \
     with points pt 5 ps 1.2 lc rgb "green" title "interpolated",\
     'profile_cells.dat' using 3:2 \
     with points pt 7 ps 1.8 lc rgb "blue" title "cell-centered"
~~~

~~~gnuplot Couette velocity error. Please note how cell-centered (non-interpolated!) datapoints are as accurate as 1e-12. On the other hand, using interpolate function to get values where they have not been computed (here $y \in [0,0.02]$) can lead to quite large artifacts.
set grid
set xlabel "|u_x-u_{exact}|"
set ylabel "y"
set title "Error with respect to analytic Couette profile"
set logscale x
set logscale y

plot 'profile_interp.dat' using 6:2 \
     with points pt 5 ps 1.2 lc rgb "green" title "interpolated error",\
		 'profile_cells.dat' using 6:2 \
     with points pt 7 ps 1.8 lc rgb "blue" title "cell-centered error"
~~~
*/

event logfile (i++; i <= 4000)
{
  double du = change (u.x, un);
  fprintf (stderr, "i=%05d t=%+6.5e du=%+6.5e\n", i, t, du);

  static FILE * fdu = fopen ("du.dat", "w");
  fprintf (fdu, "%d %+6.5e %+6.5e\n", i, t, du);
  fflush (fdu);

  if (i > 0 && du < 1e-10) {

    /*
    Cell-centered values near x = 0
    columns:
    1: x
    2: y
    3: u.x
    4: u.y
    5: u_exact
    6: |u.x - u_exact|
    */
    FILE * fp1 = fopen ("profile_cells.dat", "w");
    foreach (serial)
      if (cs[] > 0. && x >= 0 && x <= Delta/2.) { // so to take cells that are on x+ side. 
        double uexact = UBOT + (UTOP - UBOT)*y/H;
        double err = fabs(u.x[] - uexact);
        fprintf (fp1, "%+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
                 x, y, u.x[], u.y[], uexact, err);
      }
    fclose (fp1);

    /*
    Interpolated values at x = 0, y in [0,H]
    columns:
    1: x
    2: y
    3: u.x
    4: u.y
    5: u_exact
    6: |u.x - u_exact|
    */
    FILE * fp2 = fopen ("profile_interp.dat", "w");
    int np = 200;
    for (int j = 0; j <= np; j++) {
      double xx = 0.;
      double yy = H*j/(double)np;
      double ux = interpolate (u.x, xx, yy);
      double uy = interpolate (u.y, xx, yy);
      double uexact = UBOT + (UTOP - UBOT)*yy/H;
      double err = fabs(ux - uexact);

      fprintf (fp2, "%+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
               xx, yy, ux, uy, uexact, err);
    }
    fclose (fp2);

    /*
    Analytic profile at x = 0
    columns:
    1: x
    2: y
    3: u_exact
    */
    FILE * fp3 = fopen ("profile_analytic.dat", "w");
    for (int j = 0; j <= np; j++) {
      double xx = 0.;
      double yy = H*j/(double)np;
      double uexact = UBOT + (UTOP - UBOT)*yy/H;
      fprintf (fp3, "%+6.5e %+6.5e %+6.5e\n", xx, yy, uexact);
    }
    fclose (fp3);

    return 1;
  }
}