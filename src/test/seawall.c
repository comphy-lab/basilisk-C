/**
# Solitary wave overtopping a seawall

This test case seeks to reproduce the experimental data of [Hsiao and
Lin, 2010](/src/references.bib#hsiao2010) for a solitary wave
overtopping a seawall. It was also studied by [Lannes and Marche,
2014](/src/references.bib#lannes2014).  Both nonlinear and dispersive
effects are important. To illustrate this, we model the wave both with
non-dispersive [Saint-Venant equations](/src/saint-venant.h) and with
the dispersive [Green-Naghdi equations](/src/green-naghdi.h) or
[layered model](/src/layered/hydro.h) in one dimension. */

#include "grid/multigrid1D.h"
#if SAINT_VENANT
# include "saint-venant.h"
#elif ML
# include "layered/hydro.h"
# include "layered/nh.h"
scalar h;
vector u;
#else
# include "green-naghdi.h"
#endif

/**
The domain is 15 metres long. The acceleration of gravity is set to
9.81 m/s^2^. The problem is discretised with 2048 grid
points. */

int main()
{
  L0 = 15.;
  G = 9.81;
  N = 1 << 11;
#if ML
  breaking = 0.15;
  TOLERANCE = 1e-4 [*];
  CFL_H = 0.5;
#endif
  run();
}

/**
The initial conditions are for a [Green-Naghdi soliton](soliton.c) in
a water depth of 0.2 metres and a relative soliton amplitude of 0.35
(type 1 wave of [Hsiao and Lin,
2010](/src/references.bib#hsiao2010)). */

double h0 = 0.2, A = 0.35;

double sech2 (double x) {
  double a = 2./(exp(x) + exp(-x));
  return a*a;
}

double soliton (double x, double t)
{
  double c = sqrt(G*(1. + A)*h0), psi = x - c*t;
  double k = sqrt(3.*A*h0)/(2.*h0*sqrt(h0*(1. + A)));
  return A*h0*sech2 (k*psi);
}

event init (i = 0)
{
  double c = sqrt(G*(1. + A)*h0);
  foreach() {
    double eta = soliton (x - 5.9, t);
    
    /**
    Here we define the bathymetry (Fig. 1 of Hsiao and Lin). We need
    to shift the origin by 3 metres to match the shift in origin
    introduced by Hsiao and Lin when reporting computation results. */

    x += 3.;
    zb[] = (x < 10. ? 0. :
	    x > 13.6 && x <= 13.9 ? 3.6/20. + (x - 13.6)*0.076/(13.9 - 13.6) :
	    x > 13.9 && x <= 13.948 ? 3.6/20. + 0.076 :
	    x > 13.948 && x <= 14.045 ? 
	    3.6/20. + 0.076 - (x - 13.948)*(0.076 - 0.022)/(14.045 - 13.948)
	    : (x - 10.)/20.)
      - h0;
    h[] = max (0., eta - zb[]);
    u.x[] = c*eta/(h0 + eta);
  }
}

/**
We output timeseries of free surface elevation for the first 12
seconds for comparison with experimental records. */

event timeseries (i++; t <= 12) {

  /**
  We also add a quadratic bottom friction term of the form
  $$
  \partial_tu_x = - 10^{-2}|\mathbf{u}|u_x/h
  $$ */
  
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-2*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
  
  /**
  The coordinates of the various gauges are given in the legend of
  Figure 8 of Hsiao and Lin. */

  static FILE * fp0 = fopen ("g0", "w");
  fprintf (fp0, "%g %g\n", t, interpolate (eta, 5.9, 0));
  fprintf (stderr, "%g %g\n", t, interpolate (eta, 7.6, 0));
  static FILE * fp10 = fopen ("g10", "w");
  fprintf (fp10, "%g %g\n", t, interpolate (eta, 9.644, 0));
  static FILE * fp22 = fopen ("g22", "w");
  fprintf (fp22, "%g %g\n", t, interpolate (eta, 10.462, 0));
  static FILE * fp28 = fopen ("g28", "w");
  
  /**
  These gauges are on "dry" land and seem to use as reference the local
  height of the topography (rather than sealevel), presumably so that
  they measure "inundation depth" rather than wave elevation. */

  fprintf (fp28, "%g %g\n", t, (interpolate (eta, 10.732, 0)
				- interpolate (zb, 10.732, 0)));
  static FILE * fp37 = fopen ("g37", "w");
  fprintf (fp37, "%g %g\n", t, (interpolate (eta, 11.005, 0)
				- interpolate (zb, 11.005, 0)));
  static FILE * fp40 = fopen ("g40", "w");
  fprintf (fp40, "%g %g\n", t, (interpolate (eta, 11.12, 0)
				- interpolate (zb, 11.12, 0)));
}

/**
~~~gnuplot Comparison between experimental (dots) and numerical timeseries
set term svg enhanced font ",10" size 640,800
set multiplot layout 7,1
set tmargin 1.
set bmargin 1.
set yrange [-0.02:0.1]
set ytics 0,0.05,0.1
set grid
set pointsize 0.5

plot [-2:2]'../wg1' pt 7 t 'WG1', \
           '../seawallsv/g0' w l lt 4 lw 1 t '', \
           '../seawall-ml/g0' w l lt 6 lw 1 t '', \
           'g0' w l lt 3 lw 2 t ''
plot [0:10]'../wg3' pt 7 t 'WG3', \
           '../seawallsv/log' w l lt 4 lw 1 t '', \
           '../seawall-ml/log' w l lt 6 lw 1 t '', \
           'log' w l lt 3 lw 2 t ''
plot [2:12]'../wg10' pt 7 t 'WG10', \
           '../seawallsv/g10' w l lt 4 lw 1 t '', \
           '../seawall-ml/g10' w l lt 6 lw 1 t '', \
           'g10' w l lt 3 lw 2 t ''
plot [2:12]'../wg22' pt 7 t 'WG22', \
           '../seawallsv/g22' w l lt 4 lw 1 t '', \
           '../seawall-ml/g22' w l lt 6 lw 1 t '', \
           'g22' w l lt 3 lw 2 t ''
plot [2:12]'../wg28' pt 7 t 'WG28', \
           '../seawallsv/g28' w l lt 4 lw 1 t '', \
           '../seawall-ml/g28' w l lt 6 lw 1 t '', \
           'g28' w l lt 3 lw 2 t ''
plot [2:12]'../wg37' pt 7 t 'WG37', \
           '../seawallsv/g37' w l lt 4 lw 1 t '', \
           '../seawall-ml/g37' w l lt 6 lw 1 t '', \
           'g37' w l lt 3 lw 2 t ''
plot [2:12]'../wg40' pt 7 t 'WG40', \
           '../seawallsv/g40' w l lt 4 lw 1 t '', \
           '../seawall-ml/g40' w l lt 6 lw 1 t '', \
           'g40' w l lt 3 lw 2 t ''
unset multiplot
~~~

The agreement is very satisfactory for this difficult problem and can
be compared with the 2D Navier-Stokes VOF simulations of Hsiao and Lin
(Figure 8 left column), as well as with the results of [Lannes and
Marche, 2014](/src/references.bib#lannes2014) (Figure 10). Dispersive
effects are particularly clear for gauge 3 after 6 seconds and are
well reproduced by the model.

The results for the Saint-Venant solver are the thin magenta lines and
the results for the layered model are the thin dark blue lines.

We also output wave profiles at times corresponding to those of Figure
2 of Hsiao and Lin. */

event profile (t = {2.63,2.89,3.01,3.19,3.35,3.71}) {
  char name[20]; sprintf (name, "p-%g", t);
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g %g\n", x, h[], u.x[], zb[]);
  fclose (fp);
}

/**
~~~gnuplot Green-Naghdi (left column) versus Saint-Venant (right column)
reset

set term svg enhanced font ",10" size 640,800
set multiplot layout 6,2
set rmargin 2.
set lmargin 4.
set tmargin 0.
set bmargin 0.
set yrange [-0.1:0.2]
set ytics -0.1,0.1,0.2
set size ratio -1
unset key

plot [9.6:10.6]'p-2.63' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'../seawallsv/p-2.63' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'p-2.89' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'../seawallsv/p-2.89' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.01' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawallsv/p-3.01' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.19' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawallsv/p-3.19' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.35' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawallsv/p-3.35' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.71' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawallsv/p-3.71' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1

unset multiplot
~~~

The results are displayed both for the Green-Naghdi solver (left
column) and for the Saint-Venant solver (right column).  The
Green-Naghdi results are again remarkably similar to the experimental
snapshots of Hsiao and Lin (Figure 2) despite the complexity of the wave
breaking process. 

~~~gnuplot Green-Naghdi (left column) versus layered (right column)
reset

set term svg enhanced font ",10" size 640,800
set multiplot layout 6,2
set rmargin 2.
set lmargin 4.
set tmargin 0.
set bmargin 0.
set yrange [-0.1:0.2]
set ytics -0.1,0.1,0.2
set size ratio -1
unset key

plot [9.6:10.6]'p-2.63' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'../seawall-ml/p-2.63' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'p-2.89' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [9.6:10.6]'../seawall-ml/p-2.89' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.01' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawall-ml/p-3.01' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.19' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawall-ml/p-3.19' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.35' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawall-ml/p-3.35' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'p-3.71' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1
plot [10.2:11.2]'../seawall-ml/p-3.71' u 1:($2+$4) w l lc -1 lw 2, \
     '' u 1:4 w filledcu x1 lc -1

unset multiplot
~~~
*/
