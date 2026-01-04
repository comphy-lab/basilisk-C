/**
# Typical (1D) tsunami wave

The initial surface displacement is a Gaussian bump with an amplitude
of 10 metres and caracteristic width 25 km. The water depth is
constant at 7000 metres. These parameters are fairly typical of
fault-triggered tsunamis close to subduction zones.

As illustrated below, dispersive effects are important.

Note also that older versions of the non-hydrostatic multilayer solver
were affected by a grid-scale instability for this test case.

~~~gnuplot Wave shape after one hour
set xlabel 'x (km)'
set ylabel 'Amplitude (m)'
unset key
plot 'log' u ($1/1e3):2 w l
~~~

~~~gnuplot Horizontal velocity after one hour
set ylabel 'Horizontal velocity (m/s)'
plot 'log' u ($1/1e3):3 w l
~~~

~~~gnuplot Vertical velocity after one hour
set ylabel 'Vertical velocity (m/s)'
plot 'log' u ($1/1e3):4 w l
~~~
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/check_eta.h"

int main()
{
  size (2500e3);
  origin (-L0/2.);
  G = 9.81;
  N = 2 << 12;
  CFL_H = 0.5;
  run();
}

event init (i = 0)
{
  res_eta = new scalar;
  conserve_elevation();
  foreach() {
    zb[] = - 7000;
    h[] = - zb[] + 10.*exp(-sq(x/25e3));
  }
}

#if 1
event gnuplot (i += 10)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n"
	     "set xlabel 'x (km)'\n"
	     "set ylabel 'Wave height (m)'\n"
	     //	     "set yrange [-4:10]\n"
	     "set xrange [%g:%g]\n"
	     // "unset key\n"
	     , X0/1e3, (X0 + L0)/1e3);
  FILE * fp1 = fopen ("gnuplot", "w");
  foreach()
    fprintf (fp1, "%g %g %g %g\n", x, eta[], deta[], res_eta[]);
  fclose (fp1);
  fprintf (fp,
	   "set title 't = %.2f min'\n"
	   "p 'gnuplot' u ($1/1e3):3 w l t 'etap - eta',"
	   "  '' u ($1/1e3):4 w d t 'res'\n"
	   //	   "pause 1\n"
	   , t/60.);
  fflush (fp);
}
#endif

event plots (t = 3600)
{
  foreach()
    fprintf (stderr, "%g %.4g %.4g %.4g\n", x, eta[], u.x[], w[]);
  fprintf (stderr, "\n");
}
