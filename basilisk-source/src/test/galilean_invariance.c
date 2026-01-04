/**
# Advection of a rippled interface

We test the ability of the [multilayer solver](/src/layered/hydro.h)
to transport a corrugated interface with a constant velocity profile,
i.e. we check that the solver preserves the property of
[Galilean invariance](https://en.wikipedia.org/wiki/Galilean_invariance).

![Advection of a rippled interface](galilean_invariance/wave.gif)

In the absence of body forces, the interface is initialised with a
sinusoidal shape which is then swept away with a uniform velocity
profile $U$ (this moving interface is therefore not to be confused
with a [inertial-gravity propagating wave](/src/examples/tsunami.c)).

The resulting solid-body like motion is an exact solution of the
[multilayer set of equations with
the non-hydrostatic corrections](/src/layered/nh.h), with $\phi = 0$. */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

/**
We non-dimensionalise the problem with the advection velocity $U$
and the wavelength $\lambda$, and we set the amplitude of the wave
to one tenth of the channel depth. */

#define wavelength 1.
#define U 1.
#define wavenumber (2.*pi/wavelength)
#define depth 0.5
#define amplitude (depth/10.)
#define endoftime (1.*wavelength/U)

double maxwaveerror;

/**
## Main function

The boundary conditions are set to periodic. The test is done in a
zero-gravity environment. The default minmod slope limiter is turned
off to avoid blunting off the wave crests. */

int main()
{
  size (1 [0]); // space is dimensionless
  periodic (right);
  G = 0.;
  gradient = NULL;

  breaking = 0.; // fixme: without this the computation diverges for N = 256

  for (N = 64; N <= 256; N *= 2)
    for (nl = 1; nl <= 16; nl *= 2)
      run();
}

/**
## Initialisation

The initial shape of the interface is a simple sine, and the velocity field
is set to a constant velocity $U = 1$ and $w = 0$ everywhere. */

event init (i = 0)
{
  foreach() {
    double H = depth + amplitude*cos(wavenumber*x);
    foreach_layer() {
      h[] = H/nl;
      u.x[] = U;
    }
  }
  maxwaveerror = 0.;
}

/**
## Outputs

We keep track of the amplitude of the wave vs time, and for a typical
case ($N = 128$ and $\text{nl} = 4$) we export deeper checks on the
non-hydrostatic pressure $\phi$ and velocity $u$. */

event monitoring (t += endoftime/100; t = 0.; t <= endoftime)
{
  stats s = statsf(eta);
  double etaamp = fabs((s.max - s.min)/(2*amplitude) - 1);
  if (N == 128 && nl == 4) {
    char name[80];
    sprintf (name, "output-N-%d-nl-%d", N, nl);
    static FILE * fpout = fopen (name, "w");

    double absdevU = 0., absPhi = 0.;
    foreach()
      foreach_layer() {
        if (fabs(u.x[] - U) > absdevU)
	  absdevU = fabs(u.x[] - U);
	if (fabs(phi[]) > absPhi)
	  absPhi = fabs(phi[]);
      }
    fprintf (fpout, "%g %g %g %g\n", t, absPhi, absdevU/U, etaamp);
    assert (absPhi < 2e-13 && absdevU/U < 1e-14);
  }
  maxwaveerror = max(maxwaveerror, etaamp);
}

/**
The reference file of this test is just based on the overall maximum wave
amplitude error. */

event errorlog (t = end)
{
  fprintf (stderr, "%d %d %.3e\n", N, nl, maxwaveerror);
}

/**
## Movie

This is how we export the headline animated gif. */

event movie (i += 5; t <= (wavelength/U))
{
  if (N == 128 && nl == 4) {
    static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
    if (i == 0)
      fprintf (fp,
	       "set term pngcairo enhanced font ',10' size 800,250\n"
	       "set size ratio -1\n"
	       "unset key\n");
    fprintf (fp,
	     "set output 'plot%04d.png'\n"
	     "set title 't = %.2f'\n"	     
	     "p [%g:%g][0:]'-' u 1:3:2 w filledcu lt rgb \"#993498DB\", "
	     " '' u 1:2 w l lw 4 lt rgb \"#3498DB\" ",
	     i/5, t, X0, X0 + L0);
    int i = 4;
    foreach_layer()
      fprintf (fp, ", '' u 1:%d w l dt 2 lt rgb \"gray50\"", i++);
    fprintf (fp, "\n");
    foreach (serial) {
      double H = 0.;
      foreach_layer()
	H += h[];
      fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
      double z = zb[];
      foreach_layer() {
	fprintf (fp, " %g", z);
	z += h[];
      }
      fprintf (fp, "\n");
    }
    fprintf (fp, "e\n\n");
    fflush (fp);
  }
}

event moviemaker (t = end)
{
  if (N == 128 && nl == 4)
    system ("mogrify -format gif plot*.png && "
	    "gifsicle --delay 10 --loop plot*.gif > wave.gif && "
	    "rm -f plot*.*");
}

/**
# Results

For the typical case $N = 128$ and $\text{nl} = 4$ the time history
of the Poisson error $\|\phi\|_\infty$ and of any deviation of
the velocity to the setpoint $\|u-U\|_\infty$ -- which would be responsible
for phase errors -- are represented and are observed to be of the
order of the machine precision.

~~~gnuplot Left: error on $\phi$ vs time estimated as $\|\phi\|_\infty$ (Poisson error). Right: error on horizontal velocity $\|u-U\|_\infty$ vs time (phase error)
set term @SVG size 960,240 font ',10'
set size ratio 2./(1. + sqrt(5.))
set logscale y
set grid
set style line 1 \
    linecolor rgb '#774F38' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 2 \
    linecolor rgb '#E08E79' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 3 \
    linecolor rgb '#F1D4AF' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 4 \
    linecolor rgb '#ECE5CE' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.
set style line 5 \
    linecolor rgb '#C5E0DC' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.

set ylabel 'Error (phi)'
set xlabel "t"
set key left
set multiplot layout 1, 2
plot [0:1][1e-16:1e-10]'./output-N-128-nl-4' u 1:2 t 'Error on phi' with lp ls 1
set ylabel 'Error (u)'
plot [0:1][1e-16:1e-10]'./output-N-128-nl-4' u 1:3 t 'Error on u' with lp ls 2
unset multiplot
~~~

For the same case ($N = 128$ and $\text{nl} = 4$) the relative error in
wave amplitude (dissipation error) is monitored through time, and is seen
to be of the order of $10^{-4}$.

~~~gnuplot Relative error of the wave amplitude with time (dissipation error).
set term @SVG size 640,320 font ',10'
set ylabel 'Error (amplitude)'
plot './output-N-128-nl-4' u 1:4 t 'Relative error on amplitude' with lp ls 3
~~~

## Convergence

We track the value of the relative error on wave amplitude for various
number of gridpoints $N$ and number of layers $\text{nl}$. For up to
16 stacked layers the results are undistinguishable, and the error is observed
to decrease as $N^{-2}$, indicating a second-order precision.

~~~gnuplot Variation of the relative error on wave amplitude with horizontal resolution, for different (fixed) number of layers.
set xlabel 'Number of grid points'
set ylabel 'Error'
set key bottom left
set logscale x 2
set label 1 "N^{-2}" at 180,5e1*180**-2 font ',12' textcolor rgb '#E08E79' offset 0,1
filter1 = '< awk ''{ if ($2 == 1) print $1, $2, $3}'' '
filter2 = '< awk ''{ if ($2 == 2) print $1, $2, $3}'' '
filter4 = '< awk ''{ if ($2 == 4) print $1, $2, $3}'' '
filter8 = '< awk ''{ if ($2 == 8) print $1, $2, $3}'' '
filter16 = '< awk ''{ if ($2 == 16) print $1, $2, $3}'' '
plot filter1.'log' u 1:3 t "nl = 1" with linespoints linestyle 1, \
     filter2.'log' u 1:3 t "nl = 2" with linespoints linestyle 2, \
     filter4.'log' u 1:3 t "nl = 4" with linespoints linestyle 3, \
     filter8.'log' u 1:3 t "nl = 8" with linespoints linestyle 4, \
     filter16.'log' u 1:3 t "nl = 16" with linespoints linestyle 5, \
     [96:96<<2]5e1*x**-2 t '' with lines linestyle 2 lw 3
~~~
*/
