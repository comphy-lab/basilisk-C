/**
# Wave advection
A sinusidal wave with periodic conditions is advected horizontally at a constant velocity from left to right in absence of any forces. 
The wave profil is supposed to stay unchanged along the simulation. 
![Animation of the wave advection](advection_simple_patch/movie.mp4)
*/

/**
## Includes
Tihs test case uses the multilayer solver.
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "nh_patch.h"
#include "layered/remap.h"

/**
## Parameters*/
#define k_ (2.*pi)
#define h_ 0.5
#define ah  h_/10.
#define T_END 3.
#define U 1.

/**
## Main function
The default minmod slope limiter is turned off to avoid cutting off in wave crests.
*/
int main()
{
  periodic (right);
  G = 0.;
  gradient = NULL;
  
  nl = 4;
  for (N = 16; N <= 64; N *= 2) {
    run();
  }
  
  N = 128;
  for (nl = 1; nl <= 16; nl *= 2)
    run();                
}

/**
## Initialisation 
The sinusoidal wave is convected horizontally from left to right at a constant velocity $U=1$ and $w=0$. 
This initial condition is solution of the Navier-Stokes equations.
*/
scalar H_i[];
event init (i = 0)
{
  foreach() {
    double H = h_ + ah*cos(k_*x);
    foreach_layer() {
      h[] = H/nl;
      u.x[] = U;
    }
  }
  
  foreach() {
    H_i[] = zb[];
    foreach_layer()
      H_i[] += h[];
  }
}

/**
## Outputs
The profile is monitored every $t + 1/U$ to check the invarience.
Besides, the standard deviation of the pressure $\phi$ and the horyzontal velocity are extracted, as well as the total height $\eta$ at every time steps.
*/
event profiles (t += 1./U; t = 0.; t <= T_END)
{
  static FILE * fpp = fopen ("profiles", "w");
  foreach() {
    fprintf (fpp, "%g %g %g\n", x, eta[], eta[] - H_i[]);
  }
  fprintf (fpp, "\n\n");
}

event output (t += T_END/100; t = 0.; t <= T_END)
{
  char name[80];
  sprintf (name, "output-N-%d-nl-%d", N, nl);
  static FILE * fpout = fopen (name, "w");

  double devU = statsf(u.x).stddev/U;
  double devPhi = statsf(phi).stddev;
  stats s = statsf(eta);
  double etaamp = (s.max - s.min)/(2*ah) - 1;
  fprintf (fpout, "%g %g %g %g\n", t, devPhi, devU, etaamp);
}

event error (t = end)
{
  double devU = statsf(u.x).stddev/U;
  stats s = statsf(eta);
  double etaamp = (s.max - s.min)/(2*ah) - 1;
  fprintf (stderr, "%d %d %g %g\n", N, nl, devU, etaamp);
}


/**
# Movie */
event movie (i += 5)
{
  if (N == 128 && nl == 4) {
    static FILE * fp = popen ("gnuplot", "w");
    if (i == 0)
      fprintf (fp, "set term pngcairo font ',9' size 800,250;"
  	     "set size ratio -1\n");  
    fprintf (fp,
  	   "set output 'plot%04d.png'\n"
  	   "set title 't = %.2f'\n"
  	   "p [%g:%g][0.4:0.6]'-' u 1:(-1):2 w filledcu lc 3 t ''",
  	   i/5, t, X0, X0 + L0);
    fprintf (fp, "\n");
    foreach_leaf() {
      double H = 0.;
      foreach_layer()
        H += h[];
      fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
      fprintf (fp, "\n");
    }
    fprintf (fp, "e\n\n");
    fflush (fp);
  }
}

event moviemaker (t = end)
{
  if (N == 128 && nl == 4) {
    system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
  }
}

/** ##Results
The standard deviation of $\Phi$ and of the horyzontal velocity are drawn for the case N = 128, nl =4. 
At the first time-steps, Phi is not exactly null, which induce numerical dissipation and reduce slightly the horyzontal velocity.

~~~gnuplot Standard deviation of the pressure $Phi$ over time
set terminal @PNG enhanced size 320,320 font ",8"
set output 'Phi.png'
set xlabel "t"
set key left
plot \
  './output-N-128-nl-4' u 1:2 t 'Std phi' w l
~~~

~~~gnuplot Standard deviation of the horyzontal velocity $u$
set output 'Du.png'
plot \
  './output-N-128-nl-4' u 1:3 t 'Std u' w l
~~~

~~~gnuplot Relative error of the wave amplitude with time.
set output 'Amplitude.png'
plot \
  './output-N-128-nl-4' u 1:4 t 'H_max' w l
~~~

The error on the horyzontal velocity is drawn for different vertical and horyzontal resolutions. 
The velocity deviance do not depend on the horyzontal grid resolution but increases linearly with the number of layers.

~~~gnuplot Standard deviation of the horyzontal velocity $u$ as a function of horyzontal resolution.
set output 'Error_Grid.png'
set xlabel 'Number of grid points'
set ylabel 'Std U'
set logscale y
set logscale x 2
set grid
plot [10:150][1e-6:1e-4]'log' u 1:3 t "Basilisk" w p
~~~

~~~gnuplot Standard deviation of the horyzontal velocity $u$ as a function of vertical resolution.
set output 'Error_Layers.png'
set xlabel 'Number of layers'
set ylabel 'Std U'
set logscale y
set logscale x 2
set grid
plot [0.5:20][1e-6:1e-4]'log' u 2:3 t "Basilisk" w p
~~~

*/