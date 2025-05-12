/**
## MWE for surface convection
A wave profil in a domain with periodic conditions is convected with a constant velocity from left to right. The convection of surfacic concentration and the height are monitored. The maximum values should remain constant.
*/

#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../surface.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
double ak = 0.35;
#define LEVEL 7
#define layers 4
#define T_END 10.
#define DELTA_T (T_END/100.)
#define MY_TOLERANCE 1e-11

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.
#define Gamma0 1.

/**
## Main function
There is no diffusion, no adsorption, no viscosity, no gravity*/
double initial_mass = 0.;
stats STi_surf;
stats STi_area;
stats STi_H;

int main()
{
  L0 = L;
  origin (0.);
  periodic (right);
  N = 1 << LEVEL;
//  TOLERANCE = MY_TOLERANCE;
  nl = layers;
  G = g_;
  for (N = 16; N <= 256; N *= 2) {
    gradient = NULL;
    initial_mass = 0;
    run();
  }
}

/**
## Initialisation 
A wave is convected from left to right at a constant velocity. Bulk and surface concentrations are inhomogeneous*/
#include "test/stokes.h"

event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    double H = h_ + wave(x, 0);
    foreach_layer() {
      h[] = H/nl;
  	  u.x[] = 1.;
      eta[] += h[];
    }
  }
  boundary(all);
  foreach() 
    M[] = Gamma0 * area_m(point) * (1 + wave(x, 0));
  boundary({M});    
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = 0.; t <= T_END) {
  scalar area_s[];
  foreach ()
    area_s[] = area_m(point);
  boundary({area_s});
  if (i == 0) {
    STi_surf = statsf(M);
    STi_H = statsf(eta);
    STi_area = statsf(area_s);  
  }
 
  char name[80];
  sprintf (name, "surf_mass-N-%d", N);
  static FILE * fp = fopen (name, "w");
  sprintf (name, "profile-N-%d", N);
  static FILE * fpp = fopen (name, "w");
  scalar surf[];
    foreach()
      surf[] = M[]/area_m(point);
  stats ST_surf = statsf(surf); 
  stats ST_M = statsf(M);   
  stats ST_H = statsf(eta);
  stats ST_area = statsf(area_s);
  initial_mass = (initial_mass == 0) ? ST_M.sum : initial_mass;
  fprintf (fp, "%g %g %g %g %g %g %g\n", t, ST_M.sum - initial_mass, (ST_surf.max), (ST_H.max - STi_H.max)/(STi_H.max - h_), (ST_area.max - STi_area.max)/(STi_area.max-1), (ST_surf.min), (ST_H.min)/(STi_H.min - h_));
  
  foreach()
    fprintf(fpp, "%g %g %g %g %g\n", t, eta[], M[], M[]/area_s[], area_s[]);
}

event error (t = end)
{
  static FILE * err = fopen ("error.txt", "w");
  stats ST_H = statsf(eta);
  fprintf (err, "%g %g\n", N/L0, -(ST_H.max - STi_H.max)/(STi_H.max - h_));
}

/**
# Movie
To see the movement of the wave */
event movie (i += 3)
{
if(N == 128) {
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
	     "set size ratio -1\n");  
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "p [%g:%g][0:1]'-' u 1:(-1):2 w filledcu lc 3 t ''",
	   i/3, t/(k_/sqrt(k_)), X0, X0 + L0);
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
if(N == 128) {
  system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}
}


/** ##Results
~~~gnuplot Difference between total mass and initial total mass
set terminal @PNG enhanced size 320,320 font ",8"
set output 'sum.png'
set xlabel "t"
set key left
plot \
  './surf_mass-N-256' u 1:2 t 'sum' w l
~~~
Total mass is well conserved, the error is the machine error.

~~~gnuplot Relative error for the minimum and maximum surface concentrations.
set output 'Gamma.png'
plot \
  './surf_mass-N-256' u 1:3 t 'max' w l, \
  './surf_mass-N-256' u 1:6 t 'min' w l
~~~

~~~gnuplot Relative error for the maximum area.
set output 'Area.png'
plot \
  './surf_mass-N-256' u 1:5 t 'max' w l
~~~

~~~gnuplot Relative error for minimum and maximum heights.
set output 'Heights.png'
plot \
  './surf_mass-N-256' u 1:4 t 'max' w l, \
  './surf_mass-N-256' u 1:7 t 'min' w l
~~~

~~~gnuplot Relative error for the maximum height as a function of the number of grid points.
set output 'Error.png'
set xlabel 'Number of grid points'
set ylabel 'Relative error'
set logscale y
set logscale x 2
set grid
plot [5:600][1e-3:1]'error.txt' t "Basilisk" w lp
~~~
*/
