/**
# Flow over a rough embedded topography

![The growth and decay of turbulence with tracer](example_rough/mov.mp4)

~~~gnuplot A plot of the friction velocity
set size ratio -1
plot 'ustar50' matrix with image
~~~
*/
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "lambda2.h"
#include "tracer.h"
#include "diffusion.h"
scalar s[], * tracers = {s}; 
/**
Including `rough-embed.h` switches the embedded boundary to a
rough surface whose roughness length is defined via `Z_0`.
 */
double z0 = 0.01;
#define Z_0 (z0)
#include "rough-embed.h"

face vector muc[];
double ue = 0.05, cse = 0.01; //Refinement criteria

int maxlevel = 6;

int main() {
  N = 32;
  L0 = 2*pi;
  X0 = Z0 = -L0/2.;
  periodic (left);
  periodic (front);
  mu = muc;
  run();
}

event init (t = 0) {
  add_rough_scalar (s, s_val = 1); //<- set boundary conditions for s
  vertex scalar phi[];
  refine (y < 2.1 && y > 0.1 && level < maxlevel);
  foreach_vertex()
    phi[] = -1.1 + 0.5*sin(x) * cos(z*2.) + y;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach()
    u.x[] = cs[] > 0;
  boundary ({u.x, cs, fs});
}

event defaults (i++) {
  foreach_face()
    muc.x[] = fm.x[]/10000.;
  boundary ((scalar*){muc});
}

event tracer_diffusion (i++) {
  diffusion (s, dt, muc);
}

event adapt (i++) {
  cs.prolongation = refine_injection;
  boundary ({cs});
  adapt_wavelet ({cs, u}, (double[]){0.001, ue, ue, ue}, maxlevel, 1.);
  cs.prolongation = fraction_refine;
  boundary ({cs});
}

#define INDEX(i,j) (nx*i + j) // Cartesian index for nx*nx array 

event output_flux_plane (t += 10) {
  int nx = 1 << maxlevel;
  double ustar[nx*nx*3]; //quasi 3D array initialized with high values
  for (int i = 0; i < nx; i++) 
    for (int j = 0; j < nx; j++)
      ustar[INDEX(i,j)] = nodata; 
  foreach() { 
    if (cs[] > 0. && cs[] < 1. && level == maxlevel) {
      coord n = facet_normal (point, cs, fs), p;
      double alpha = plane_alpha (cs[], n);
      plane_area_center (n, alpha, &p);
      int xi = point.i - GHOSTS, zj = point.k - GHOSTS;
      double v[dimension];
      ustar[INDEX(xi, zj)] = u_tau (point, cs, n ,p, (scalar*){u}, v);
    }
  }
#if _MPI //We should reduce all copies of ustar to a single one
  MPI_Barrier (MPI_COMM_WORLD);
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, ustar, nx*nx, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  else
    MPI_Reduce (ustar, NULL, nx*nx, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  
#endif
  //Write it to some space-delimited asci file
  if (pid() == 0) {
    char fname[99];
    sprintf (fname, "ustar%g", t);
    FILE * fp = fopen (fname, "w");
    for (int i = 0; i < nx; i++) { 
      for (int j = 0; j < nx; j++)
	fprintf (fp, "%g ", ustar[INDEX(i,j)]);
      fputc ('\n', fp);
    }
    fclose (fp);
  }
}

event movie (t += 0.2 ; t < 50) {
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});
  foreach()
    l2[] *= cs[];
  view (ty = -0.2, theta = 0.3, phi = 0.4);
  draw_vof ("cs", "fs", fc = {0.1,0.7,0.2});
  squares ("s", n = {0, 0 ,1}, alpha = -L0/2., min = -0.1, max = 1.1);
  cells(n = {0, 0 ,1}, alpha = -L0/2.);
  isosurface ("l2", -0.05);
  save("mov.mp4");
}
/**
We make profiles of the velocity components as a function of distance
to the surface.
 */
#include "../profile6.h"
#include "../frac-dist.h"
event profile_maker (t += 10) {
  char file[99];
  sprintf (file, "prof%g", t);
  vertex scalar phi[];
  distance_to_surface (cs, fs, phi = phi);
  profiles (rough_list, phi, rf = 0.5,
	    fname = file, min = 0.05, max = 2.6);
}
/**
   ~~~gnuplot x-velocity profile
   set yr [0 : 2.3]
   set xr [0 : 1.1]
   set key box top left
   set ylabel 'Distance to surface'
   set xlabel 'Velocity: x-component'
   plot 'prof20' u 2:1 w l lw 2, 'prof50' u 2:1 w l lw 2 ,	\
   'prof90' u 2:1 w l lw 2
   ~~~
*/
