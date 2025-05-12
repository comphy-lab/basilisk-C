/**
# Adaptive Ray caster

This program sends pixel coordinates to stdout, reads the correspoding
RGB codes from `stdin`, and writes a ppm image to `stderr`.
 */
#include "utils.h"
scalar ref[], r[], g[], b[], * rgb = {r, g, b};

void refine_ref (Point point, scalar s) {
  foreach_child()
    s[] = nodata;
}

int maxlevel = 9;
double tol = 5;

int main(int argc, char ** argv) {
  if (argc > 1)
    maxlevel = atoi(argv[1]);
  if (argc > 2)
    tol = atof(argv[2]);
  int Ni = 1 << maxlevel;
  init_grid (1 << (maxlevel - 2));
  ref.refine = refine_ref;
  foreach() 
    ref[] = nodata;
  int count = 1;
  while (count > 0) {
    count = 0;
    foreach(reduction(+:count)) 
      if (ref[] == nodata) 
	count++;
    fwrite (&count, 1, sizeof(int),  stdout);
    fflush (stdout);
    int s = sizeof(double);
    if (count > 0) {
      foreach_leaf() { //Z-order sequence
	if (ref[] == nodata) {
	  fwrite (&x, s, 1, stdout);
	  fwrite (&y, s, 1, stdout);
	}
      }
      fflush (stdout);
      foreach_leaf() { //Identical Z-order sequence
	if (ref[] == nodata) {
	  for (scalar s in rgb) {
	    unsigned char c;
	    read (STDIN_FILENO, &c, 1);
	    s[] = c;
	  }
	}
	ref[] = 0.;
      }
      boundary (rgb);
      int d = min (depth() + 1, maxlevel);
      while(adapt_wavelet ({r,g,b}, (double[]){tol, tol, tol}, d, 99).nf);
    }
  }
  fprintf (stderr, "P6\n%d %d\n%d\n", Ni, Ni, 255);
  for (int j = 0; j < Ni; j++) {
    double yp = L0 - (j + 0.5)*L0/Ni;
    for (int i = 0; i < Ni; i++) {
      double xp = X0 + (i + 0.5)*L0/Ni;
      unsigned char px[3];
      Point point = locate (xp, yp);
      if (level == maxlevel) {
	px[0] = r[];
	px[1] = g[];
	px[2] = b[];
      } else {
	px[0] = interpolate_linear (point, (struct _interpolate){r, xp, yp});
	px[1] = interpolate_linear (point, (struct _interpolate){g, xp, yp});
	px[2] = interpolate_linear (point, (struct _interpolate){b, xp, yp});
      }
      fwrite (px, 1, 3, stderr);
    }
  }
  fflush (stderr);
}
