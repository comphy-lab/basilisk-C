/**
# Adaptive Ray caster

This program sends pixel coordinates to stdout, reads the correspoding
RGB codes from `stdin`, and writes a ppm image to `stderr`.
 */
#include "utils.h"
#include "my_vertex.h"
#include "adapt_field.h"
vector u;
void refine_r (Point point, scalar s) {
  fine (s, 0, 0, 0) = s[];
  fine (s, 2, 0, 0) = s[1];
  fine (s, 0, 2, 0) = s[0,1];
  fine (s, 2, 2, 0) = s[1,1];
  fine (s, 1, 0, 0) = nodata;
  fine (s, 1, 1, 0) = nodata;
  fine (s, 0, 1, 0) = nodata;
  fine (s, 1, 2, 0) = nodata;
  fine (s, 2, 1, 0) = nodata;
}

void refine_gb (Point point, scalar s) {
  fine (s, 0, 0, 0) = s[];
  fine (s, 2, 0, 0) = s[1];
  fine (s, 0, 2, 0) = s[0,1];
  fine (s, 2, 2, 0) = s[1,1];
}

double interpv (scalar s, double xp, double yp) {
  Point point = locate (xp, yp);
  if (point.level > 0) {
    double w00 = (xp - (x - Delta/2.))*(yp - (y - Delta/2.));
    double w10 = ((x + Delta/2.) - xp)*(yp - (y - Delta/2.));
    double w01 = (xp - (x - Delta/2.))*((y + Delta/2.) - yp);
    double w11 = ((x + Delta/2.) - xp)*((y + Delta/2.) - yp);
    return (w00*s[1,1] + w10*s[0,1] + w01*s[1] + w11*s[])/sq(Delta);
  }
  else return 0;
}

void get_w (scalar * list, scalar w) {
  scalar * wl = list_clone (list);
  scalar s, ws;
  for (s, ws in list,wl)
    wavelet (s, ws);
  for (int l = depth() - 2; l < depth(); l++)
    foreach_coarse_level (l) { 
      double maxw = 0.;
      foreach_child() {
	double a = 0;
	for (scalar ws in wl) 
	  a += fabs(ws[]);
	maxw = max(maxw, a);
      }
      w[] = maxw;
    }
  foreach()  
    w[] = coarse(w,0,0,0);
}

int maxlevel = 9;
double tol = 5;

int main (int argc, char ** argv) {
  if (argc > 1)
    maxlevel = atoi(argv[1]);
  if (argc > 2)
    tol = atof(argv[2]);
#if _OPENMP
  int threads = npe();
#endif
  int s = sizeof(double);
  int Ni = 1 << maxlevel;
  init_grid (1 << (maxlevel - 2));
  scalar w[];
  scalar r[], g[], b[], * rgb = {r, g, b};
  for (scalar s in rgb) {
    s.refine = refine_gb;
    s.prolongation = refine_vert;
    s.restriction = restriction_vert;
  }
  r.refine = refine_r;
  foreach_vertex()
    r[] = nodata;
  int count = 1;
  while (count > 0) {
    count = 0;
    foreach_vertex(reduction(+:count)) 
      if (r[] == nodata) 
	count++;
    fwrite (&count, 1, sizeof(int),  stdout);
    fflush (stdout);
    if (count > 0) {
#if _OPENMP
      omp_set_num_threads(1);
#endif
      foreach_vertex() {
	if (r[] == nodata) {
	  fwrite (&x, s, 1, stdout);
	  fwrite (&y, s, 1, stdout);
	}
      }
      fflush (stdout);
      foreach_vertex() {
	if (r[] == nodata) {
	  for (scalar s in rgb) {
	    unsigned char c;
	    read (STDIN_FILENO, &c, 1);
	    s[] = c;
	  }
	}
      }
#if _OPENMP
      omp_set_num_threads(threads);
#endif
      if (depth() == maxlevel - 2)
	multigrid_restriction(rgb);
      get_w (rgb, w);
      adapt_field (w, tol, 0, maxlevel, 99, list = rgb);
    }
  }
  fprintf (stderr, "P6\n%d %d\n%d\n", Ni, Ni, 255);
  for (int j = 0; j < Ni; j++) {
    double yp = L0 - (j + 0.5)*L0/Ni;
    for (int i = 0; i < Ni; i++) {
      double xp = (i + 0.5)*L0/Ni;
      unsigned char px[3];
      px[0] = interpv (r, xp, yp);
      px[1] = interpv (g, xp, yp);
      px[2] = interpv (b, xp, yp);
      fwrite (px, sizeof(unsigned char), 3, stderr);
    }
  }
  fflush (stderr);
}
