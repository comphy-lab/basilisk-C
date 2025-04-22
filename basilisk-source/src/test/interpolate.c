/* interpolation */

#include "utils.h"

scalar v[];

int main (int argc, char ** argv)
{
  size (1.[0]);
  origin (-0.5, -0.5);
  for (int n = 8; n <= 64; n *= 2) {
    init_grid (n);

    refine (level == log2(n) && sq(x) + sq(y) > sq(0.49));
    refine (level <= log2(n) + 1 && sq(x) + sq(y) > sq(0.55));

    foreach()
      v[] = cos(2.*pi*x)*cos(2.*pi*y);

    //    FILE * fp = fopen("error","w");
    double emax = 0.;
    int ni = 4*n + 7;
    double Delta = 1./ni;
    for (int i = 0; i < ni; i++) {
      double x = X0 + Delta*i;
      for (int j = 0; j < ni; j++) {
	double y = Y0 + Delta*j;
	double e = fabs (cos(2.*pi*x)*cos(2.*pi*y) - interpolate (v, x, y));
	//	fprintf (fp, "%g %g %g\n", x, y, e);
	if (e > emax) emax = e;
      }
      //      fprintf (fp, "\n");
    }
    //    fclose (fp);

    fprintf (stderr, "%d %g\n", n*4, emax);
    if (n == 16)
      output_cells (stdout);
  }
}
