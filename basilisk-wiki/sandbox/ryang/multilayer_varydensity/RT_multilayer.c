/**
# Rayleigh-Taylor instability
*/
#include "grid/multigrid1D.h"
#define NL (nl/2)
#define rho(l) ((l) < NL ? 1.[-3,0,1] : 1.367)
#define G_PHI G
#define G_ETA 0
#include "layered/hydro.h"
#include "layered/nh.h"

#define HALF 1
#include "remap1.h"

const double H0 = 0.108;
const double hi = 0.054;
double ai;

/** Function for writing fields at time t.
*/
int writefields_1D (double t, const char *suffix) {
    char filename1[50], filename2[50], filename3[50], filename4[50];
    sprintf (filename1, "field/ux_%s_t%g", suffix, t);
    sprintf (filename2, "field/uz_%s_t%g", suffix, t);
    sprintf (filename3, "field/h_%s_t%g", suffix, t);
    sprintf (filename4, "field/phi_%s_t%g", suffix, t);
    FILE * fux = fopen (filename1, "w");
    FILE * fuz = fopen (filename2, "w");
    FILE * fh = fopen (filename3, "w");
    FILE * fphi = fopen (filename4, "w");
    foreach() {
      foreach_layer() {
        fprintf (fux, "%g %g \n", x, u.x[]);
        fprintf (fuz, "%g %g \n", x, w[]);
        fprintf (fh, "%g %g \n", x, h[]);
        fprintf (fphi, "%g %g \n", x, phi[]);
      }
    }
    fflush (fux);
    fclose (fux);
    fflush (fuz);
    fclose (fuz);
    fflush (fh);
    fclose (fh);
    fflush (fphi);
    fclose (fphi);
    return 0;
}

int main()
{
  size (0.108 [1]);
  //origin (0.);
  //periodic (right);
  //dimension (nx = hi, ny = H0);
  N = 128;
  nl = 256;
  G = 0.74*9.81;
  DT = 1.e-3;
  theta_H = 1.;
  nu = 5e-6;
#if 0
  lambda_b[] = {HUGE}; // free slip
#endif
ai = 0.0;
  //max_slope = 1.;
  run();
}

/**
The initial layer positions. */

event init (i = 0)
{
  geometric_beta (1./3., true);
  foreach() {
    foreach_layer()
      h[] = point.l < NL ? (hi + 0.0012*cos(2.*pi*x/0.054))/NL : (H0 - hi - 0.0012*cos(2.*pi*x/0.054))/NL;
  }
  vertical_remapping(h, tracers);
}

event field_log (t += 0.01) {
  char *suffix = "matrix";
  writefields_1D (t, suffix);
}

event end (t = 0.4) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}
