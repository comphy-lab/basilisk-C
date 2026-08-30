/**
# Testing `probability_distribution_1D()` on a flat tanh interface

`c(x,y) = 0.5*(1 + tanh(y/eps))` is uniform in x and monotonically
increasing in y, with a smooth but sharp transition of width `eps` at
y = 0. Since y is uniform on [Y0, Y0+L0] and c(y) is a known, invertible
function of y, the CDF has a closed form:

  CDF(c) = (y(c) - Y0)/L0,  y(c) = eps*atanh(2c-1)

and the checker derives the expected, bin-averaged PDF from it as
(CDF(c+db) - CDF(c))/db, rather than the continuum density
eps/(4*L0*c*(1-c)) -- the two only agree in the limit of many grid cells
per bin (see test_histograms.py).

*/

#define LEVEL 10

#include "run.h"
#include "acastillo/output_fields/histograms1D.h"

double eps = 0.15;

int main(){

  L0 = 1.0;
  X0 = Y0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  scalar c[];
  foreach()
    c[] = 0.5*(1. + tanh(y/eps));
  boundary({c});

  int nbin = 110;
  double p_range[nbin], p_pdf[nbin], p_sum[nbin];
  probability_distribution_1D(c, 0, 1, nbin, p_range, p_pdf, p_sum,
                               "histogram.asc", true);

  /**
  Verify the stored PDF and CDF against the closed form above. The report
  goes to stderr, i.e. to the `log` diffed against `test_histograms.ref`. */
  if (pid() == 0)
    system ("python3 ../test_histograms.py --histogram histogram.asc "
            "--eps 0.15 --Y0 -0.5 --L0 1.0 1>&2");
}
