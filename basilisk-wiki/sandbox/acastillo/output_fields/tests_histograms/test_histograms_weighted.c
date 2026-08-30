/**
# Testing `probability_distribution_1D_weighted()` on a flat tanh interface

Same field as `test_histograms.c`, `c(x,y) = 0.5*(1 + tanh(y/eps))`, but the
histogram is restricted to a region-of-interest `f = 1` for `y < yf`, `f = 0`
otherwise. Since c(y) is monotonically increasing and uniform in x, the
weighted CDF has a closed form too:

  CDF(c) = (y(c) - Y0)/(yf - Y0),  clamped to [0,1],  y(c) = eps*atanh(2c-1)

i.e. the same closed form as the unrestricted case, just renormalized by the
restricted region's extent (yf - Y0) instead of the full domain L0.

*/

#define LEVEL 10

#include "run.h"
#include "acastillo/output_fields/histograms1D.h"

double eps = 0.15;
double yf = 0.2;

int main(){

  L0 = 1.0;
  X0 = Y0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  scalar c[], f[];
  foreach(){
    c[] = 0.5*(1. + tanh(y/eps));
    f[] = (y < yf);
  }
  boundary({c, f});

  int nbin = 110;
  double p_range[nbin], p_pdf[nbin], p_sum[nbin];
  probability_distribution_1D_weighted(f, c, 0, 1, nbin, p_range, p_pdf,
                                        p_sum, "histogram_weighted.asc", true);

  /**
  Verify the stored CDF against the closed form above. The report goes to
  stderr, i.e. to the `log` diffed against `test_histograms_weighted.ref`. */
  if (pid() == 0)
    system ("python3 ../test_histograms.py --histogram histogram_weighted.asc "
            "--eps 0.15 --Y0 -0.5 --L0 0.7 1>&2");
}
