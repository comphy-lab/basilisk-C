/**
# Testing `probability_distribution_2D_weighted()` on two independent tanh fields

Same fields as `test_histograms2D.c`, `u(x,y) = 0.5*(1 + tanh(x/epsu))` and
`v(x,y) = 0.5*(1 + tanh(y/epsv))`, but the histogram is restricted to a
region-of-interest `f = 1` for `x < xf`, `f = 0` otherwise. Restricting the x
range renormalizes the u-marginal's CDF over (X0, xf) instead of the full
domain (X0, X0+L0); the v-marginal is unaffected since f does not depend on
y and v's distribution over the (unrestricted) y-range is unchanged. The
joint (bin-averaged) PDF still factors as

  pdf(u,v) = pdf_u(u) * pdf_v(v)

with

  CDF_u(u) = (x(u) - X0)/(xf - X0),  clamped to [0,1]
  CDF_v(v) = (y(v) - Y0)/L0

The stored joint PDF can be visualized as a heatmap
~~~gnuplot Joint probability density function, restricted to x < xf
set xlabel 'u'
set ylabel 'v'
set pm3d map
set palette rgbformulae 33,13,10
splot "histogram2D_weighted.asc" u 2:3:4 notitle
~~~

*/

#define LEVEL 9

#include "run.h"
#include "acastillo/output_fields/histograms2D.h"

double epsu = 0.15;
double epsv = 0.2;
double xf = 0.2;

int main(){

  L0 = 1.0;
  X0 = Y0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  scalar u[], v[], f[];
  foreach(){
    u[] = 0.5*(1. + tanh(x/epsu));
    v[] = 0.5*(1. + tanh(y/epsv));
    f[] = (x < xf);
  }
  boundary({u, v, f});

  int nx = 55, ny = 55;
  double p_xrange[nx], p_yrange[ny], p_pdf[nx*ny];
  probability_distribution_2D_weighted(f, u, v, 0, 1, 0, 1, nx, ny, p_xrange,
                                        p_yrange, p_pdf,
                                        "histogram2D_weighted.asc", true);

  /**
  Verify the stored joint PDF against the closed-form product of marginals.
  The report goes to stderr, i.e. to the `log` diffed against
  `test_histograms2D_weighted.ref`. */
  if (pid() == 0)
    system ("python3 ../test_histograms2D.py --histogram histogram2D_weighted.asc "
            "--epsu 0.15 --epsv 0.2 --X0 -0.5 --Y0 -0.5 --Lu 0.7 --Lv 1.0 1>&2");
}
