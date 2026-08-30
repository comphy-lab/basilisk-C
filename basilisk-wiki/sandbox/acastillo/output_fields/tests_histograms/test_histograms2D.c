/**
# Testing `probability_distribution_2D()` on two independent tanh fields

`u(x,y) = 0.5*(1 + tanh(x/epsu))` and `v(x,y) = 0.5*(1 + tanh(y/epsv))` each
depend on a different coordinate, so on a square domain uniformly sampled in
(x,y) they are statistically independent. Each marginal has the same
closed-form CDF as in `test_histograms.c`,

  CDF_u(u) = (x(u) - X0)/L0,  x(u) = epsu*atanh(2u-1)
  CDF_v(v) = (y(v) - Y0)/L0,  y(v) = epsv*atanh(2v-1)

so the joint (bin-averaged) PDF factors as

  pdf(u,v) = pdf_u(u) * pdf_v(v)

with each marginal bin-averaged PDF given by (CDF(b+db) - CDF(b))/db, as in
the 1D test.

The stored joint PDF can be visualized as a heatmap
~~~gnuplot Joint probability density function
set xlabel 'u'
set ylabel 'v'
set pm3d map
set palette rgbformulae 33,13,10
splot "histogram2D.asc" u 2:3:4 notitle
~~~

*/

#define LEVEL 9

#include "run.h"
#include "acastillo/output_fields/histograms2D.h"

double epsu = 0.15;
double epsv = 0.2;

int main(){

  L0 = 1.0;
  X0 = Y0 = -L0 / 2;
  N = 1 << LEVEL;
  init_grid(N);

  scalar u[], v[];
  foreach(){
    u[] = 0.5*(1. + tanh(x/epsu));
    v[] = 0.5*(1. + tanh(y/epsv));
  }
  boundary({u, v});

  int nx = 55, ny = 55;
  double p_xrange[nx], p_yrange[ny], p_pdf[nx*ny];
  probability_distribution_2D(u, v, 0, 1, 0, 1, nx, ny, p_xrange, p_yrange,
                               p_pdf, "histogram2D.asc", true);

  /**
  Verify the stored joint PDF against the closed-form product of marginals.
  The report goes to stderr, i.e. to the `log` diffed against
  `test_histograms2D.ref`. */
  if (pid() == 0)
    system ("python3 ../test_histograms2D.py --histogram histogram2D.asc "
            "--epsu 0.15 --epsv 0.2 --X0 -0.5 --Y0 -0.5 --L0 1.0 1>&2");
}
