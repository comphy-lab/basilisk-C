/**
# Test for refinement in 2D

~~~gnuplot Average error
set logscale x 2
set logscale y
set grid
set xr [ 17:1000]
set xlabel 'N'
set ylabel 'L_1 error'
set key box bottom left
set size square
plot 'out' u 1:2 t '1^2 stencil', '' u 1:4 t '2^2 stencil',\
 '' u 1:6 t '3^2 stencil', '' u 1:8 t '4^2 stencil', '' u 1:10 t '5^2 stencil',\
0.4*x**-1 lw 2  t'1st order', 2e4*x**-5 lw 2 t '5th order'
~~~

~~~gnuplot Max. error
set logscale x 2
set logscale y
set grid
set xr [ 17:1000]
set xlabel 'N'
set ylabel 'L_{\inf} error'
set key box bottom left
set size square
plot 'out' u 1:3 t '1^2 stencil', '' u 1:5 t '2^2 stencil',\
 '' u 1:7 t '3^2 stencil', '' u 1:9 t '4^2 stencil', '' u 1:11 t '5^2 stencil',\
8*x**-1 lw 2  t'1st order', 1e6*x**-5 lw 2 t '5th order'
~~~
 */

#include "grid/multigrid.h"
#include "higher-order.h"
#include "utils.h"

scalar s[], w[];
int main () {
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2.;
  for (int l = 4; l < 10 ; l++) {
    init_grid (1 << l);
    foreach()
      s[] = exp (-(sq(x) + sq(y) + sq(z)));
    boundary ({s});
    
    printf ("%d ", N);
    
    s.prolongation = refine_1st;
    wavelet (s, w);
    norm f = normf (w);
    printf ("%g %g ", f.avg, f.max);
    
    s.prolongation = refine_2nd;
    wavelet (s, w);
    f = normf (w);
    printf ("%g %g ", f.avg, f.max);
    
    s.prolongation = refine_3rd;
    wavelet (s, w);
    f = normf (w);
    printf ("%g %g ", f.avg, f.max);
   
    s.prolongation = refine_4th;
    wavelet (s, w);
    f = normf (w);
    printf ("%g %g ", f.avg, f.max);
   
    s.prolongation = refine_5th;
    wavelet (s, w);
    f = normf (w);
    printf ("%g %g\n", f.avg, f.max);

    fflush (stdout);
  }
}
