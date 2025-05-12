/**
# A test

This is a test for the `adapt_field()` function. On this page the aim is to reproduce the behaviour of `adapt_wavelet()`
*/

#include "grid/bitree.h"
#include "adapt_field.h"

scalar s[], w[];

int main(){
  FILE * fp1 = fopen ("adapt_field", "w");
  FILE * fp2 = fopen ("adapt_wavelet", "w");
  FILE * fp3 = fopen ("adapt_field_c", "w");
  FILE * fp4 = fopen ("adapt_wavelet_c", "w");
  L0 = 8;
  X0 = -L0/2;
  init_grid (1 << 3);
  do {
    foreach()
      s[] = exp(-sq(x));
    boundary({s});
  } while (adapt_wavelet ({s}, (double[]){0.015}, 8).nf);
  foreach()
    fprintf(fp2, "%g\t%g\n", x, s[]);
  refine (level < 3);
  unrefine (level >= 3);
  do {
    foreach()
      s[] = exp(-sq(x));
    boundary({s});
    wavelet (s, w);
  } while (adapt_field (w, 0.015, 0.01, 8).nf);
  foreach()
    fprintf(fp1, "%g\t%g\n", x, s[]);
/**

## Results

~~~gnuplot The grid structures are identical
set xr [-4:4]
set yr [-0.1:1.1]
set xlabel 'x'
set ylabel 's[ ]'
plot 'adapt_field' pt 4, 'adapt_wavelet'
~~~

we can also check the coarsening behaviour:
*/
  refine (level < 12);
  do {
    foreach()
      s[] = exp(-sq(x));
    boundary({s});
  } while (adapt_wavelet ({s}, (double[]){0.015}, 8).nc);
  foreach()
    fprintf(fp4, "%g\t%g\n", x, s[]);
  refine (level < 12);
  do {
    foreach()
      s[] = exp(-sq(x));
    boundary({s});
    wavelet (s, w);
  } while (adapt_field (w, 0.015, 0.01, 8).nc);
  foreach()
    fprintf(fp3, "%g\t%g\n", x, s[]);
/**
~~~gnuplot The grid structures are identical
set xr [-4:4]
set yr [-0.1:1.1]
set xlabel 'x'
set ylabel 's[ ]'
plot 'adapt_field_c' pt 4, 'adapt_wavelet_c'
~~~

*/


}