/**
# Convergence of the face fraction field

We wish to assess the convergence of the face fraction calculation by
comparing the face fraction field, reconstructed from the volume
fraction, with the benchmark value obtained from the intialization. */

#include "fractions.h"
#include "utils.h"

scalar c[], f[];
face vector cs[], fs[];

#define circle(x,y,z,R) (sq(R) - sq(x) - sq(y) - sq(z))

int main (void) {
  for (int maxlevel = 6; maxlevel <= 9; maxlevel++) {
    init_grid (1 << maxlevel);
    solid (c, cs, circle (x, y, z, 0.3));
    fraction (f, circle (x, y, z, 0.3));

    face_fraction (f, fs);

    int np = 0;
    double norm1 = 0.;
    foreach_face (reduction(+:norm1) reduction(+:np)) {
      np++;
      norm1 += fabs (cs.x[] - fs.x[]);
    }
    norm1 /= np;

    fprintf (stderr, "%d %.4g\n", N, norm1);
  }
}

/**
~~~gnuplot 2D case
set xlabel "Resolution"
set ylabel "Error"

set xr[2**5:2**10]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "log" u 1:2 via a,b

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)

plot "log" title "norm1", \
     f(x) w l title ftitle(a, b)
~~~
*/
