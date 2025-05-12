/**					
# Compact upwind advection

The compact upwind advection scheme is tested for functions with
varying degrees of smoothness. 

1. noise
2. Heaviside function 
3. Tent function
4. Compact squared cosine
5. Compact cubed cosine
6. Gaussian 

~~~gnuplot Initial `solutions`
set key right outside
set xlabel 'x'
set xr [-6:4]
set yr [ -1.1:1.1]
plot 'log' u 1:2 t 'noise', '' u 1:3 w l t 'Heaviside',  '' u 1:4 w l  t 'Triangle', '' u 1:5 w l t 'cos^2', '' u  1:6 w l t 'cos^3', '' u 1:7  w l t 'Gaussian'
~~~

~~~gnuplot Solutions at t = 1
plot 'out' u 1:2 t 'noise', '' u 1:3 w l t 'Heaviside',  '' u 1:4 w l  t 'Triangle', '' u 1:5 w l t 'cos^2', '' u  1:6 w l t 'cos^3', '' u 1:7  w l t 'Gaussian'
~~~

~~~gnuplot Total variance
reset
set key right outside
set xlabel 't'
plot 'data' u 1:2 t 'noise' w l, 'data' u 1:4 w l t 'Heaviside', '' u 1:6 w l  t 'Triangle',\
 '' u 1:8 w l t 'cos^2', '' u  1:10 w l t 'cos^3', '' u 1:12  w l t 'Gaussian'
~~~

~~~gnuplot Total variation
set logscale y
plot 'data' u 1:3 t 'noise' w l, 'data' u 1:5 w l t 'Heaviside', '' u 1:7 w l  t 'Triangle',\
 '' u 1:9 w l t 'cos^2', '' u  1:11 w l t 'cos^3', '' u 1:13  w l t 'Gaussian'
~~~
*/

#include "grid/multigrid1D.h"
#include "higher-order.h"
#include "lsrk.h"

vector u[];

void advection (scalar * sl, scalar * dsl) {
  scalar s, ds;
  for (s, ds in sl, dsl) {
    vector grad[];
    compact_upwind ({s}, {grad}, u);
    foreach()
      ds[] = -u.x[]*grad.x[];
  }
}

scalar n[], H[], tr[], sqcos[], cbcos[], smt[];
scalar * tracers = {n, H, tr, sqcos, cbcos, smt};

int main() {
  periodic (left);
  L0 = 10;
  X0 = -6;
  N = 64;
  DT = L0/N;
  run();
}

event init (t = 0) {
  foreach() {
    u.x[] = 1;
    n[] = noise();
    H[]  = fabs(x - 1) < 1 ? 1 : 0;
    tr[] = x < -1 ? 0 : x < 0 ? 1 + x : x < 1 ? 1 - x : 0;
    sqcos[] = fabs(x + 1) < pi/2. ? sq(cos(x + 1)) : 0;
    cbcos[] = fabs(x + 2) < pi/2. ? cube(cos(x + 2)) : 0;
    smt[] = exp(-sq(x + 3));
    fprintf (stderr, "%g ", x);
    for (scalar s in tracers)
      fprintf (stderr, "%g ", s[]);
    fprintf (stderr, "\n");
  }
}

event analysis (i++) {
  static FILE * fp = fopen("data", "w"); 
  fprintf (fp, "%g ", t);
  for (scalar s in tracers) {
    double V = 0, VR = 0;
    foreach() {
      V += sq(s[]);
      VR += fabs(s[] - s[-1]);
    }
    fprintf (fp, "%g %g ", V, VR);
  }
  fprintf (fp, "\n");
}

event advance (i++) {
  dt = dtnext (DT);
  A_Time_Step (tracers, dt, advection);
  boundary (tracers);
}

event stop (t = 1) {
  foreach() {
    printf ("%g ", x);
    for (scalar s in tracers)
       printf ("%g ", s[]);
    printf ("\n");
  }
}
