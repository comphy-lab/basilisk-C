/**
# Advection of noise

It can be useful to be aware of the limiting behaviour of the ["B-C-G"
scheme](/src/bcg.h). As such we setup an ABC-flow and advect random noise.

![Evolution of a near-zero contour](conservation/s.mp4)

Inspiration in taken from [this](https://github.com/microhh/microhh/tree/master/cases/conservation) $\mu$-hh test case.
 */
#include "grid/multigrid3D.h"
#include "advection.h"
#include "view.h"

scalar s[], *tracers = {s};

int main(){
  foreach_dimension()
    periodic (left);
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0/2;
  N = 64;
  run();
}

event init (t = 0) { //A-B-C flow
  foreach_face(x)
    uf.x[] = cos(y) + sin(z);
  foreach_face(y)
    uf.y[] = sin(x) + cos(z);
  foreach_face(z)
    uf.z[] = cos(x) + sin(y);
  foreach()
    s[] = noise();
}

event mov (t += 0.02; t <= 10) {
  isosurface ("s", 0.02);
  save ("s.mp4");
}

event quantitative (i += 2) 
  printf("%g %g %g\n", t, fabs(statsf(s).sum), normf(s).rms);

/**
The scheme is conservative for the first order moment, but dissipates
the second order moment.  

~~~gnuplot
set xlabel 'time' 
set logscale y
plot 'out' u 1:2 w l lw 2 t 'abs sum' , 'out' u 1:3 w l lw 3 t 'rms'
~~~

 */
