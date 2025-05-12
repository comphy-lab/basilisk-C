/**
# 6th-order accurate finite-volume advection in 2D

![Test case](fv62d/s.mp4)

~~~gnuplot
set xr [12:192]
set logscale x 2
set logscale y 
set xlabel 'N'
set ylabel 'L_1 error'
set size square
set grid
plot 'out' t 'data', 2e8*x**(-6) t '6th order'
~~~
 */
#define RKORDER (4)
#include "lsrk.h"
#include "higher-order.h"
#include "utils.h"

double u_x (double x, double y) {
  return  -y;
}

double u_y (double x, double y) {
  return x;
}

void advfv6 (scalar * sl, scalar * dsdtl) {
  boundary (sl);
  scalar s = sl[0], dsdt = dsdtl[0];
  face vector sf[];
  scalar sv[];
  vector flx[];
  // Get scalar-averaged value on faces
  cell_to_face (s, sf);
  // Get Scalar value on vertices
  face_to_vertex (sf, sv);
  // Flux vector value at vertex
  foreach() {
    foreach_dimension()
      flx.x[] = sv[]*u_x(x - Delta/2, y - Delta/2);
  }
  boundary ((scalar*){flx});
  // Get flux average on faces, reuse `sf`
  vertex_vector_to_face (flx, sf);
  boundary_flux ({sf});
  // The cell-averaged tendency is the flux divergence
  foreach() {
    dsdt[] = 0;
    foreach_dimension()
      dsdt[] -= (sf.x[1] - sf.x[0])/Delta;
  }
}

scalar s[];

int main () {
  foreach_dimension()
    periodic (left);
  L0 = 20;
  X0 = Y0 = -L0/2;
  for (N = 16; N <= 256; N *= 2) {
    DT = 2./(N);
    run();
  }
}

event init (t = 0) {
  foreach() {
    double a = 0;
    foreach_child() foreach_child() foreach_child()
      a += exp(-(sq(x - 2) + sq(y)));
    s[] = a/64.;
  }
}

event mov (t += pi/20) {
  output_ppm (s, file = "s.mp4", n = 300, min = -0.1, max = 1.1);
}

event advance (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({s}, dt, advfv6);
}

event stop (t = 2.*pi) {
  double e = 0;
  scalar es[];
  foreach() {
    double a = 0;
    foreach_child() foreach_child() foreach_child()
      a += exp(-(sq(x - 2) + sq(y)));
    es[] = (s[] - a/64.);
    e +=  sq(Delta)*fabs(es[]);
  }
  output_ppm (es, file =  "err.mp4", n = 512);
  printf ("%d %g\n", N, e);
  return 1;
}
