/**
# Poiseuille flow

Due to the usage of 3rd-order accurrate boundaries, The parabolic flow
profile is resolved exactly.

![Quasi steady](poiseuille/ux.mp4)

~~~gnuplot error
set xr [4:64]
set yr [-1e-15:1e-15]
set grid
set xlabel 'N'
set logscale x 2
set ylabel 'Error'
plot 'out' t 'Data'
~~~
 */
#define NOSLIP_TOP (1)
#define NOSLIP_BOTTOM (1)
#include "nsf4t.h"
scalar * tracers = NULL;

double dpdx = 0.08, muv = 1e-2; //Normalized maximum speed

double Poiseuille (double x, double y) {
  return -(dpdx/(2*muv)*y*(y - L0));
}

int main() {
  periodic (left);
  const vector av[] = {dpdx, 0.};
  const scalar muc[] = muv;
  nu = muc;
  a = av;
  for (N = 8; N <= 32; N *= 2)
    run();
}

event init (t = 0) {
  foreach_face(x)
    u.x[] = Gauss6_x(x, y, Delta, Poiseuille);
  boundary ((scalar*){u});
}

event mov (t += 0.5) {
  output_ppm (u.x, file = "ux.mp4", n = 300, min = 0, max = 1.0);
}
	       
event stop (t = 100) {
   double e = 0;
   foreach_face(x)
     e += dv()*fabs(u.x[] - Gauss6_x(x, y, Delta, Poiseuille));
   printf ("%d %g\n", N, e); 
}
