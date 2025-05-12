/**
# Solving an integral equation

We aim to find the integral of

$$ S = \int e^{-x^2} \mathrm{d}x = \frac{\sqrt{\pi}}{2}\mathrm{erf}(x) + C $$

We use `solve.h` with various discretizatoins.

~~~gnuplot
set xlabel 'x'
set key left top
set grid
plot 'out' u 1:3 t 'Left handed', '' u 1:4 t 'Right handed', '' u 1:5 t 'Central ', '' u 1:6 w l lw 2 t 'erf(x)'
~~~

Note that the `Central` stencil solution is not invertible, and hence
causes a crash if we would call it...

Further, the right-handed stencil has poor convergence due to the
direction of iterations.

 */
#include "grid/multigrid1D.h"
#include "solve.h"
double b;
int main() {
  L0 = 6;
  X0 = -L0/2.;
  init_grid (N);
  scalar s[], S[], S1[], S2[];
  b = 0.5*sqrt(pi);
  S[left] = dirichlet (-b);  //LH stencil
  S1[right] = dirichlet (b); //RH stencil
  S2[left] = dirichlet (-b);
  S2[right] = dirichlet (b);
  foreach() {
    s[] = exp(-sq(x));
    S[] = S1[] = S2[] = b*x/(L0/2.);
  }
  boundary (all);
  solve (S, (S[] - S[-1])/Delta, s);   // LH
  solve (S1, (S1[1] - S1[])/Delta, s); // RH
  // solve (S2, (S2[1] - S2[-1])/(2.*Delta), s); // Central stencil is not invertible
  foreach()
    printf ("%g %g %g %g %g %g\n", x, s[], S[], S1[], S2[], b*erf(x));
}
