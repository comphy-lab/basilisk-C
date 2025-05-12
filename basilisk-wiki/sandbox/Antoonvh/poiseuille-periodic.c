/**
# A Poiseuille flow test case

On this page we test the accuracy of the *Stokes* solver to represent
a periodic Poiseuille flow over a domain $x \in [-1, 1]$ forced with $\frac{1}{\rho}\frac{\mathrm{d}P}{\mathrm{d}y} = 1$:

$$u_y = 0.5 \left( 1-x \right) ^2$$
 */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

/**
   A macro to obtain the cell-averaged analytical solution:
*/

#define sol(x) (0.5*(1. - sq(x) - sq(Delta)/12.))

int main(){
  L0 = 2.;
  X0 = -L0/2.;
  periodic (top);
  stokes = true;
  TOLERANCE = 1e-5;
  /**
     We set Neumann conditions corresponding to the analytical
     solution and run with increasing grid resolutions. 
  */
  u.t[left] = neumann(-1);
  u.t[right] = neumann(-1);
  for (N = 2; N <= 64; N *= 2)
    run();
}

event init (t = 0) {
  DT = 0.01;
  const face vector g[] = {0., 1.};
  a = g;
  const face vector muc[] = {1.,1.};
  mu = muc;
  foreach()
    u.y[] = sol(x);
}

event profile (t = 6) {
  static FILE * fp = fopen("output", "w");
  double e = 0;
  foreach()
    e += fabs(u.y[] - sol(x)) * sq(Delta);
  fprintf (fp, "%d %g \n", N, e);
}

/**

The 2nd order solver is able to *exactly* represent the second order
polynomial solution:

~~~gnuplot mind the Y-axis
set logscale x 2
set xr [ 1:128]
set yr [-0.00000001:0.00000001]
set xlabel 'N'
set ylabel 'Error'
plot 'output' u 1:2
~~~

*/