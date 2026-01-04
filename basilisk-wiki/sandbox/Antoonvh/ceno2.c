/**
# Implicit (compact) ENO

For the development of a implicit (compact) WENO scheme, we first
investigate an implicit ENO scheme.  We choose between two compact 3-rd order
accurate one-sided (left and right) schemes and a compact 4th order accurate
centered scheme, depending on a a local smoothness indicator.

~~~gnuplot Non smooth data
set xlabel 'x'
set ylabel 's'
plot 'out' u 1:3
~~~

~~~gnuplot Derivative
set xlabel 'x'
set ylabel 'ds/dx'
plot 'out' u 1:4
~~~
*/
#include "grid/multigrid1D.h"
#include "solve.h"

// (WENO) Smoothness Indicator of Jiang and Shu
#define beta0 (13./12.*sq(s[-1] - 2*s[] + s[1]) + 0.25*sq(s[1] - s[-1]))
#define beta(i) (13./12.*sq(s[-1 + i] - 2*s[i] + s[1 + i]) + 0.25*sq(3*s[] - 4*s[i] + 3*s[2*i]))

// Compact schemes LHS
#define CS (0.25*(ds[-1] + ds[1]) + ds[])
#define RCS (ds[] + 2*ds[1])
#define LCS (ds[] + 2*ds[-1])

// Compact schemes RHS
#define RC (3./(4.*Delta)*(s[1 ] - s[-1 ]))
#define RRS ((-15./6.*s[0] + 2.*s[1] + s[2]/2.)/Delta)
#define LRS ((15./6.*s[0] - 2.*s[-1] - s[-2]/2.)/Delta)

int main() {
  L0 = 2;
  X0 = -1;
  TOLERANCE = 1e-5;
  periodic (left);
  init_grid (N);
  scalar s[], ds[], betaf[];
  // Set signal
  foreach() 
    s[] = x > 0 ? pi/2*x - pi/2 : (cos(pi*x/2));
  // Compute bias 
  foreach() {
    double Betam1 = beta(-1);
    double Beta0 = beta0;
    double Beta1 = beta(1);
    if ((Beta0 < 2*Betam1 && beta0 < 2*Beta1))
      betaf[] = 0;
    else if (Betam1 < Beta1)
      betaf[] = -1;
    else if (Betam1 > Beta1)
      betaf[] = 1;
    else
      betaf[] = 0;
  }
  // The implicit scheme may not converge if we do not ...
  // ... remove stencil "islands": Right sided -> <- Left sided
  scalar to_0[];
  foreach() {
    to_0[] = 0;
    if (betaf[]) {
      if (fabs(betaf[(int)betaf[]] + betaf[]) < 1e-3)
	to_0[] = 1;
    }
  }
  // Bias on levels
  multigrid_restriction({to_0, betaf});
  betaf.restriction = no_restriction;
  foreach_cell()
    if (to_0[])
      betaf[] = 0;

  // Initial guess
  foreach() {
    ds[] = (s[1] - s[-1])/(2.*Delta);
    if (betaf[] == 1)
      ds[] = (s[1] - s[])/(Delta);
    if (betaf[] == -1)
      ds[] = (s[] - s[-1])/(Delta);
  }

  double crit = 0.1; // Tuning this value does not really do much (0 < crit < 1) 
  mgstats ceno = solve(ds,
	((betaf[] > crit) ? RCS : ((betaf[] < -crit) ? LCS : CS)) ,
	((betaf[] > crit) ? RRS : ((betaf[] < -crit) ? LRS : RC)));
  foreach()
    printf ("%g %g %g %g\n", x, betaf[], s[], ds[]);
  
  printf ("# %d %d\n", ceno.i, ceno.nrelax);
}
