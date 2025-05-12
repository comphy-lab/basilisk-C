/**
# Implicit Burgers equation solver at $\mathrm{CFL} = 100$

~~~gnuplot Solutions obtained with dt $\in$ 0.001, 0.01, 0.1, 1, showing as sharpend and smoothened solutions, respectively
set grid
set size square
set key top left
plot 'out' w l lw 2 t 'Solutions', exp(-x**2) + 1 t 'Initialized'
~~~

We are also curious about the convergence costs and plot the total number of times a cell was "relaxed" (see below) to obtain the finial solution.

~~~gnuplot
set logscale xy
set xr [0.5: 1500]
set ylabel 'Number of cell relaxations'
set xlabel 'Timesteps to solution'
plot 'log' , x*5000 t 'linear'
~~~

The linear line corresponds to about 10 relaxations for every leaf cell per step. 

*/
#include "grid/multigrid1D.h"
#include "utils.h"
#include "run.h"

scalar u[];

double TOLERANCE = 1e-3;
  
int nrelaxs = 0;
int main() {
  L0 = 10;
  X0 = -5;
  periodic (left);
  N = 1 << 9;
  for (DT = 0.001; DT <= 1; DT *= 10)
    run();
}

event init (t = 0) {
  nrelaxs = 0;
  foreach()
    u[] = exp(-sq(x)) + 1.;
}

event timestep (i++) {
  dt = dtnext(DT);
}

/**
## Implicit integration

The implicit scheme,

$$u_{n+1}[\ ] = u_n[\ ] - \mathrm{d}t \ u_{n+1}[\ ]\frac{u_{n+1}[1] - u_{n+1}[-1]}{2\Delta},$$

Can be solved iteratively, asuming that the gradient is relatively steady,

$$u^{k+1}_{n+1} = \frac{u_n[\ ]}{1 + \mathrm{d}t\frac{u^k_{n+1}[1] - u^k_{n+1}[-1]}{2\Delta}}.$$

Apart from the "dont devide by zero condition",

$$\mathrm{d}t\frac{u^k_{n+1}[1] - u^k_{n+1}[-1]}{2\Delta} \neq -1,$$

the iteration does not converge unless there is a sufficiently good initial guess ($k = 0$) *and* sufficient underrelaxation is applied. Fortuntely, the convergence can be accelerated by using the approximation from a coarser level.
*/

event advance (i++, last) {
  scalar ut[];
  foreach()
    ut[] = u[];
  restriction ({u, ut});
  CFL = dt*normf(u).max/(L0/N);
  int minlevel = min(depth(),
		     depth() - (int)log(CFL)/log(2) - 1) ;
  assert (minlevel > 0);
  for (int l = minlevel; l <= depth(); l++) {
    int m = 0;
    double resp = HUGE;
    do {
      if (m++ == 0 && l > minlevel) { 
	boundary_level ({u}, l - 1);
	foreach_level(l)
	  u[] = bilinear(point, u);
      }
      boundary_level ({u}, l);
      CFL = 0;
      foreach_level(l) 
	CFL = max (dt*fabs(u[])/Delta, CFL);
      double alpha = min(.55, 1.2/CFL); // Maybe very small
      foreach_level(l) {
	nrelaxs++;
	u[] = ((1 - alpha)*u[] +
	       alpha*ut[]/(1 + (dt*(u[1] - u[-1])/(2*Delta))));
      }
      resp = -HUGE;
      foreach_level(l)
	resp = max(resp, fabs(ut[] -
			      dt*u[]*(u[1] - u[-1])/(2*Delta) -
			      u[]));
    } while (resp/dt > TOLERANCE/(1 << (depth() - l)) 
	     && m < 100);
  }
}

event stop (t = 1) {
  putchar('\n');
  foreach()
    printf ("%g %g %g %d\n", x, u[], log(dt), nrelaxs);
  fprintf (stderr, "%d %d\n", i, nrelaxs);
  return 1;
}

