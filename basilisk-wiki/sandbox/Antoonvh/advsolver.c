/**
# A simple solver for $\frac{\partial c}{\partial t} = -\frac{\partial c}{\partial x}$

We employ an equidistant grid and the generic Runge-Kutta time-integrator to solve a one dimensional advection equation.
*/
#define BGHOSTS 2
#include "grid/multigrid1D.h"
#include "runge-kutta.h"

/**
   A marco for the analytical solution is `DEFINE`d, two $c$ fields are
   declared and the corresponding boundary conditions are set.
*/

#define SOLUTION(t) (exp(-sq(x + 0.5 - t)))
scalar c2[], c4[];
double t; //Global for time dependent boundary conditions
c2[left] = dirichlet(SOLUTION(t));
c2[right] = dirichlet(SOLUTION(t));
c4[left] = dirichlet(SOLUTION(t));
c4[right] = dirichlet(SOLUTION(t));

/**
These functions compute the right-hand-side of the time evolution
equation. `dcdt2` and `dcdt4` are function with a 2nd and 4th order
spatial accuracy, respectively. Their implementations may seem a bit
specific, but it mush match the syntax imposed by more generic
formulation in `runge-kutta.h`.
*/

static void dcdt2 (scalar * s, double t, scalar *kl){//2nd order accurate
  scalar m  = s[0];
  scalar k = kl[0];
  boundary(all);
  foreach()
    k[] = -(m[1] - m[-1]) / (2*Delta);
}

static void dcdt4 (scalar * s, double t, scalar *kl){//4th order accurate
  scalar m  = s[0];
  scalar k = kl[0];
  boundary(all);
  foreach()
    k[] = -(-m[2] + 8*m[1] - 8*m[-1] + m[-2]) / (12*Delta);
}

/**
We set the case up and solve for different grids.
 */
int main(){
  L0 = 10;
  X0 = -L0/2.;
  double DT = 1./1000.;
  FILE * fpc = fopen("convergence", "w");
  for (int j = 5; j <= 9; j++){
    init_grid(1 << j);
    foreach(){ //Initialize
      c2[] = SOLUTION(0);
      c4[] = c2[];
    }
    for (t = 0; t < 1; t += DT){ //Time loop
      runge_kutta ({c2}, t, DT, dcdt2, 4);
      runge_kutta ({c4}, t, DT, dcdt4, 4);
    }
    char fname[99];
    sprintf(fname, "data%d", j);
    FILE * fpd = fopen(fname, "w");
    double err2 = 0;
    double err4 = 0;
    foreach(){
      err2 += fabs(c2[] - SOLUTION(1)) * Delta;
      err4 += fabs(c4[] - SOLUTION(1)) * Delta;
      fprintf(fpd, "%g\t%g\t%g\t%g\t%g\n", x, c2[], c4[], fabs(c2[] - SOLUTION(1)), fabs(c4[] - SOLUTION(1)));
    }
    fprintf(fpc, "%d\t%g\t%g\n", j, err2, err4);
  }
}

/**
## Results

The convergence is OK
~~~gnuplot
set xr [4.5:9.5]
set logscale y
set ylabel 'L1 norm'
set xlabel 'refinement level'
plot 'convergence' u 1:2 t '2nd order solver' ,\
     'convergence' u 1:3 t '4rd order solver'
~~~
*/