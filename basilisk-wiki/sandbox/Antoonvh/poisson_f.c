/**
# Errors for the Poisson solver.

This page concerns the convergence of the Poisson solver for $f(x,y)$,

$$\nabla^2 f = s.$$

Since [sine-rich
solutions](http://gerris.dalembert.upmc.fr/gerris/tests/tests/poisson.html) are
*not* suitable for testing, we design an other test case where $s$ is
a dipole,

$$s = e^{-(x - 2)^2 - \left(\frac{y}{2}\right)^2} - e^{-\left(\frac{x
+ 2.}{2}\right)^2 - (y - 1)^2} $$

The solution ($f$) only exists at the mercy of boundary conditions for
$f$.
*/
#include "utils.h"
#include "poisson.h"
scalar s[], f[];

#define SOURCE (exp(-sq(x - 2.) - sq(y/2.)) - exp(-sq((x + 2.)/2.) - sq(y - 1)))

f[top] = dirichlet (0.);
f[bottom] = dirichlet (0.);
/**
The cell-averaged source term is approximated with 64 points.
*/
double source (Point point) {
  double a = 0;
  foreach_child() foreach_child() foreach_child() a += SOURCE;
  return a/pow(pow(2., (double)dimension), 3.);
}
/**
# A reference solution

It may be hard to find an analytical expression for the
cell-averaged solution $f$ and hence we compute a reference solution
via a superior-resolution grid. The solution on all relevant levels
is then stored in an array.
*/
#define OFFSET(l) ((1 << (2*(l)))*((l) + 1))
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i+_O) + (1 << l)*(j+_O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level)) 

int maxlevel = 9;         //The maximum resolution for analysis
int reference_level = 11;  //The superior resolution for reference 

double * sol;

void obtain_reference_solution() {
  sol = (double*)malloc (OFFSET (maxlevel)*sizeof(double));
  init_grid (1 << reference_level);
  foreach() 
    s[] = source (point);
  poisson (f, s);
  restriction ({f});
  foreach_cell()
    if (level <= maxlevel)
      sol[INDEX] = f[];
  output_ppm (s, file = "s.png", n = 350);
  output_ppm (f, file = "f.png", n = 350);
}
/**
 The test system looks like this:

![The source term ($s$)](poisson_f/s.png)

![The Solution ($f$)](poisson_f/f.png)

## The setup
   
We set a large domain, a small `TOLERANCE` and a periodic solution in
the left-right direction. Then a reference solution is computed.
*/
int main() {
  L0 = 20.;
  X0 = Y0 = -L0/2;
  TOLERANCE = 1e-7;
  periodic (left);
  obtain_reference_solution();
/**
Solutions are approximated on increasingly refined
grids. Furthermore, some error data is written.
*/
  FILE * fp = fopen ("error_data", "w");
  for (int l = 4; l <= maxlevel; l++) {
    init_grid (1 << l);
    foreach() 
      s[] = source (point);
    poisson (f, s);
    boundary ({f});
    scalar dfdx[];
    foreach()
      dfdx[] = (f[1] - f[-1])/(2*Delta);
    boundary ({dfdx});
    double err = 0;
    foreach() { 
      double e = f[] - sol[INDEX];
      err += fabs(e)*sq(Delta);
      fprintf (fp, "%d\t%g\t%g\t%g\n", l, e/sq(Delta), s[],
	       (dfdx[0,1] - dfdx[0,-1])/(2*Delta));
    }
    printf ("%d\t%g\n", N, err);
    if (l == maxlevel) {
      output_ppm (s, file = "s.png", n = 512);
      output_ppm (f, file = "f.png", n = 512);
    }
  }
  fclose (fp);
  free (sol);
}

/**
## Results 

We check the converegence by plotting the total error
versus `N`.

~~~gnuplot It's seconder order
set xr [10:800]
set logscale xy
set xlabel 'N'
set ylabel 'Total error' 
plot 'out' t 'data', 100*x**(-2) t 'second order convergence'
~~~

This motivates to plot the scaled error against the cell averaged
value of $s$ (`s[]`). 

~~~gnuplot Error behaviour amongst levels
reset
set term pngcairo
set output 'plot1.png'
set xlabel 'e/sq(Delta)'
set ylabel 'level + s[]' 
set grid
plot 'error_data' u 2:($1+$3):1 palette t 'data', -x*16.5 + 7 lw 2 t 'linear'
~~~

The error could maybe also scale with
$\frac{\mathrm{d}^2f}{\mathrm{d}x\mathrm{d}y}$,

~~~gnuplot Error behaviour amongst levels
reset
set term pngcairo
set output 'plot2.png'
set xlabel 'e/sq(Delta)'
set ylabel 'level + ddf/dxdy[]' 
set grid
plot 'error_data' u 2:($1+$4):1 palette t 'data'
~~~

but is does not.  
*/
