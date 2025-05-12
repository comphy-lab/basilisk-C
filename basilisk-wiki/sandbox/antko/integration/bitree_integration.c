/**
# Integral of an oscillatory function

We are interested in computing the integral of a `nasty' function that is
fastly oscillating (lots of destructive interferences) and whose frequency
drifts in the domain so as to approach the grid spacing.

Let's take $\sin(x^3)$ for example in the range $x \in \left[0, 3 \pi^{1/3}\right]$.
This function looks like this:

~~~gnuplot The function represented in the entire range and its sampled version on regular 128 points grid.
set term @SVG size 1280,320
unset tics
unset key
set samples 10000
plot [0: 3 * pi**(1./3.)] [-1.2:1.2] sin(x**3) lw 2 lt rgb "#3498DB", \
     'sampled.txt' u 1:2 pt 7 ps 0.75 lc rgb "#FF5733"
~~~

The integral $I$ of this function in this range has an analytical expression
involving exponential integrals, and which evaluates numerically to
$I \simeq 0.46375354024413091536$.

In the following we test, compare and monitor different strategies for
integral computation.

#### Headers and declarations

For (cubic) spline interpolation and Gauss-Kronrod quadratures we will rely
on the GNU Scientific Library.

*/

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#pragma autolink -lgsl -lgslcblas

/**

Our starting point will be to use a regular 128-point grid ($N = 2^7$), covering
the domain $x \in \left[0, 3 \pi^{1/3}\right]$. We enforce the fact that our
function vanishes on the edges of the domain with the dirichlet directive.

*/
int LEVEL = 7;

scalar sampled_function[];
sampled_function[left] = dirichlet(0.);
sampled_function[right] = dirichlet(0.);

/**

# Integration strategies

We test four strategies. In the first two, we suppose that we have only access
to the function value at the (128) regularly spaced collocation points. This is
expected to be fast and not very precise. In the remaining two, we will look
at the precision gain brought about using refinement, i.e. additional evaluations
of the function at intermediary points.

More precisely we use:

1. A simple trapezoidal integration rule on the regular grid,
2. A cubic spline interpolation using the original data points to perform the integration,
3. A bitree adaptative remeshing procedure to refine the integral's value (still using the trapezoidal rule), and
4. A Gauss-Kronrod quadrature in each of the 128 cells.

*/

void simple_trapezoidal_integration ();
void spline_interp_and_integrate();
void adaptive_bitree_trapezoidal_integration();
void gauss_kronrod_integration();

int main ()
{
  simple_trapezoidal_integration();
  spline_interp_and_integrate();
  adaptive_bitree_trapezoidal_integration ();
  gauss_kronrod_integration();
}

/**
#### Toolbox

We first define our nasty function $f(x) = \sin(x^3)$ and a reference value
for the integral (computed with Mathetica using Levin's rule, 20 digits
working precision and successive subdivisions to reach an accuracy goal
of twenty digits).

*/

trace
double oscillating_func (double pos, void * params) {
  long unsigned int *neval = (long unsigned int *) params;
  double y = pow(pos, 3.);
  (*neval)++;
  return sin(y);
}

void sample_on_grid (scalar field, void * params) {
  foreach()
    field[] = oscillating_func (x, params);
  boundary ({field});
}

const double reference_value = 0.46375354024413091536;

/**
## Trapezoidal integration, no adaptation

We first suppose that the function is known at the initial collocation points
and that there is no way to refine our knowledge of the function at
intermediary points (because e.g. it is too costly).

We implement a simple trapezoidal rule to estimate $I$ using this constant
grid.

Remember that each sampled value is at the center of the cell, so that we have
two truncated trapezoids in each cell :


$$I_\mathrm{trapezoid} = \sum_\mathrm{cells} \left(\frac{f[-1] + f[]}{2} + f[]\right) \frac{\Delta}{4} + \left(\frac{f[1] + f[]}{2} + f[]\right)
\frac{\Delta}{4}$$

*/

/**
~~~gnuplot Left : whenever the function is smooth, the trapezoidal rule captures very well the area under the curve. Right : As soon as the function is more *wiggly* the trapezoids lose track of the area.
set term @SVG size 800,320
set tics
unset key
set samples 10000
set multiplot layout 1, 2 ;
plot [0: 1.] [-1.2:1.2] sin(x**3) lw 2 lt rgb "#3498DB", \
     'sampled.txt' u 1:2 pt 7 ps 0.75 lc rgb "#FF5733",\
     'sampled.txt' u 1:2 with filledcurves y1=0 lt rgb "#DDFF5733",\
     'sampled.txt' u 1:2 with lines dt 2 lt rgb "#FF5733"
plot [3 * pi**(1./3.)-1.: 3 * pi**(1./3.)] [-1.2:1.2] sin(x**3) lw 2 lt rgb "#3498DB", \
     'sampled.txt' u 1:2 pt 7 ps 0.75 lc rgb "#FF5733",\
     'sampled.txt' u 1:2 with filledcurves y1=0 lt rgb "#DDFF5733",\
     'sampled.txt' u 1:2 with lines dt 2 lt rgb "#FF5733"
unset multiplot
~~~
*/
trace
double trapz (scalar field) {
  double sum = 0.;
  foreach()
    sum += (0.5 * (field[-1] + field[1]) + 3. * field[]) * Delta * 0.25;
  return sum;
}

trace
void simple_trapezoidal_integration () {
  origin (0.);
  double xmax = 3. * pow(pi, 1./3.);
  size (xmax);
  init_grid (1 << LEVEL);
  long unsigned int call_to_function = 0;
  sample_on_grid (sampled_function, &call_to_function);
  double trapezoidal_approx = trapz (sampled_function);
  FILE * fpsample = fopen ("sampled.txt", "w+");
  for(double x = xmax / 256; x < xmax; x += xmax/128) {
    fprintf (fpsample, "%.8g %.8g\n", x, interpolate(sampled_function,x));
  }
  fclose (fpsample);
  printf ("--------------------------------------------\n");
  printf ("> Simple trapezoidal integration\n");
  printf ("Number of calls to function : %lu\n", call_to_function);
  printf ("Reference value : %.15g\n", reference_value);
  printf ("Trapezoidal approximation : %.15g\n", trapezoidal_approx);
  printf ("Absolute error : %.4e\n", fabs(trapezoidal_approx-reference_value));
  printf ("Relative error : %.6g %%\n", 100 * fabs((trapezoidal_approx-reference_value)/reference_value));
}

/**
This simple implementation yields a result with less than a percent of error,
which is not too bad, clearly sufficient for a quick estimate, but not for
precise computation:

`--------------------------------------------`\
`> Simple trapezoidal integration`\
`Number of calls to function : 128`\
`Reference value : 0.463753540244131`\
`Trapezoidal approximation : 0.459787485747487`\
`Absolute error : 3.9661e-03`\
`Relative error : 0.855207 %`\
`--------------------------------------------`

Comparatively to others methods, the regular trapezoidal rule is super fast,
as expected:

| `calls` | `total` | `self`| `% total` | `function` |
|   ---:  |  ---:   |  ---: |      ---: |      :--- |
| `1`     | `0.00`  | `0.00`| `0.1%`    | `simple_trapezoidal_integration():bitree_integration.c:165` |
| `3`     | `0.00`  | `0.00`| `0.7%`    | `init_grid():/Users/antko/basilisk/src/grid/tree.h:1710`    |
| `2`     | `0.00`  | `0.00`| `0.4%`    | `trapz():bitree_integration.c:144`                          |

*/

/**
## Cubic spline interpolation and integration (no adaptation)

We next move on to a cubic spline interpolation and integration procedure.
To do so, we make use of the GNU Scientific Library framework.

*/

trace
void spline_interp_and_integrate () {
  int N = 1 << LEVEL; // For comparison purposes, we reinitialize the grid
  init_grid (N);
  double xmin = 0.;
  double xmax = 3. * pow(pi, 1./3.);

  long unsigned int call_to_function = 0;
  sample_on_grid (sampled_function, &call_to_function);

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  const gsl_interp_type *spline_type= gsl_interp_cspline;
  gsl_spline *spline = gsl_spline_alloc (spline_type, N+2);

  double *x_values = calloc(N+2, sizeof(double));
  double *f_values = calloc(N+2, sizeof(double));

/**

We export the position and values of the function in a single vector.
This is *so not* Basilisk friendly.

*/
 int count = 0;
 for(double x = xmax / 256; x < xmax; x += xmax/128) {
    count++;
    x_values[count] = x;
    f_values[count] = interpolate(sampled_function,x);
  }

/**

and we add the end points

*/

  x_values[0] = xmin;
  x_values[N+1] = xmax;
  f_values[0] = f_values[N+1] = 0.;

/**

The spline interpolation and integration are obtained with these single calls

*/

  gsl_spline_init (spline, x_values, f_values, N+2);

  double spline_approx = gsl_spline_eval_integ(spline, xmin, xmax, acc);

  FILE * fpspline = fopen ("spline.txt", "w+");
  for (double x = 0. ; x <= xmax ; x += xmax/10000) {
    fprintf (fpspline, "%.8g %.8g\n", x, gsl_spline_eval(spline, x, acc));
  }
  fclose (fpspline);

  printf ("--------------------------------------------\n");
  printf ("> Cubic spline interpolation and integration\n");
  printf ("Number of calls to function : %lu\n", call_to_function);
  printf ("Reference value : %.15g\n", reference_value);
  printf ("Spline approximation : %.15g\n", spline_approx);
  printf ("Absolute error : %.4e\n", fabs(spline_approx-reference_value));
  printf ("Relative error : %.6g %%\n", 100 * fabs((spline_approx-reference_value)/reference_value));

/**

Finally we make some cleanup

*/

  free (x_values);
  free (f_values);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

/**
~~~gnuplot The cubic spline allow to nicely recover the variation of the function inbetween the collocation points.
set term @SVG size 1240,320
set tics
unset key
set samples 10000
plot [3 * pi**(1./3.)-1.: 3 * pi**(1./3.)] [-1.2:1.2] sin(x**3) lw 2 lt rgb "#3498DB", \
     'sampled.txt' u 1:2 pt 7 ps 0.75 lc rgb "#FF5733",\
     'spline.txt' u 1:2 with filledcurves y1=0 lt rgb "#DDFF5733",\
     'spline.txt' u 1:2 with lines dt 2 lt rgb "#FF5733"
~~~
*/

/**
From the previous graph we could expect a dramatic increase in precision, but actually
the gain is "only" of one order in magnitude, for a comparable cost :

`--------------------------------------------`\
`> Cubic spline interpolation and integration`\
`Number of calls to function : 128`\
`Reference value : 0.463753540244131`\
`Spline approximation : 0.463281377523503`\
`Absolute error : 4.7216e-04`\
`Relative error : 0.101813 %`\
`--------------------------------------------`

| `calls` | `total` | `self`| `% total` | `function` |
|   ---:  |  ---:   |  ---: |      ---: |      :--- |
| `1`     | `0.00`  | `0.00`| `1.0%`    | `spline_interp_and_integrate():bitree_integration.c:290` |

*/


/**

## Adaptive bitree trapezoidal integration

Here we make use of Basilisk adaptive remeshing capabilities to refine the
sampling where needed, thanks to `adapt_wavelet()`

*/

trace
void adaptive_bitree_trapezoidal_integration () {
  long unsigned int call_to_function = 0;
  printf ("--------------------------------------------\n");
  printf ("> Adaptive bitree trapezoidal integration\n");
  printf("Cycles of adaptation to reach tol = 1e-15\n");
  astats s;
  int count = 0;
  int cell_change = -1;
  while (cell_change != 0) {
    s = adapt_wavelet ({sampled_function}, (double []){1e-15}, maxlevel = LEVEL+4, minlevel = LEVEL-4); // note: I used maxlevel = LEVEL+12 for producing the logs
    printf ("# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
    sample_on_grid (sampled_function, &call_to_function);
    count++;
    cell_change = s.nf;
  }
  printf("Needed %i cycles.\n", count);
  double trapezoidal_approx = trapz (sampled_function);

  #if 0
  FILE * fpadaptive = fopen ("adaptive.txt", "w+");
  foreach() {
    if (x >= (3 * pow(pi,1./3.)-0.1))
      fprintf (fpadaptive, "%.8g %.8g\n", x, sampled_function[]);
  }
  fclose (fpadaptive);
  #endif

  printf ("Number of calls to function : %lu\n", call_to_function);
  printf ("Reference value : %.15g\n", reference_value);
  printf ("Trapezoidal approximation : %.15g\n", trapezoidal_approx);
  printf ("Absolute error : %.4e\n", fabs(trapezoidal_approx-reference_value));
  printf ("Relative error : %.6g %%\n", 100 * fabs((trapezoidal_approx-reference_value)/reference_value));
}

/**

The plot of the refined points (not shown) exhibits a so large density of points
that they cover entirely the analytical representation. More information
can be obtained from these logs:

`Note : these logs were obtained with maxlevel = LEVEL+12 (impossible to reach on the server)`\
`--------------------------------------------`\
`> Adaptive bitree trapezoidal integration`\
`Cycles of adaptation to reach tol = 1e-15`\
`# refined 128 cells, coarsened 0 cells`\
`# refined 256 cells, coarsened 0 cells`\
`# refined 512 cells, coarsened 0 cells`\
`# refined 1024 cells, coarsened 0 cells`\
`# refined 2048 cells, coarsened 0 cells`\
`# refined 4096 cells, coarsened 0 cells`\
`# refined 8192 cells, coarsened 0 cells`\
`# refined 16384 cells, coarsened 0 cells`\
`# refined 32768 cells, coarsened 0 cells`\
`# refined 65536 cells, coarsened 0 cells`\
`# refined 131072 cells, coarsened 0 cells`\
`# refined 262144 cells, coarsened 0 cells`\
`# refined 0 cells, coarsened 0 cells`\
`Needed 13 cycles.`\
`Number of calls to function : 1572608`\
`Reference value : 0.463753540244131`\
`Trapezoidal approximation : 0.463753539905172`\
`Absolute error : 3.3896e-10`\
`Relative error : 7.30903e-08 %`\
`--------------------------------------------`

we see that the gain in precision is of 7 orders of magnitude w.r.t. the
trapezoidal method, but that there is a massive price to pay with **1.5m**
evaluations of function.

This high computational price is also seen in the monitoring:

| `calls` | `total` | `self`| `% total` | `function` |
|   ---:  |  ---:   |  ---: |      ---: |      :--- |
| `1575680`     | `0.14`  | `0.14`| `29.5%`    | `oscillating_func():bitree_integration.c:105` |
| `28`     | `0.12`  | `0.12`| `23.9%`    | `boundary():/Users/antko/basilisk/src/grid/cartesian-common.h:368` |
| `13`     | `0.24`  | `0.12`| `24.5%`    | `adapt_wavelet():/Users/antko/basilisk/src/grid/tree-common.h:292` |
| `1`     | `0.48`  | `0.10`| `19.9%`    | `adaptive_bitree_trapezoidal_integration():bitree_integration.c:385` |

Note: to be fair, this procedure is not the only calling `oscillating_func()` but it
accounts for more than 99.8% of them.
*/

/**
## Gauss-Kronrod

Eventually we test a [Gauss-Kronrod](https://en.wikipedia.org/wiki/Gauss–Kronrod_quadrature_formula)
quadrature formula on each of the initial cells. This type of formula yields exact results for
polynomials of high order (the meaximum order depending on the number of points involved). The advantage
of this procedure is that additional points can be recruited if the target precision is not met,
without requiring to compute again the function value at the previous points.

*/

trace
void gauss_kronrod_integration () {
  int N = 1 << LEVEL;
  init_grid (N);
  long unsigned int call_to_function = 0;
  sample_on_grid (sampled_function, &call_to_function);
  gsl_function F;
  F.function = &oscillating_func;
  F.params = &call_to_function;
  double result, abserr;
  long unsigned int neval;
  double gauss_kronrod_approx = 0.;

  foreach() {
    gsl_integration_qng(&F, x-Delta/2., x+Delta/2., 1e-15, 0, &result, &abserr, &neval);
    gauss_kronrod_approx += result;
  }

  printf ("--------------------------------------------\n");
  printf ("> Gauss-Kronrod quadrature\n");
  printf ("Number of calls to function : %lu\n", call_to_function);
  printf ("Reference value : %.15g\n", reference_value);
  printf ("Gauss-Kronrod approximation : %.15g\n", gauss_kronrod_approx);
  printf ("Absolute error : %.4e\n", fabs(gauss_kronrod_approx-reference_value));
  printf ("Relative error : %.6g %%\n", 100 * fabs((gauss_kronrod_approx-reference_value)/reference_value));
  printf ("Last number of evals : %lu\n", neval);
  printf ("Last absolute error : %.15g\n", abserr);
  printf ("--------------------------------------------\n");
}

/**

Here are the results:

`--------------------------------------------`\
`> Gauss-Kronrod quadrature`\
`Number of calls to function : 2816`\
`Reference value : 0.463753540244131`\
`Gauss-Kronrod approximation : 0.463753540244131`\
`Absolute error : 3.8858e-16`\
`Relative error : 8.37898e-14 %`\
`Last number of evals : 21`\
`Last absolute error : 2.69239643679262e-16`\
`--------------------------------------------`

| `calls` | `total` | `self`| `% total` | `function` |
|   ---:  |  ---:   |  ---: |      ---: |      :--- |
| `1`     | `0.00`  | `0.00`| `0.0%`    | `gauss_kronrod_integration():bitree_integration.c:423` |

Note also that this procedure accounts for $\simeq$ 0.2% of the function calls,
corresponding to $\simeq$ 0.05% of the total computational time. This is way
better than the bitree computation.

# Conclusion

This is a clear win for the Gauss-Kronrod procedure which allows, for a reasonable
number of calls to the function to gain 13 orders of magnitude in precision
(and actually reach the machine precision!).
*/
