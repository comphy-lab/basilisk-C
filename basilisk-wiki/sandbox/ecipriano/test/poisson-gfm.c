/**
# Test cases for the Poisson--Helmholtz solver with Ghost Fluid Method

We reproduce the test cases proposed by [Liu et al. 2000](#liu2000boundary) to
verify the implementation of the GFM approach to impose interface boundary
conditions in the solution of Poisson--Helmholtz equations.
*/

#include "grid/multigrid1D.h"
#include "poisson-gfm.h"
#include "fractions.h"
#include "run.h"

int maxlevel = 11;

scalar a[], f[];

int main (void) {
  L0 = 1 [0];
  init_grid (1 << maxlevel);

  /**
  ## Case 1: Laplace equation

  We solve the following Laplace equation:
  $$
    \nabla^2 a = 0
  $$

  with:
  $$
    \left[a\right]_\Gamma = 1
  $$
  $$
  \left[\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = 0
  $$
  */

  fraction (f, x - 0.5*L0);

  a[left] = dirichlet (0.);
  a[right] = dirichlet (2.);

  scalar b[];
  foreach()
    b[] = 0.;

  face vector ajump[];
  foreach_face()
      ajump.x[] = 1.;

  poisson_ghost (a, b, f = f, ajump = ajump,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:1]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x <= 0.5) ? x : x + 1

  plot "log" index 0 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */

  /**
  ## Case 2: Laplace equation

  We solve the following Laplace equation:
  $$
    \nabla^2 a = 0
  $$

  with:
  $$
    \left[a\right]_\Gamma = 0
  $$
  $$
  \left[\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = 1
  $$
  */

  a[left] = dirichlet (0.);
  a[right] = dirichlet (1.5);

  foreach() {
    a[] = 0.;
    b[] = 0.;
  }

  face vector bjump[];
  foreach_face()
    bjump.x[] = 1.;

  poisson_ghost (a, b, f = f, bjump = bjump,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:1]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x <= 0.5) ? x : 2.*x - 0.5

  plot "log" index 1 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */

  /**
  ## Case 3: Poisson equation

  We solve the following Poisson equation:
  $$
    \nabla\cdot\left(\beta\nabla a\right) = \begin{cases}
      -2e^x\sin (x) & \text{ if } x < 1 \\
      4(\cos^2(x) - \sin^2(x)) & \text { if } x \geq 1
    \end{cases}
  $$

  with:
  $$
    \left[a\right]_\Gamma = e^x\cos (x) - (\sin^2(x) - \cos^2(x))
  $$
  $$
  \left[\beta\nabla a\cdot\mathbf{n}_\Gamma\right]_\Gamma = -
    e^x(\cos (x) - \sin (x)) - 4\sin (x) \cos (x)
  $$
  */

  L0 = 2 [0];
  init_grid (1 << maxlevel);

  a[left] = dirichlet (exp(x)*cos(x));
  a[right] = dirichlet (sq(sin(x)) - sq(cos(x)));

  foreach()
    f[] = (x < 1.) ? 1. : 0.;

  foreach() {
    a[] = 0.;
    b[] = (x < 1.) ? (-2.*exp(x)*sin(x)) : (4.*(sq(cos(x)) - sq(sin(x))));
  }

  face vector betac[];
  foreach_face()
    betac.x[] = 1.;

  foreach_face()
    ajump.x[] = (exp(x)*cos(x) - sq(sin(x)) + sq(cos(x)));

  foreach_face()
    bjump.x[] = -(exp(x)*(cos(x) - sin(x)) - 4.*sin(x)*cos(x));

  poisson_ghost (a, b, f = f, ajump = ajump, bjump = bjump, alpha = betac,
      tolerance = TOLERANCE/sq(dt), nrelax = 4);

  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel))
    fprintf (stderr, "%g %g %g\n", x, interpolate (f, x), interpolate (a, x));
  fprintf (stderr, "\n\n");

  /**
  ~~~gnuplot
  reset
  set xr[0:2]
  set grid
  set xlabel "x"
  set ylabel "a"
  set key top left

  f(x) = (x <= 1.) ? (exp(x)*cos(x)) : (sin(x)**2 - cos(x)**2)

  plot "log" index 2 u 1:3 w l lw 2 t "Numerical", \
       f(x) w p pt 34 t "Exact"
  ~~~
  */
}

/**
## References

~~~bib
@article{liu2000boundary,
  title={A boundary condition capturing method for Poisson's equation on irregular domains},
  author={Liu, Xu-Dong and Fedkiw, Ronald P and Kang, Myungjoo},
  journal={Journal of computational Physics},
  volume={160},
  number={1},
  pages={151--178},
  year={2000},
  publisher={Elsevier}
}

@article{liu2017second,
  title={A second order ghost fluid method for an interface problem of the Poisson equation},
  author={Liu, Cheng and Hu, Changhong},
  journal={Communications in Computational Physics},
  volume={22},
  number={4},
  pages={965--996},
  year={2017},
  publisher={Cambridge University Press}
}
~~~
*/

