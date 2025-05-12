/**
# Convergence on a tree grid

We consider the convergence of the numerical errors as more cells are
applied according to the `adapt_wavelet()` algorithm. We compute the
tendency for the diffusion equation (Laplacian), $$\frac{\partial s}{\partial t} =
\nabla^2 s,$$ for a Gaussian blob.

~~~gnuplot L_1 error norm
set logscale x 2
set logscale y
set ylabel 'L_1 Error'
set xlabel 'sqrt (Cells)' 
set key bottom left box
set size square
set grid
set xr [6:800]
set yr [0.001:5]
plot 'out' pt 7 t 'Adapt wavelet()', 'log' pt 5 t 'NxN', 500*x**-2 t 'second order' 
~~~

~~~gnuplot Max error
set logscale x 2
set logscale y
set ylabel 'L_{inf} Error'
set xlabel 'sqrt (#Cells)' 
set key bottom left box
set size square
set grid
set xr [6:800]
set yr [0.0001:5]
plot 'out' u 1:3 ps 1 pt 7 t 'Adapt wavelet()', 'log' u 1:3 ps 1 pt 5 t 'NxN', 500*x**-2 t 'second order' 
~~~

Whoa, Lets check some meshes:

![$\approx 30^2$ cells,](diffusion_error_fix/s889.png)

![$\approx 344^2$ cells.](diffusion_error_fix/s118030.png)

*/
#include "view.h"
#include "higher-order.h"
scalar s[];

/**
The analytical solution
 */
double t0 = 0.25;
#define sqR (sq(x - 0.1234) + sq(y + 0.1357))
#define Gaussian (exp(-sqR/(4*(t + t0)))/(4*pi*(t+t0)))
#define Gaussiand (-(Gaussian*(-sqR + 4*(t + t0)))/(4*sq(t+t0)))

double eta = 0.01;

double bell (double x, double y) {
  return Gaussian;
}

void get_errors (double * emp, double * etp) {
  double em = 0, et = 0;
  foreach (reduction (max:em) reduction (+:et)) {
    double ddt = 0;
    foreach_dimension()
      ddt += (s[-1] - 2*s[0] + s[1])/sq(Delta);
    double err = fabs(ddt - Gaussiand);
    et += sq(Delta)*err;
    if (err > em)
      em = err;
    }
  emp[0] = em;
  etp[0] = et;
}

int main() {
  L0 = 20;
  X0 = Y0 = -L0/2.;
  // N x N grid
  for (N = 8; N <= 512; N *= 2) {
    init_grid (N);
    foreach() 
      s[] = Gauss6 (x, y, Delta, bell);
    double em, et;
    get_errors (&em, &et);
    fprintf (stderr, "%g %g %g\n", sqrt(grid->tn), et, em);
  }
  // Using adapt_wavelet
  for (eta = 0.01; eta > 1e-5; eta /= 1.25) {
    init_grid(8);
    s.prolongation = refine_bilinear;
    s.dirty = true;
    do {
      foreach() 
	s[] = Gauss6 (x, y, Delta, bell);
    } while (adapt_wavelet ({s}, (double[]){eta}, 99).nf);
    s.prolongation = refine_4th;
    s.dirty = true;
    boundary ({all});
    double em, et;
    get_errors (&em, &et);
    printf ("%g %g %g %d\n", sqrt(grid->tn), et, em, depth());
    view (width = 1024, height = 1024, fov = 19);
    squares ("s");
    cells();
    char fname[99];
    sprintf (fname, "s%d.png", (int)grid->tn);
    save (fname);
  }
}
