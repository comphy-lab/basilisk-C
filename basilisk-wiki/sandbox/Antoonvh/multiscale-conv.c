/**
![An approximate solution to the multiscale Poisson problem studied herein](multiscale-conv/a.png)

## Convergence for the gradient of a solution to a Poisson problem

Lets look at the convergence using various refinement methods.

~~~gnuplot $L_1$ error norm
set size square
set xr [28:1224]
set yr [0.004:10]
set logscale x 2
set logscale y
set xlabel 'N'
set ylabel 'Error'
set key right outside 
set grid
a = 5
col_equi = '#000000'
col_wav = '#AA01AA'
col_met = '#01AAAA'
plot 'err_equi' u 1:a lt rgb col_equi t 'Equidistant',\
for [i = 9:12] 'err_wav'.i u 1:a lt rgb col_wav t 'Wavlet ML = '.i ,\
for [c = 0:3] i=2**c 'err_met'.i u 1:a lt rgb col_met t 'metric norm = '.i
~~~

~~~gnuplot Max error 
set size square
set xr [28:1224]
set yr [0.001:4]
set logscale x 2
set logscale y
set xlabel 'N'
set ylabel 'Error'
set key right outside 
set grid
a = 6
col_equi = '#000000'
col_wav = '#AA01AA'
col_met = '#01AAAA'
plot 'err_equi' u 1:a lt rgb col_equi t 'Equidistant',\
for [i = 9:12] 'err_wav'.i u 1:a lt rgb col_wav t 'Wavlet ML = '.i ,\
for [c = 0:3] i=2**c 'err_met'.i u 1:a lt rgb col_met t 'metric norm = '.i
~~~

It seems that the metric-based adaption is a solid performer. Can it be improved?
*/

#include "higher-order.h"
#include "poisson.h"
#include "../prouvost/AMR_tools/amr.h"
#include "utils.h"

#define GN (3)

double sigmas[GN] = {.2, .5, 1};
double As[GN] = {1,2,4};
double xc, yc, A, sigma;
double source_2D (double x, double y) {
  return A*4*(sq(x - xc) + sq(y - yc) - sq(sigma))*exp((-sq(x - xc) - sq(y - yc))/sq(sigma))/(sq(sq(sigma)));
}

double solution (double x, double y) {
  return A*exp((-sq(x - xc) - sq(y - yc))/sq(sigma));
}

double grad_x (double x, double y) {
  return -2*A*(x - xc)/sq(sigma)*exp((-sq(x - xc) - sq(y - yc))/sq(sigma));
}

double grad_y (double x, double y) {
  return -2*A*(y - yc)/sq(sigma)*exp((-sq(x - xc) - sq(y - yc))/sq(sigma));
}

scalar a[], b[];

void init_problem (void) {
  foreach() {
    a[] = b[] = 0;
    for (int i = 0; i < GN; i++)
      for (int j = 0; j < GN; j++) {
	xc = X0 + (i + 0.75)*L0/(GN + .5);
	yc = Y0 + (j + 0.75)*L0/(GN + .5);
	sigma = sigmas[i];
	A = As[j];
	b[] += Gauss6 (x, y, Delta, source_2D);
      }
  }
}

void print_errors (FILE * fp) {
  scalar sol[];
  face vector sol_grad[];
  foreach() 
    sol[] = 0;
  foreach_face()
    sol_grad.x[] = 0;
  for (int i = 0; i < GN; i++)
    for (int j = 0; j < GN; j++) {
      xc = X0 + (i + 0.75)*L0/(GN + .5);
      yc = Y0 + (j + 0.75)*L0/(GN + .5);
      sigma = sigmas[i];
      A = As[j];
      foreach() 
	sol[] += Gauss6 (x, y, Delta, solution);
      foreach_face()
	sol_grad.x[] += Gauss6_x(x, y, Delta, grad_x);
    }
  double err = 0, me = 0, err_grad = 0, merr_grad = 0;
  foreach(reduction(+:err), reduction(max:me)) {
    double le = fabs (a[] - sol[]);
    err += sq(Delta)*le;
    if (le > me)
      me = le;
  }
  foreach_face(reduction(+:err_grad), reduction(max:merr_grad)) {
    double le = fabs(face_gradient_x(a, 0) - sol_grad.x[]);
    err_grad += sq(Delta)*le;
    if (le > merr_grad)
      merr_grad = le;
  }
  fprintf (fp, "%g %d %g %g %g %g\n", sqrt(grid->tn), grid->maxdepth,
	   err, me, err_grad, merr_grad);
  printf ("%g %d %g %g %g %g\n", sqrt(grid->tn), grid->maxdepth,
	  err, me, err_grad, merr_grad);
}

int main() {
  TOLERANCE = 1e-5;
  L0 = 30;
  foreach_dimension() {
    a[left] = dirichlet (0);
    a[right] = dirichlet (0);
  }
  /**
First, we study the convergence on increasingly refined equidistant meshes
   */
  FILE * fp_equi = fopen("err_equi", "w");
  for (maxlevel = 8; maxlevel <= 10; maxlevel++) {
    init_grid (1 << maxlevel);
    init_problem();
    poisson (a, b);
    print_errors (fp_equi);
  }
  fclose (fp_equi);
  output_ppm (a, file = "a.png", n = 600, min = -2, max = 2);
  /**
Next, we study the convergence for the wavelet-based adaptation with
increasingly strict refinement criteria on grids constrained to
varying maximum depths.
   */
  for (maxlevel = 9; maxlevel <= 12; maxlevel++) {
    char fname[99];
    sprintf (fname, "err_wav%d", maxlevel);
    FILE * fp_wav = fopen (fname, "w");
    init_grid (32);
    for (double eps = 5e-2; eps > 1e-4; eps /= 1.5) {
      // Find mesh
      do {
	init_problem();
	a.prolongation = refine_3rd;//improves solution
	poisson (a, b);
	a.prolongation = refine_bilinear; //refinement indicator
      } while (adapt_wavelet ({a}, (double[]){eps}, maxlevel).nf > grid->tn/50);
      a.prolongation = refine_3rd;
      init_problem();
      poisson (a, b);
      print_errors (fp_wav);
    }
    fclose (fp_wav);
  }
  /**
The metric-based refinement functions by Lucas Prouvost are convinient
to adopt. We study the convergence with decreasing `AMReps` values for
various norms.  */
  for (user_norm = 1; user_norm <= 9; user_norm *= 2) {
    init_grid (64);
    char fname[99];
    sprintf (fname, "err_met%d", user_norm);
    FILE * fp_met = fopen (fname, "w");
    for (AMReps = 1e-2; AMReps > 1e-6*pow(1.5, user_norm); AMReps /= 2) {
      //find mesh
      do {
	init_problem();
	a.prolongation = refine_3rd;//improves solution
	poisson (a, b);
	a.prolongation = refine_bilinear; //refinement indicator?
      } while (adapt_metric ({a}).nf > grid->tn/50);
      a.prolongation = refine_3rd;
      init_problem();
      poisson (a, b);
      print_errors (fp_met);
    }
    fclose (fp_met);
  }
}
