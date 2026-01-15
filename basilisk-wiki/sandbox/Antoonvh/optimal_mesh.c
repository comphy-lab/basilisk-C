/**
# Optimal bitree meshes for the second derivative

We seek to minimize a given $L_p$ error norm of a second-order
accurate estimate for the second derivative for a given number of
cells. The $L_p$ norm for some approximated field $f$ in $D$ dimensions is given by,

$$L_p = \left(\sum_{i}^{\text{cells}} |f_{a,i} - f_{e,i}|^{p}\Delta_i^D \right)^{1/p},$$

where $f_{a,i}$ and $f_{e,i}$ are the analytical and estimated value
for the $f$ field in the $i$-th cell. To minimize $L_P$, we must
minimize the sum in the base. The contribution of the $i$-th cell is then ($\epsilon_i$)

$$\epsilon_i = |f_{a,i} - f_{e,i}|^{p}\Delta_i^D$$

For a second-order accurate second derivative in 1D,

$$\epsilon_i \propto |\frac{\mathrm d^4 f}{\mathrm d x^4}\Delta_i^2|^{p}\Delta_i.$$

This fourth derivative can be estimated using a 4-th order accurate
multilevel interpolation (wavelet: `refine_4th`). For $p \geq 0$,
$L_p$ is minimized when $\epsilon_i$ is equidistributed over the mesh.

## Results

Results for a Gaussian bump:

### Example mesh statistics

~~~gnuplot Number of cells per level for $N \approx 512$
set style data histogram
set style histogram cluster gap 2
set xlabel 'Level'
set ylabel 'Cells'
set style fill solid border -1
set key top left
set xr [2.1:11]
set yr [0:512]
plot 'level_1' u 2:xtic(1) t 'Optimal L_1', 'level_2' u 2:xtic(1) t 'Optimal L_2', 'level_4' u 2:xtic(1) t 'Optimal L_4', 'equi_level' u 2:xtic(1) t 'Equi', 'wavelet_level' u 2:xtic(1) t 'Default Wavelet'
~~~

### Errors

~~~gnuplot $L_1$ norm
reset
set size square
set logscale x 2
set logscale y
set xr [10: 1500]
set ylabel 'L_1 norm'
set grid
plot 'optimal_1', 'optimal_2', 'optimal_4', 'equi', 'default_wavelet'
~~~

~~~gnuplot $L_2$ norm
set size square
set logscale x 2
set logscale y
set xr [10: 1500]
set grid
set xlabel 'Cells'
set ylabel 'L_2 norm'
set key bottom left
plot 'optimal_1' u 1:3, 'optimal_2' u 1:3, 'optimal_4' u 1:3, 'equi' u 1:3, 'default_wavelet' u 1:3
~~~

~~~gnuplot $L_4$ norm
set size square
set logscale x 2
set logscale y
set xr [10: 1500]
set grid
set xlabel 'Cells'
set ylabel 'L_4 norm'
set key bottom left
plot 'optimal_1' u 1:4, 'optimal_2' u 1:4, 'optimal_4' u 1:4, 'equi' u 1:4, 'default_wavelet' u 1:4
~~~

The optimal meshes all seem to perform well, Even for the Error norms
that they are not concererd with. They are all much better
than the default wavelet-based approach.
 */
#include "grid/bitree.h"
vector u; // :S
#include "adapt_field.h"
#include "higher-order.h"

double Gaussian (double x) {
  return (exp(-sq(x)));
}

double flux_div (double x, double Delta) {
  double xl = x - Delta/2.;
  double xr = x + Delta/2.;
  return (-2*xr*exp(-sq(xr)) + 2*xl*exp(-sq(xl)))/Delta;
}

double diagnose_error (scalar s, int norm) {
  // Give all methods accurate ghost cells
  s.prolongation = refine_5th;  // This does not do much
  tree->dirty = true;  // This does not do much
  double err = 0;
  foreach(reduction(+:err)) {
    double errl = flux_div(x, Delta) - (s[1] - 2*s[] + s[-1])/sq(Delta);
    err += Delta*pow(fabs(errl), norm);
  }
  return pow(err, 1./(double)norm);
}

int main() {
  L0 = 10;
  X0 = -L0/2.;
  int norms[3] = {1, 2, 4};
 
  // Optimal?
  for (int j = 0 ; j < 3; j++) {
    int norm = norms[j];
    int expo = -2 + -4*norm;
    double zmin = pow(10, expo);
    char fname[99];
    sprintf (fname,"optimal_%d", norm);
    FILE * fp = fopen(fname, "w");
    for (double zeta = 1e-2; zeta > zmin; zeta /= 3) {
      init_grid(1024);
      periodic(left);
      scalar s[], w[];
      s.prolongation = refine_4th;
      astats r = {0};
      do {
	foreach()
	  s[] = Gauss6 (x, Delta, Gaussian);
	wavelet (s, w);
	foreach_cell() 
	  w[] = pow(fabs(w[])/sq(Delta), (double)norm)*pow(Delta, dimension);
	r = adapt_field (w, zeta, zeta/5., 99);
      } while (r.nf || r.nc);
      fprintf (fp,"%ld %g %g %g\n", grid->tn, diagnose_error(s, norms[0]),
	       diagnose_error(s, norms[1]), diagnose_error(s, norms[2]));
      // Output grid
      if (grid->tn <= 512 && grid->tn > 350) {
	sprintf (fname, "level_%d", norm);
	FILE * fp1 = fopen (fname, "w");
	int cells[99] = {0};
	foreach(reduction(+:cells[:grid->depth+1]))
	  cells[level]++;
	for (int i = 1; i <= grid->depth; i++) 
	  fprintf(fp1, "%d %d\n", i, cells[i]);
	fclose (fp1);
      }
    }
    fclose (fp);
   
  }
  // wavelet
  FILE * fp = fopen("default_wavelet","w");
  for (double zeta = 1e-2; zeta > 1e-5; zeta /= 2) {
    init_grid(1024);
    periodic(left);
    scalar s[];
    astats r = {0};
    do {
      foreach()
	s[] = Gauss6 (x, Delta, Gaussian);
      boundary ({s});
      r = adapt_wavelet ({s}, (double[]){zeta}, 99);
    } while (r.nf || r.nc);
    fprintf (fp, "%ld %g %g %g\n", grid->tn, diagnose_error(s, norms[0]),
	     diagnose_error(s, norms[1]), diagnose_error(s, norms[2]));
    // Output grid
      if (grid->tn <= 550 && grid->tn > 350) {
	FILE * fp1 = fopen ("wavelet_level", "w");
	int cells[99] = {0};
	foreach(reduction(+:cells[:grid->depth+1]))
	  cells[level]++;
	for (int i = 1; i <= grid->depth; i++) 
	  fprintf(fp1, "%d %d\n", i, cells[i]);
	fclose (fp1);
      }
  }
  fclose(fp);
   // equidistant
  fp = fopen("equi", "w");
  for (N = 32; N <= 1024; N *= 2) {
    init_grid(N);
    periodic(left);
    scalar s[];
    foreach()
	s[] = Gauss6 (x, Delta, Gaussian);
    fprintf (fp, "%ld %g %g %g %g\n", grid->tn, diagnose_error(s, norms[0]),
	     diagnose_error(s, norms[1]), diagnose_error(s, norms[2]));
     // Output grid
      if (grid->tn == 512) {
	FILE * fp1 = fopen ("equi_level", "w");
	int cells[99] = {0};
	foreach(reduction(+:cells[:grid->depth+1]))
	  cells[level]++;
	for (int i = 1; i <= grid->depth; i++) 
	  fprintf(fp1, "%d %d\n", i, cells[i]);
	fclose (fp1);
      }
  }
  
}
