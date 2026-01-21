/**
# Optimal bitree meshes for an advection-diffusion equation

We seek to minimize a given $L_p$ error norm for the estimated tendency in the equation:
$$\frac{\partial s}{\partial t}= -U\frac{\partial s}{\partial x} + \kappa\frac{\partial^2 s}{\partial x^2}$$
Using second order accurate approximatations for the patial derivatives. 

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
that they are not concererd with. Interestingly, the wavelet-based approach seems to do perform rather well too.
 */
#include "grid/bitree.h"
vector u; // :S
#include "adapt_field.h"
#include "higher-order.h"

double U = 1;
double kap = 1;

double Gaussian (double x) {
  return (exp(-sq(x)));
}

double Gaussian_1st (double x) {
  return (-2*x*exp(-sq(x)));
}

double flux_div (double x, double Delta) {
 double xl = x - Delta/2.;
  double xr = x + Delta/2.;
  double diffl = -kap*Gaussian_1st(xl);
  double diffr = -kap*Gaussian_1st(xr);
  double advl = U*Gaussian(xl);
  double advr = U*Gaussian(xr);
  return (diffl + advl - diffr - advr)/Delta;
}

static inline void refine_Gaussian (Point point, scalar s) {
  foreach_child() {
    s[] = Gauss6(x, Delta, Gaussian);
  }
}
/**
## Note on Conservative flux restriction

The optimal meshes utilize the fact that the leading-order error term
for the fluxes (mostly) cancel out when computing the flux'
divergence. At resolution boundaries, this is not the case, as the
accurate high-level flux is restricted to the lower level cells. This
introduces an addional errors that can become dominant on our "simply
optimized" meshes. Indeed, Using the less accurate coarse-level flux
for the coarse levels yields lower errors on these "optimal"
meshes. As such, we do not use the conservative formulation, eventough we
do consider our variables to be cell-averages (for accurate
restriction).

~~~gnuplot The conservation error converges at second order
reset
set size square
set xr [10:800]
set grid
set ylabel 'Conservation error'
set xlabel 'Cells'
set logscale y
set logscale x 2
plot 'conservation_1','conservation_2','conservation_4', 10*x**(-2)
~~~
 */
double diagnose_error (int norm, double * c = NULL) {
  //  accurate ghost cells (refine_5th also works)
  scalar s[];
  s.prolongation = refine_Gaussian;  
  foreach()
    s[] = Gauss6(x, Delta, Gaussian);
  face vector flx[];
  foreach_face() 
    flx.x[] = U*(s[] + s[-1])/2. - kap*(s[] - s[-1])/(Delta); // Not to use!
  
  double err = 0;
  double conservation = 0;
  foreach(reduction(+:err) reduction(+:conservation)) {
    double tendency = (-U*(s[1] - s[-1])/(2*Delta) + kap*(s[1] - 2*s[] + s[-1])/sq(Delta));
    //double tendency = (flx.x[0] - flx.x[1])/Delta; // Not to use! 
    if (c != NULL)
      conservation += tendency * Delta;
    double errl = flux_div(x, Delta) - tendency;
    err += Delta*pow(fabs(errl), norm);
  }
  if (c != NULL)
    *c = conservation;
  return pow(err, 1./(double)norm);
}
int main() {
  L0 = 10;
  X0 = -L0/2.;
  int norms[3] = {1, 2, 4};
 
  // Optimal?
  for (int j = 0 ; j < 3; j++) {
    int norm = norms[j];
    int expo = -4 + -4*norm;
    double zmin = pow(10, expo);
    char fname[99];
    sprintf (fname,"optimal_%d", norm);
    FILE * fp = fopen(fname, "w");
    sprintf (fname,"conservation_%d", norm);
    FILE * fp3 = fopen(fname, "w");
      
    for (double zeta = 1e-2; zeta > zmin; zeta /= 3) {
      init_grid(1024);
      periodic(left);
      scalar s[], w[], eps[];
      eps.restriction = no_restriction;
      astats r = {0};
      do {
	foreach() 
	  s[] = Gauss6 (x, Delta, Gaussian);
	s.prolongation = refine_4th;
	wavelet (s, w);
	foreach_cell() 
	  eps[] = (-2./9.*kap*(w[]/sq(Delta)));
	s.prolongation = refine_3rd;
	wavelet (s, w);
	foreach_cell() {
	  eps[] += (-4./9.*-U*((child.x > 0 ? -1 : 1)*w[]/(Delta)));
	  eps[] = pow(fabs(eps[]), (double)norm)*Delta;
	}
	
	r = adapt_field (eps, zeta, zeta/5., 99);
      } while (r.nf || r.nc);
      double c = 0;
      fprintf (fp,"%ld %g %g %g\n", grid->tn, diagnose_error(norms[0], &c),
	       diagnose_error(norms[1], &c), diagnose_error(norms[2], &c));
      fprintf (fp3,"%ld %g\n", grid->tn, fabs(c));
      
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
    fclose (fp3);
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
    fprintf (fp, "%ld %g %g %g\n", grid->tn, diagnose_error(norms[0]),
	     diagnose_error(norms[1]), diagnose_error(norms[2]));
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
    fprintf (fp, "%ld %g %g %g\n", grid->tn, diagnose_error(norms[0]),
	     diagnose_error( norms[1]), diagnose_error(norms[2]));
     // Output grid
      if (grid->tn == 512) {
	FILE * fp1 = fopen ("equi_level", "w");
	int cells[99] = {0};
	foreach(reduction(+:cells[:grid->depth + 1]))
	  cells[level]++;
	for (int i = 1; i <= grid->depth; i++) 
	  fprintf(fp1, "%d %d\n", i, cells[i]);
	fclose (fp1);
      }
  }
  
}
