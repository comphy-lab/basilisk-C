
/**
# Higher-order accurate face-value prolongation  

 Here we test the methods for 4th order face refinement and Solenoidal
 refinement.
 */
#include "higher-order.h"
#include "poisson.h" //For projection()
#include "utils.h"   

face vector fv[];
scalar s[];

double max_divergence (face vector fv);//Prototype

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2;
  foreach_dimension() {
    periodic (left);
    fv.x.prolongation = refine_face_4_x;
  }
  fv.x.refine = refine_face_solenoidal;
  
  /**
     First, the prolongation method is tested:
     
     ~~~gnuplot 4th order prolongation
     set logscale xy 2
     set xr [9:4000]
     set yr [1e-12:1e-2]
     set grid
     set xlabel 'N'
     set ylabel 'L_1 Error'
     set size square
     plot 'out' u 1:2, 1000*x**(-4)
     ~~~
  */
  for (int l = 4; l < 12; l++) {
    init_grid (1 << l);
    foreach_face() {
      coord f = {y, -x};
      fv.x[] =  f.x*exp(-sq(x) - sq(y));
    }
    boundary ((scalar*){fv});
    wavelet (fv.x, s);
    fprintf (stdout, "%d %g\n", N, normf(s).avg); 
  }
  /**
     Next, the Solenoidal refinement test is performed
     
     Accuracy:
     
     ~~~gnuplot 5th order solenoidal refinement
     set logscale xy 2
     set xr [9:4000]
     set yr [1e-13:1e-2]
     set grid
     set size square
     set xlabel 'N'
     set ylabel 'L_1 Error'
     plot 'log' u 1:2, 500000*x**(-5)
     ~~~
     
     Divergence:

     ~~~gnuplot The divergence is within the machine precision
     reset
     set logscale x
     set xr [9:4000]
     set yr [-1e-11:1e-11]
     set grid
     set xlabel 'N'
     set ylabel 'Max. absolute divergence'
     set size square
     plot 'log' u 1:4
     ~~~
     
     
  */
  NITERMIN = 25;
  for (int l = 4; l < 12; l++) {
    init_grid (1 << l);
    foreach_face() {
      coord f = {y, -x};
      fv.x[] =  f.x*exp(-sq(x) - sq(y));
    }
    project (fv, s);
    double md = max_divergence (fv); //approx zero
    double * fva = malloc (2*(sq(N))*sizeof(double));
    long int j = 0;
    foreach_face()
      fva[j++] = fv.x[]; //Array of reference values.
    
    unrefine (level == l - 1);
    refine (level < l);

    double err = 0;
    j = 0;
    foreach_face()
      err += fabs(fv.x[] - fva[j++]);
    fprintf (stderr, "%d %g %g %g\n",
	     N, err/sq(N), md, max_divergence (fv));
    free (fva);
  }
}

double max_divergence (face vector fv) {
  double max = 0; 
  foreach() {
    double d = 0;
    foreach_dimension()
      d += fv.x[1] - fv.x[];
    d /= Delta;
    if (fabs(d) > max)
      max = fabs(d);
  }
  return max;
}
