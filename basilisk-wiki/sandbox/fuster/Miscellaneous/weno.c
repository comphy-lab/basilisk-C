/**
# Discrete derivative with 3rd- and 5th-order WENO

*/

// we need two layers of ghost cells for the 5-points stencil near boundaries
#define BGHOSTS 2

#include "grid/multigrid1D.h"
#include "utils.h"

double sigma = 0.1;

// Note dir is 1 or -1 depending on the direction

foreach_dimension()
static double weno3_x (Point point, scalar f, int dir)
{
  static double epsilon = 1.e-7;
  static double coeff[2][2] = {
    {-1./2., 3./2.},
    {1./2., 1./2.}
  };    
  static double weights[2] = {1./3.,2./3.};

  double p[2], beta[2], ar[2];

  beta[0] = sq(f[] - f[-dir]);
  beta[1] = sq(f[dir]  - f[]);
  
  p[0] = coeff[0][0]*f[-dir] + coeff[0][1]*f[];
  p[1] = coeff[1][0]*f[]     + coeff[1][1]*f[dir];
  
  double arsum = 0.;
  for (int i = 0; i < 2; i++) {
    ar[i] = weights[i]/sq(epsilon + beta[i]);
    arsum += ar[i];
  }
  
  double fweno = 0.;
  for (int i = 0; i < 2; i++)
    fweno += ar[i]/arsum*p[i];    

  return fweno;
}


// Note dir is 1 or -1 depending on the direction

foreach_dimension()
static double weno5_x (Point point, scalar f, int dir) 
{
  static double epsilon = 1.e-7;  
  static double coeff[3][3]= {
    {-1./6., 5./6., 1./3.},
    {1./3., 5./6., -1./6.},
    {11./6., -7./6., 1./3.}
  };
  static double weights[3] = {3./10.,3./5.,1./10.};

  double p[3], beta[3], ar[3];
  
  p[0] = coeff[0][0]*f[2*dir] + coeff[0][1]*f[dir]  + coeff[0][2]*f[];
  p[1] = coeff[1][0]*f[dir]   + coeff[1][1]*f[]     + coeff[1][2]*f[-dir];
  p[2] = coeff[2][0]*f[]      + coeff[2][1]*f[-dir] + coeff[2][2]*f[-2*dir];

  beta[0] = 13./12.*sq(f[2*dir] - 2.*f[dir] + f[]) + sq(f[2*dir] - 4.*f[dir] + 3.*f[])/4.;
  beta[1] = 13./12.*sq(f[1] - 2.*f[] + f[-1]) + sq(f[-1] - f[1])/4.;
  beta[2] = 13./12.*sq(f[] - 2.*f[-dir] + f[-2*dir]) + sq(3.*f[] - 4.*f[-dir] + f[-2*dir])/4.;
  
  double arsum = 0.;
  for (int i = 0; i < 3; i++) {
    ar[i] = weights[i]/sq(epsilon + beta[i]);
    arsum += ar[i];
  }
  
  double fweno = 0;
  for (int i = 0; i < 3; i++)
    fweno += ar[i]/arsum*p[i];

  return fweno;
}

void get_weno_derivative (scalar f, int order, scalar df)
{
  double (* weno)(Point, scalar, int) =
    order == 3 ? weno3_x :
    order == 5 ? weno5_x :
    NULL;
  foreach()
    df[] = (weno (point, f, 1) - weno (point, f, -1))/Delta;
  boundary ({df});
}

int main() {
  double x0 = 0.5;
  
  /**
  We replace the constant value of *n* with a loop going through
  16, 32, 64, 128 and 256. */

  for (int n = 16; n <= 256; n *= 2) {
    init_grid (n);
    
    scalar f[], dfexact[];
    foreach() {
      f[] = erf((x-x0)/sigma);
      dfexact[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-x0)/sigma,2));
    }
    boundary ({f, dfexact});
    
    scalar dw5[], dw3[];
    get_weno_derivative (f, 3, dw3);
    get_weno_derivative (f, 5, dw5);
    
    if (n == 32) {
      FILE * fp = fopen("solution.dat","w");
      foreach()
        fprintf (fp, "%g %g %g %g \n", x, dw5[], dw3[], dfexact[]);
      fclose(fp);
    }
      
    scalar e3[], e5[];
    foreach() {
      e5[] = dw5[] - dfexact[];
      e3[] = dw3[] - dfexact[];
    }
    
    printf ("%d %g %g\n", n, normf(e5).max, normf(e3).max);
    
    free_grid();
  }
}

/**

~~~gnuplot Numerical derivative approximation
set xlabel 'x'
set ylabel 'df/dx'
sigma=0.1
x0=0.5
f(x)=2./sigma/sqrt(pi)*exp(-((x-x0)/sigma)**2)
p "solution.dat" u 1:2 t 'Weno 5th' w p, \
  "solution.dat" u 1:3 t 'Weno 3rd' w p, \
  f(x) t 'Exact' w l
~~~ 

~~~gnuplot Convergence: The derivative converges with one order less than the approximation of the face values
set logscale
set xlabel '1/dx'
set ylabel 'normf(e).max'
fit [3:]a*x+b 'out' u (log($1)):(log($2)) via a,b
fit [3:]a1*x+b1 'out' u (log($1)):(log($3)) via a1,b1
set xtics 8,2,512
set cbrange [1:2]
plot [8:512]'out' u 1:2 pt 7 t '5th-order WENO', \
            exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a), \
            'out' u 1:3 pt 7 t '3rd-order WENO', \
            exp(b1)*x**a1 t sprintf("%.0f/n^{%4.2f}", exp(b1), -a1)
~~~ 

*/  
