/**
# Discrete derivative 3rd-, 4rd and 5th-order 
*/

// we need two layers of ghost cells for the 5-points stencil near boundaries
#define BGHOSTS 2

#include "grid/multigrid.h"
#include "utils.h"

typedef struct { 
  double * data;
  int r[3];
  int np[3];
} MData;

/* Test function */
double sigma = 0.1;
double x0 = 0.5;

double fexact (double x) {

   return erf((x-x0)/sigma);

}

// algorithm taken from http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
void add_combinations (int v[], double data[], int start, int n, int k, int maxk, double *sum) {

  int i;
  if (k > maxk) {
    double f = 1.;
    for (i=1; i<=maxk; i++) 
       f *= data[v[i]-1];
    *sum += f;
    return;
  }

  for (i=start; i<=n; i++) {
    v[k] = i;
    add_combinations (v, data, i+1, n, k+1, maxk, sum);
  }
}

// Obtains the analytical integral of (x-x0)*..(x-xN) in the range [x0:x1]
double analytic_integral( double vals[], int porder, double x0, double x1 ) {

  int v[100];

  double integral = 0;
  for (int i=0; i <= porder-1; i++) {
    double coef = 0.;
    add_combinations (v, vals, 1, porder, 1,porder-i, &coef);
    for (int j=i; j > i; j--)
        coef *= j;
    coef *= pow(-1,porder-i);
    integral += coef*(pow(x1, i+1) - pow(x0, i+1))/(i+1);
  }
  integral += (pow(x1, porder+1) - pow(x0, porder+1))/(porder+1);
    
  return integral;
}


double eval_poly ( double x0[], int np, double x ) {

    double Li = 1.;
    for (int i=0; i< np; i++) 
        Li *= (x - x0[i]);
    return Li;
}

double * eval_legendre1D (double x0[], int np, double x) {

    double * coefs = (double *) malloc(np * sizeof(double));
    for (int i=0; i< np; i++) {
        double xi = x0[i];
        x0[i] = x0[np-1];
        double denom = eval_poly(x0, np-1, xi);
        double num   = eval_poly(x0, np-1, x);
        x0[i] = xi;
        coefs[i] = num/denom;
    }
    return coefs;

}

double * integrate_legendre1D (double x0[], int np, double x1, double x2) {

    double * coefs  = (double *) malloc(np * sizeof(double));
    for (int i=0; i< np; i++) {
        double xi = x0[i];
        x0[i] = x0[np-1];
        double denom = eval_poly(x0, np-1, xi);
        double num = analytic_integral( x0, np-1, x1, x2 );
        x0[i] = xi;
        coefs[i] = num/denom;
    }
    return coefs;
}

/* Integrate the legendre polynomial */
MData Init_Matrix_Integration( double x1, double x2, int np, int r) {

  MData M;
  for (int d=0; d < dimension; d++){
      M.np[d] = np; 
      M.r[d] = r;
  }
  for (int d=dimension; d < 3; d++){
      M.np[d] = 1; 
      M.r[d] = 0;
  }

  M.data = (double *) malloc(M.np[0]*M.np[1]*M.np[2] * sizeof(double)); 

  double *coefs[3];
  for (int d=0; d < 3; d++){

    double * x0 = (double*) malloc(M.np[d] * sizeof(double));
    for (int i=0; i< M.np[d]; i++) 
      x0[i]=i-M.r[d];

    coefs[d] = integrate_legendre1D (x0, M.np[d], x1, x2);
    
    free(x0); 

  }

  for (int i=0; i< M.np[0]; i++) 
      for (int j=0; j< M.np[1]; j++) 
          for (int k=0; k< M.np[2]; k++) 
              M.data[(i*M.np[1] + j)*M.np[2] + k] = coefs[0][i]*coefs[1][j]*coefs[2][k];

  for (int d=0; d < 3; d++)
      free(coefs[d]);

  return M;
}

double * derive_legendre1D (double x0[], int np, double x) {

    double *coefs = (double *) malloc(np * sizeof(double));
    for (int i=0; i< np; i++) {
        double xi1 = x0[i];
        x0[i] = x0[np-1];
        double denom = eval_poly(x0, np-1, xi1);
        double num = 0.;
        for (int j=0; j< np-1; j++) {
            double xi2 = x0[j];
            x0[j] = x0[np-2];
            num += eval_poly(x0, np-2, x);
            x0[j] = xi2;
        }
        x0[i] = xi1;
        coefs[i] = num/denom;
    }
    return coefs;
}

/* Derivate of the legendre polynomial */
MData Init_Matrix_Derivate( double x[], int np, int r ) {

  MData M;
  for (int d=0; d < dimension; d++){
      M.np[d] = np; 
      M.r[d] = r;
  }
  for (int d=dimension; d < 3; d++){
      M.np[d] = 1; 
      M.r[d] = 0;
  }

  M.data = (double *) malloc(M.np[0]*M.np[1]*M.np[2] * sizeof(double)); 

  double * x0 = (double*) malloc(M.np[0] * sizeof(double));

  for (int i=0; i< M.np[0]; i++) 
      x0[i]=i-M.r[0];

  double * coefs[3];
  coefs[0] = derive_legendre1D (x0, M.np[0], x[0]);
  coefs[1] = eval_legendre1D (x0, M.np[1], x[1]); 
  coefs[2] = eval_legendre1D (x0, M.np[2], x[2]); 

  free(x0); 

  for (int i=0; i< M.np[0]; i++) 
      for (int j=0; j< M.np[1]; j++) 
          for (int k=0; k< M.np[2]; k++) 
              M.data[(i*M.np[1] + j)*M.np[2] + k] = coefs[0][i]*coefs[1][j]*coefs[2][k];

  for (int d=0; d < 3; d++)
      free(coefs[d]);

  return M;
}

foreach_dimension()
static double eval_x (Point point, scalar f, MData * M) {

    double sum = 0.;
    for (int i=0; i< M->np[0]; i++) 
        for (int j=0; j< M->np[1]; j++) 
            for (int k=0; k< M->np[2]; k++) 
                sum += M->data[(i*M->np[1] + j)*M->np[2] + k]*f[-M->r[0] + i,-M->r[1] + j, -M->r[2] + k];

    return sum;
}


// This function gives the cell centered averaged value
// given a field variable using a stencil of length order/2
void Init_avg_variable ( scalar f, scalar favg, int order ) 
{
  
  MData LM = Init_Matrix_Integration( -1./2, 1./2, order, order/2 );

  foreach() 
    favg[] = eval_x(point, f, &LM);

  boundary({favg});

  free(LM.data);
}

void Deriv_variable (scalar f, vector df, MData *M)
{
    foreach() 
        foreach_dimension () 
            df.x[] = eval_x (point, f, M)/Delta;
    
    boundary ((scalar *){df});
}

int main() {

  double xder[3] = {0.,0.,0.};
  MData MDer3 = Init_Matrix_Derivate( xder, 3, 1 );
  MData MDer5 = Init_Matrix_Derivate( xder, 5, 2 );
  xder[0] = 1./2;
  MData MDer4 = Init_Matrix_Derivate( xder, 4, 1 );

  for (int n = 16; n <= 256; n *= 2) {
  init_grid (n);
  
  scalar fe[];
  vector dfe[], dfef[];
  foreach() {
    fe[] = fexact(x);
    dfe.x[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-x0)/sigma,2));
    dfe.y[] = 0.;
    dfef.x[] = 2./sigma/sqrt(M_PI)*exp(-pow((x + Delta/2. -x0)/sigma,2));
    dfef.y[] = 0.;
  }
  boundary({fe});
  boundary((scalar *){dfe, dfef});
  
  scalar favg3[], favg4[], favg5[], dfe3x[], dfe4x[], dfe5x[];
  Init_avg_variable ( fe, favg3, 3 );
  Init_avg_variable ( fe, favg4, 4 );
  Init_avg_variable ( fe, favg5, 5 );
  Init_avg_variable ( dfe.x, dfe3x, 3 );
  Init_avg_variable ( dfef.x, dfe4x, 4 );
  Init_avg_variable ( dfe.x, dfe5x, 5 );

  vector dw5[], dw4[], dw3[];
  Deriv_variable (favg3, dw3, &MDer3);
  Deriv_variable (favg4, dw4, &MDer4);
  Deriv_variable (favg5, dw5, &MDer5);
  
  if (n == 32) {
      FILE * fp = fopen("solution.dat","w");
      foreach()
          fprintf (fp, "%g %g %g %g %g \n", x, x+Delta/2., dw5.x[], dw4.x[], dw3.x[]);
      fclose(fp);
  }

  scalar e3[], e4[], e5[];
  foreach() {
      e5[] = dw5.x[] - dfe5x[];
      e4[] = dw4.x[] - dfe4x[];
      e3[] = dw3.x[] - dfe3x[];
  }

  printf ("%d %g %g %g\n", n, normf(e5).max, normf(e4).max, normf(e3).max);

  free_grid();

  }

  free(MDer3.data); free(MDer4.data); free(MDer5.data);
}

/**

~~~gnuplot Numerical derivative approximation
set xlabel 'x'
set ylabel 'df/dx'
sigma=0.1
x0=0.5
f(x)=2./sigma/sqrt(pi)*exp(-((x-x0)/sigma)**2)
p "solution.dat" u 1:3 t 'deriv 5th' w p, \
  "solution.dat" u 2:4 t 'deriv 4rd' w p, \
  "solution.dat" u 1:5 t 'deriv 3rd' w p, \
  f(x) t 'Exact' w l
~~~ 

~~~gnuplot Convergence: The derivative converges with one order less than the approximation of the face values
set logscale
set xlabel '1/dx'
set ylabel 'normf(e).max'
fit [3:]a*x+b 'out' u (log($1)):(log($2)) via a,b
fit [3:]a1*x+b1 'out' u (log($1)):(log($3)) via a1,b1
fit [3:]a2*x+b2 'out' u (log($1)):(log($4)) via a2,b2
set xtics 8,2,512
set cbrange [1:2]
plot [8:512]'out' u 1:2 pt 7 t '5th-order deriv', \
            exp(b)*x**a t sprintf("%.0f/n^{%4.2f}", exp(b), -a), \
            'out' u 1:3 pt 7 t '4rd-order deriv', \
            exp(b1)*x**a1 t sprintf("%.0f/n^{%4.2f}", exp(b1), -a1), \
            'out' u 1:4 pt 7 t '3rd-order deriv', \
            exp(b2)*x**a2 t sprintf("%.0f/n^{%4.2f}", exp(b2), -a2)
~~~ 

*/
