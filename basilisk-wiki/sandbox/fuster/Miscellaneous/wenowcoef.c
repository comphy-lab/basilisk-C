/**
# Discrete derivative with 3rd- and 5th-order WENO
*/

// we need two layers of ghost cells for the 5-points stencil near boundaries
#define BGHOSTS 2

#include "grid/multigrid1D.h"
#include "utils.h"

typedef struct { 
  int k;
  int np;
  double xi;
  double ***p;
  double *weights;
  double ***beta;
} Mweno;

double sigma = 0.1;
double x0 = 0.5;

double fexact (double x) {

   return erf((x-x0)/sigma);

}

void get_matrix_coef (double **p, int k, int r) {

  for (int i=0; i< k; i++) 
     for (int j=0; j< k; j++)  
         p[i][j] = (pow((double) -r+i+1./2,j+1)-pow((double) -r+i-1./2,j+1))/(j+1);

}

void get_matrix_coef_exact (double **p, int k, int r) {

  for (int i=0; i< k; i++) 
     for (int j=0; j< k; j++) 
         p[i][j] = pow((double) -r+i,j);

}

double get_pointval_coef ( double **M, double x, int c, int k ) {

    double eval = 0.;
    for (int i=0; i< k; i++) 
        eval += M[i][c]*pow(x,i);
    return eval;

}

void reset_matrix_coef ( double **M, int k ) {

    for (int i=0; i< k; i++) 
      for (int j=0; j< k; j++) 
        M[i][j]=0.;

}


void derive_matrix_coef ( double **M, double **Md, int nder, int k ) {

    for (int i=0; i< k-nder; i++) 
      for (int j=0; j< k; j++) {
        double c = 1.;
        for (int n=0; n< nder; n++)
          c *= i+nder-n;
        Md[i][j]=c*M[i+nder][j];
      }

    for (int i=k-nder; i< k; i++) 
      for (int j=0; j< k; j++) 
          Md[i][j] = 0.;

}

void integrate_matrix_coef ( double **M, double **MI, double x1, double x2, int k ) {

    for (int i=0; i< k; i++) 
      for (int j=0; j< k; j++) 
        MI[i][j]=M[i][j]*(pow(x2,i+1)-pow(x1,i+1))/(i+1);

}

// This function gives the "order" approximation of the 
// cell centered averaged value given an exact field variable field
void init_avg_variable ( scalar f, scalar favg, int order ) 
{
 double **P =  matrix_new (order, order, sizeof(double));
 get_matrix_coef_exact (P, order, (order-1)/2);
 matrix_inverse (P, order, 1e-10);
 integrate_matrix_coef(P, P, -1./2, 1./2, order);

 foreach() {
     favg[] = 0.;
     for (int i=0; i< order; i++) 
         for (int j=0; j< order; j++) 
             favg[] += P[i][j]*f[-(order-1)/2+j];
 }

 matrix_free (P); 
 boundary({favg});

}

void destroy_Mweno ( Mweno *Mw ) {
 
 for (int r=0; r < Mw->np; r++) {
   matrix_free (Mw->p[r]); 
   matrix_free (Mw->beta[r]); 
 }
 free(Mw->p); free(Mw->weights); free(Mw->beta); 
}

void init_weno_polynomials ( Mweno * Mw) {

  /* Fit polynomials p_r(x) */
  for (int r=0; r < Mw->np; r++) { 
    Mw->p[r] =  matrix_new (Mw->k, Mw->k, sizeof(double));
    get_matrix_coef (Mw->p[r], Mw->k, r);
    matrix_inverse (Mw->p[r], Mw->k, 1e-10);
  }

}
 
void get_weno_weights ( Mweno * Mw, double xi) {

 Mw->xi = xi;
 int order = Mw->np + Mw->k - 1; 

 /* Fit high order polynomial */
 double **P =  matrix_new (order, order, sizeof(double));
 get_matrix_coef (P, order, (order-1)/2);
 matrix_inverse (P, order, 1e-10);

 /* Obtain weights for each polynomial*/
 // code least-squares method?
  double cfirst = get_pointval_coef(Mw->p[Mw->np-1], xi, 0, Mw->k);
  double cfull = get_pointval_coef(P, xi, 0, order);
  Mw->weights[Mw->np-1] = cfull/cfirst;

  if ( order == 3 ) 
    Mw->weights[0] = 1. - Mw->weights[1];
  else if ( order == 4 ) {
    double cpol0 = get_pointval_coef(Mw->p[0], xi, 1, Mw->k);
    cfull = get_pointval_coef(P, xi, 3, order);
    Mw->weights[0] = cfull/cpol0;
  }
  else if ( order == 5 ) {
    //last point
    double cpol0 = get_pointval_coef(Mw->p[0], xi, 2, Mw->k);
    cfull = get_pointval_coef(P, xi, 4, order);
    Mw->weights[0] = cfull/cpol0;
    
    //complementary weight
    Mw->weights[1] = 1. - Mw->weights[0] - Mw->weights[2];
  }

 matrix_free(P);

}

/* Obtain smooth indicators (defined for each polynomial irrespective of xi!) */
void get_weno_beta ( Mweno * Mw) {

 double **Dp     =  matrix_new (Mw->k, Mw->k, sizeof(double));
 for (int r=0; r < Mw->np; r++) {
    Mw->beta[r]    =  matrix_new (Mw->k, Mw->k, sizeof(double));
    reset_matrix_coef ( Mw->beta[r], Mw->k );
    for (int l=1; l < Mw->k-1; l++) {
        derive_matrix_coef ( Mw->p[r], Dp, l, Mw->k );
        integrate_matrix_coef ( Dp, Dp, -1./2, 1./2, Mw->k );
        for (int i=0; i< Mw->k; i++) 
          for (int j=0; j< Mw->k; j++) 
            Mw->beta[r][i][j] += Dp[i][j];
    }
 }
 matrix_free (Dp); 

}

void get_Mweno ( double xi, Mweno *Mw ) {

  init_weno_polynomials(Mw);
  get_weno_weights(Mw, xi);
  get_weno_beta(Mw);

}

void init_Mweno ( Mweno *Mw, int k, int order, double xi ) {
 
  Mw->k = k;
  Mw->np = order - k + 1;
  Mw->p = malloc(Mw->np * sizeof(double **));
  Mw->weights =  (double*) malloc(Mw->np * sizeof(double)); 
  Mw->beta = malloc(Mw->np * sizeof(double **));

  get_Mweno ( xi, Mw);
}

static double weno_x (Point point, scalar f, double xi, Mweno *Mw)
{

 static double epsilon = 1.e-7;

 double * beta = (double*) malloc(Mw->k * sizeof(double));

 for (int r=0; r < Mw->np; r++) {
     beta[r] = 0.;
     for (int i=0; i< Mw->k; i++) 
         for (int j=0; j< Mw->k; j++) 
             beta[r] += Mw->beta[r][i][j]*f[-j+i]*pow(xi,i);
     beta[r] = fabs(beta[r]);
 }

 double asum = 0.;

 if (xi/Mw->xi == 1) {
     for (int r=0; r < Mw->np; r++) {
         beta[r] = Mw->weights[r]/sq(epsilon + beta[r]);
         asum += beta[r];
     }
 } 
 else {
     for (int r=0; r < Mw->np; r++) {
         beta[r] = Mw->weights[Mw->np-1-r]/sq(epsilon + beta[r]);
         asum += beta[r];
     }
 }

 double fweno = 0.;
 for (int r=0; r < Mw->np; r++) {
     double px = 0.;
     for (int i=0; i < Mw->k; i++) 
        px += get_pointval_coef(Mw->p[r], xi, i, Mw->k)*f[-r+i];
     fweno += beta[r]/asum*px; 
 }

 free(beta);

 return fweno;

}

void get_weno_derivative (scalar f, scalar df, Mweno *Mw)
{
    foreach() {
        df[] = (weno_x (point, f, 1./2, Mw) - weno_x (point, f, -1./2, Mw))/Delta;
    }
    
    boundary ({df});
}

int main() {

  Mweno Mw3, Mw5;
  init_Mweno ( &Mw3, 2, 3, 1./2);
  init_Mweno ( &Mw5, 3, 5, 1./2);

  for (int n = 16; n <= 256; n *= 2) {
  init_grid (n);
  
  scalar fe[], dfe[];
  foreach() {
    fe[] = fexact(x);
    dfe[] = 2./sigma/sqrt(M_PI)*exp(-pow((x-x0)/sigma,2));
  }
  boundary({fe, dfe});
  
  scalar favg3[], favg5[], dfe3[], dfe5[];
  init_avg_variable ( fe, favg3, 3 );
  init_avg_variable ( fe, favg5, 5 );
  init_avg_variable ( dfe, dfe3, 3 );
  init_avg_variable ( dfe, dfe5, 5 );
  
  scalar dw5[], dw3[];
  get_weno_derivative (favg3, dw3, &Mw3);
  get_weno_derivative (favg5, dw5, &Mw5);

  
  if (n == 32) {
      FILE * fp = fopen("solution.dat","w");
      foreach()
          fprintf (fp, "%g %g %g %g \n", x, dw5[], dw3[], dfe[]);
      fclose(fp);
  }

  scalar e3[], e5[];
  foreach() {
      e5[] = dw5[] - dfe5[];
      e3[] = dw3[] - dfe3[];
  }

  printf ("%d %g %g\n", n, normf(e5).max, normf(e3).max);
  
  //foreach() 
  //  printf("%g %g %g %g \n", x, fe[], favg1[], favg2[]);

  free_grid();

  }

  destroy_Mweno ( &Mw3 ); destroy_Mweno ( &Mw5 ); 

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
