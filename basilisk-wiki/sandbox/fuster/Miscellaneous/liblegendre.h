/**
# Library of functions based on Legendre polynomials
*/

// we need two layers of ghost cells for the 5-points stencil near boundaries
#define BGHOSTS 2

#include "grid/multigrid.h"
#include "utils.h"

typedef struct { 
  double *** data;
  int r[3], np[3];
} LData;


/* allocation/deallocation operations */
void allocate_LData(LData * M) {
  
  M->data =(double ***) malloc(M->np[0] * sizeof(double **));

  for(int i=0; i < M->np[0]; i++) {
    M->data[i]=(double **)malloc(M->np[1] * sizeof(double *));
    for(int j=0; j < M->np[1]; j++)
      M->data[i][j]=(double *)malloc(M->np[2] * sizeof(double));
  }

}

void deallocate_LData (LData * M) {

  for(int i=0;i < M->np[0];i++)
  {
    for(int j=0;j < M->np[1];j++)
      free(M->data[i][j]);
    free(M->data[i]);
  }

  free(M->data);

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
void Init_Matrix_Integration( LData * M, double x1[3], double x2[3]) {


  double *coefs[3];
  for (int d=0; d < 3; d++){

    double * x0 = (double*) malloc(M->np[d] * sizeof(double));
    for (int i=0; i< M->np[d]; i++) 
      x0[i]=i-M->r[d];

    coefs[d] = integrate_legendre1D (x0, M->np[d], x1[d], x2[d]);
    
    free(x0); 

  }

  for (int i=0; i< M->np[0]; i++) 
      for (int j=0; j< M->np[1]; j++) 
          for (int k=0; k< M->np[2]; k++) 
            M->data[i][j][k] = coefs[0][i]*coefs[1][j]*coefs[2][k];

  for (int d=0; d < 3; d++)
      free(coefs[d]);

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
void Init_Matrix_Derivate( LData * M, double x[]) {

  double * x0 = (double*) malloc(M->np[0] * sizeof(double));

  for (int i=0; i< M->np[0]; i++) 
      x0[i]=i-M->r[0];

  double * coefs[3];
  coefs[0] = derive_legendre1D (x0, M->np[0], x[0]);
  coefs[1] = eval_legendre1D (x0, M->np[1], x[1]); 
  coefs[2] = eval_legendre1D (x0, M->np[2], x[2]); 

  free(x0); 

  for (int i=0; i< M->np[0]; i++) 
      for (int j=0; j< M->np[1]; j++) 
          for (int k=0; k< M->np[2]; k++) 
            M->data[i][j][k] = coefs[0][i]*coefs[1][j]*coefs[2][k];

  for (int d=0; d < 3; d++)
      free(coefs[d]);

}

LData Matrix_Derivate( double x[], int np, int r ) {

  LData M;
  for (int d=0; d < dimension; d++){
      M.np[d] = np; 
      M.r[d] = r;
  }
  for (int d=dimension; d < 3; d++){
      M.np[d] = 1; 
      M.r[d] = 0;
  }

  allocate_LData(&M);

  Init_Matrix_Derivate( &M, x); 

  return M;
}

foreach_dimension()
static double eval_x (Point point, scalar f, LData * M) {

    double sum = 0.;
    for (int i=0; i< M->np[0]; i++) 
        for (int j=0; j< M->np[1]; j++) 
            for (int k=0; k< M->np[2]; k++) 
              sum += M->data[i][j][k]*f[-M->r[0] + i,-M->r[1] + j, -M->r[2] + k];

    return sum;
}


// This function gives the cell centered averaged value
// given a field variable using a stencil of length order/2
void Init_avg_variable ( scalar f, scalar favg, int order ) 
{
  
  LData LM;
  for (int d=0; d < dimension; d++){
      LM.np[d] = order; 
      LM.r[d] = order/2;
  }
  for (int d=dimension; d < 3; d++){
      LM.np[d] = 1; 
      LM.r[d] = 0;
  }

  allocate_LData(&LM);

  double x1[3] = {-1./2, -1./2, -1./2};
  double x2[3] = { 1./2,  1./2,  1./2};
  Init_Matrix_Integration( &LM, x1, x2);

  foreach() 
    favg[] = eval_x(point, f, &LM);

  boundary({favg});

  deallocate_LData(&LM);
}

void Variable_Gradient (scalar f, vector df, LData *M)
{
    foreach() 
        foreach_dimension () 
            df.x[] = eval_x (point, f, M)/Delta;
    
    boundary ((scalar *){df});
}