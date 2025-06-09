/**
## Double least-squares problem

Consider the system of linear equations,

$$A\mathbf{x} = \mathbf{b},$$

with known (full rank) matrix $A$. We wish to find the vector
$\mathbf{b}$ for which $\mathbf{x}^*$ is the least-squares solution of
the system. Under the constraints that,

$$P\mathbf{x}^* =\mathbf{a}$$
$$\mathbf{b} \text{such that,}|\mathbf{b} - \mathbf{b}_i|\ \text{is minimized}$$

with $C$ and $\mathbf{a}$, a matrix and a vector specifying the linear
condition on $\mathbf{x}^*$, and $\mathbf{b}_i$ some reference
righthandside vector.

Using the left inverse:
$\mathbf{x}^* = (A^TA)^{-1}A^T\mathbf{b} = D\mathbf{b}$

and,

$$P\mathbf{x}^* = PD(mathbf{b}_i + \Delta\mathbf{b}) = E(mathbf{b}_i + \Delta\mathbf{b}) =\mathbf{a}_i + E\delta\mathbf{b} =\mathbf{a}$$
Such that the least squears ($|\Delta\mathbfb}|#) problem is,
$$E\delta\mathbf{b} =\mathbf{a} - \mathbf{a}_i.$$
*/
//mat mul
extern void dgemm_(char *TRANSA, char *TRANSB,
                   int *M, int *N, int *K,
                   double *ALPHA, double *A, int *LDA,
                   double *B, int *LDB,
                   double *BETA, double *C, int *LDC);

// LU
extern void dgetrf_(int *m, int *n, double *a, int *lda,
                    int *ipiv, int *info);
// Inverse
extern void dgetri_(int *n, double *a, int *lda,
                    int *ipiv, double *work, int *lwork, int *info);
//least squares
extern void dgels_(char *trans, int *m, int *n, int *nrhs,
                   double *a, int *lda, double *b, int *ldb,
                   double *work, int *lwork, int *info);

void update_double_ls (double * A, int m, int n, // Matrix and its size
		       double * b,                //rhs (b_i in) and delta b_i out.
		       double * P, int k, double * a) { //contraint matrix, its nr of contraints, and rhs
  double alpha = 1, beta = 0;;
  int one = 1;
  char trans = 'T';
  char ntrans = 'N';
  double x[n];
  int ipiv[n], info;
  int lwork = -1;
  double work_size;
  double C[n*n], D[n*m];
  double E[k*m];
  // Compute (C = A^TA)
  dgemm_(&trans, &ntrans, &n, &n, &m, &alpha, A, &m, A, &m, &beta, C, &n);
  // C^{-1}
  dgetrf_(&n, &n, C, &n, ipiv, &info);
  dgetri_(&n, C, &n, ipiv, &work_size, &lwork, &info);
  lwork = (int)work_size;
  double *work = malloc(lwork * sizeof(double));
  // Step 2: Compute inverse
  dgetri_(&n, C, &n, ipiv, work, &lwork, &info);
  //  D = C^{-1} A^T to find D
  //printf ("merkar 1\n");
  dgemm_(&ntrans, &trans, &n, &m, &n, &alpha, C, &n, A, &m, &beta, D, &n);
  //printf ("merkar 2 %c %c\n", ntrans, ntrans);
  // Compute initial `x`, corresponding to the initial rhs vector
  dgemm_(&ntrans, &ntrans, &n, &one, &m, &alpha, D, &n, b, &m, &beta, x, &n);
  // compute a_i
  double ai[k];
  dgemm_(&ntrans, &ntrans, &k, &one, &n, &alpha, P, &k, x, &n, &beta, ai, &k);
  for (int i = 0; i < k; i++)
    a[i] -= ai[i];
  
  // The leastsquares-updates problem:
  // Compute E
  dgemm_(&ntrans, &ntrans, &k, &m, &n, &alpha, P, &k, D, &n, &beta, E, &k);
  // Least squares
  lwork = -1;
  for (int i = 0; i < k; i++)
    b[i] = a[i];
  dgels_(&ntrans, &k, &m, &one, E, &k, b, &m, &work_size, &lwork, &info);
  lwork = (int)work_size;
  work = (double *)realloc(work, lwork * sizeof(double));
  dgels_(&ntrans, &k, &m, &one, E, &k, b, &m, work, &lwork, &info);
  free (work);
}

/**
## Test

*[A 4-eq 3 unknown with 2 constraints tests](fpm_test_double_ls.c)

## Usage

* [Inverting the FPM Laplacian](fpm_poisson.h)

*/
