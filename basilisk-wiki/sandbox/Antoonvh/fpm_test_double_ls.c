/**
# Testing the double least squares method

$$A\mathbf{x} = \mathbf{b}$$

with

$$x_1 + x_2 = 1;$$
$$x_1 + x_2 + x_3 = 2$$

You can see the ["result" here](fpm_test_double_ls.c/out)

$\mathbf{b}$ os somewhat close to its initial vector and $\mathbf{x}$
satisfies the conditions.
*/

#include "fpm_double_ls.h"

int main() {
  int m = 4;
  double b[m];
  double bv = 1;
  for (int i = 0; i < m; i++)
    b[i] = bv;
  int n = 3;
  double x[n];
  double A[m*n];
  for (int i = 0; i < n*m; i++)
    A[i] = noise();
  int c = 2;
  double C[6] = {1, 0.5,  // col 1
		 1, 0.5,  // col 2
		 0, 0.5}; // col 3
  double a[c];
  for (int i = 0; i < c; i++)
    a[i] = 1;
  
  update_double_ls (A, m, n, b, C, c, a);
  // rhs
  for (int i = 0; i < m; i++) {
    b[i] += bv;
    printf ("b[%d] = %g\n",i, b[i]);
  }

  //Check the proposed rhs
  int lwork = -1;
  double work_size;
  int one = 1;
  int info;
  dgels_("N", &m, &n, &one, A, &m, b, &m, &work_size, &lwork, &info);
  lwork = (int)work_size;
  double * work = (double *)malloc(lwork * sizeof(double));
  dgels_("N", &m, &n, &one, A, &m, b, &m, work, &lwork, &info);
  putc ('\n', stdout);

  for (int i = 0; i < n; i++)
    printf ("x[%d] = %g\n",i, b[i]);
}
