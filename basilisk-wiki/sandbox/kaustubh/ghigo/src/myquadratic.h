#include "utils.h"

typedef struct {
  coord o; //
  int n, m;
#if dimension == 2 // m = a[0] + a[1]*x + a[2]*y + a[3]*x*y + a[4]*x^2 + a[5]*y^2
  double ** M, rhs[6], a[6];
#else // dimension == 3, m = a[0] + a[1]*x + a[2]*y + a[3]*z
      //                     a[4]*xy + a[5]*xz + a[6]*yz +
      //                     a[7]*x^2 + a[8]*y^2 + a[9]*z^2
  double ** M, rhs[10], a[10];
#endif // dimension
  bool linear; // Force linear extrapolation instead of quadratic
} QuadraticFit;

#if dimension == 2
int neigh = 3; // Number of neighbors to use above matrix dimension to improve matrix conditioning
#else // dimension == 3
int neigh = 5;
#endif // dimension
double pivtol = 1.e-10; // Tolerance for matrix inversion

static void quadratic_fit_init (QuadraticFit * p, coord o, bool linear)
{
  // Origin
  foreach_dimension()
    p->o.x = o.x;
#if dimension == 2
  p->n = 6;
  p->m = 6; // Dimension of reduced matrix if necessary
#else // dimension == 3
  p->n = 10;
  p->m = 10; // Dimension of reduced matrix if necessary
#endif // dimension
  // Matrix
  p->M = (double **) matrix_new (p->n, p->n, sizeof(double));  
  for (int i = 0; i < p->n; i++) {
    for (int j = 0; j < p->n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
  // Linear or quadratic extrapolation
  p->linear = linear;
}

static void quadratic_fit_add (QuadraticFit * p, coord o, double m)
{
  /**
  We update here the coefficients in the lower triangle matrix. */
  
#if dimension == 2
  // Relative coordinates to the origin
  double x1 = o.x - p->o.x, y1 = o.y - p->o.y;
  // Useful variables
  double x2 = x1*x1, y2 = y1*y1, xy = x1*y1;
  double x3 = x1*x2, y3 = y1*y2, x2y = x1*xy, xy2 = y1*xy;
  double x4 = x1*x3, y4 = y1*y3, x3y = x1*x2y, xy3 = y1*xy2, x2y2 = xy*xy;
  // Fill the matrix X^TX
  p->M[0][0] += 1.;
  p->M[1][0] += x1; p->M[1][1] += x2;
  p->M[2][0] += y1; p->M[2][1] += xy;   p->M[2][2] += y2;
  p->M[3][0] += xy; p->M[3][1] += x2y;  p->M[3][2] += xy2;  p->M[3][3] += x2y2;
  p->M[4][0] += x2; p->M[4][1] += x3;   p->M[4][2] += x2y;  p->M[4][3] += x3y;  p->M[4][4] += x4;
  p->M[5][0] += y2; p->M[5][1] += xy2;  p->M[5][2] += y3;   p->M[5][3] += xy3;  p->M[5][4] += x2y2; p->M[5][5] += y4;
  p->rhs[0]  += m;  p->rhs[1]  += x1*m; p->rhs[2]  += y1*m; p->rhs[3]  += xy*m; p->rhs[4]  += x2*m; p->rhs[5] += y2*m;
#else // dimension == 3
  // Relative coordinates to the origin
  double x1 = o.x - p->o.x, y1 = o.y - p->o.y, z1 = o.z - p->o.z;
  // Useful variables
  double x2 = x1*x1, y2 = y1*y1, z2 = z1*z1, xy = x1*y1, xz = x1*z1, yz = y1*z1;
  double x3 = x1*x2, y3 = y1*y2, z3 = z1*z2, x2y = x1*xy, x2z = x1*xz, xy2 = y1*xy, y2z = y1*yz, xz2 = z1*xz, yz2 = z1*yz, xyz = x1*y1*z1;
  double x4 = x1*x3, y4 = y1*y3, z4 = z1*z3, x3y = x1*x2y, x3z = x1*x2z, xy3 = y1*xy2, y3z = y1*y2z, xz3 = z1*xz2, yz3 = z1*yz2, x2y2 = xy*xy, x2z2 = xz*xz, y2z2 = yz*yz, xyz2 = xy*z2, xy2z = xz*y2, x2yz = x2*yz;
  // Fill the matrix X^TX
  p->M[0][0] += 1.;
  p->M[1][0] += x1; p->M[1][1] += x2;
  p->M[2][0] += y1; p->M[2][1] += xy;  p->M[2][2] += y2;
  p->M[3][0] += z1; p->M[3][1] += xz;  p->M[3][2] += yz;  p->M[3][3] += z2;
  p->M[4][0] += xy; p->M[4][1] += x2y; p->M[4][2] += xy2; p->M[4][3] += xyz; p->M[4][4] += x2y2;
  p->M[5][0] += xz; p->M[5][1] += x2z; p->M[5][2] += xyz; p->M[5][3] += xz2; p->M[5][4] += x2yz; p->M[5][5] += x2z2;
  p->M[6][0] += yz; p->M[6][1] += xyz; p->M[6][2] += y2z; p->M[6][3] += yz2; p->M[6][4] += xy2z; p->M[6][5] += xyz2; p->M[6][6] += y2z2;
  p->M[7][0] += x2; p->M[7][1] += x3;  p->M[7][2] += x2y; p->M[7][3] += x2z; p->M[7][4] += x3y;  p->M[7][5] += x3z;  p->M[7][6] += x2yz; p->M[7][7] += x4;
  p->M[8][0] += y2; p->M[8][1] += xy2; p->M[8][2] += y3;  p->M[8][3] += y2z; p->M[8][4] += xy3;  p->M[8][5] += xy2z; p->M[8][6] += y3z;  p->M[8][7] += x2y2; p->M[8][8] += y4;
  p->M[9][0] += z2; p->M[9][1] += xz2; p->M[9][2] += yz2; p->M[9][3] += z3;  p->M[9][4] += xyz2; p->M[9][5] += xz3;  p->M[9][6] += yz3;  p->M[9][7] += x2z2; p->M[9][8] += y2z2; p->M[9][9] += z4;
  p->rhs[0]  += m;  p->rhs[1] += x1*m; p->rhs[2] += y1*m; p->rhs[3] += z1*m; p->rhs[4] += xy*m;  p->rhs[5] += xz*m;  p->rhs[6] += yz*m;  p->rhs[7] += x2*m;  p->rhs[8] += y2*m;  p->rhs[9] += z2*m;
#endif // dimension  
}

static void quadratic_fit_set (QuadraticFit * p, const int nc)
{
  /**
  We set the symmetric coefficients (upper triangle matrix) */

  for (int i = 0; i < p->n; i++)
    for (int j = i + 1; j < p->n; j++)
      p->M[i][j] = p->M[j][i];
  
  /**
  For degenereate cases, when not enough neighbors are available, we
  use a linear fit or simple injection and therefore set the quadratic
  (and linear for injection) matrix coefficient to 0. 

  Note that we switch to a linear fit if the number of neighbors
  available is smaller than *dimension + neigh*. This is to avoid
  inverting a matrix with a large condition number, which would
  generate a inaccurate solution. An alternative would be to use a
  $QR$ decomposition instead of matrix inversion. */

  // Injection
  if (nc < (dimension + 1))
    p->m = 1;
  // Linear interpolation/extrapolation
  else if (nc < p->n + neigh || p->linear)
    p->m = dimension + 1;

  if (p->m != p->n)
    for (int i = 0; i < p->n; i++)
      for (int j = p->m; j < p->n; j++)
	p->M[i][j] = p->M[j][i] = 0.;
}

static double quadratic_fit_solve (QuadraticFit * p, const int nc)
{
  quadratic_fit_set (p, nc);
  double pivmin = matrix_inverse (p->M, p->m, pivtol);
  if (pivmin)
    for (int i = 0; i < p->n; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < p->n; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else /* this may be a degenerate/isolated interface fragment */  
    for (int i = 0; i < p->n; i++)
      p->a[i] = 0.;
  matrix_free (p->M);
  return pivmin;
}

