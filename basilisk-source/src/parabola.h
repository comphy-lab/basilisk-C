#include "utils.h"

#define PARABOLA_FIT_CENTER_WEIGHT .1

// Define this to use a x^iy^j polynomial with i = 0...NP-1, j = 0...NP-1
// #define NP 3

typedef struct {
  coord o;
#if dimension == 2 /* y = a[0]*x^2 + a[1]*x + a[2] */
  coord m;
  double ** M, rhs[3], a[3];
#else /* 3D z = a[0]*x^2 + a[1]*y^2 + a[2]*x*y + a[3]*x + a[4]*y + a[5] */
  double t[3][3];
# ifdef NP
  double ** M, rhs[NP*NP], a[NP*NP];
# else
  double ** M, rhs[6], a[6];
# endif
#endif /* 3D */
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  foreach_dimension()
    p->o.x = o.x;
#if dimension == 2
  foreach_dimension()
    p->m.x = m.x;
  normalize (&p->m);
  int n = 3;
#else /* 3D */
  double max;
  coord nx = {0., 0., 0.}, ny, nz;
  int d = 0;

  foreach_dimension()
    nz.x = m.x;
  normalize (&nz);
  max = sq(nz.x);
  /* build a vector orthogonal to nz */
  if (sq(nz.y) > max) { max = sq(nz.y); d = 1; }
  if (sq(nz.z) > max) d = 2;
  switch (d) {
  case 0: nx.x = - nz.z/nz.x; nx.z = 1.0; break;
  case 1: nx.y = - nz.z/nz.y; nx.z = 1.0; break;
  case 2: nx.z = - nz.x/nz.z; nx.x = 1.0; break;
  }
  normalize (&nx);

  /* build a second vector orthogonal to nx and nz */
  foreach_dimension()
    ny.x = nz.y*nx.z - nz.z*nx.y;

  /* transformation matrix from (i,j,k) to (nx, ny, nz) */
  p->t[0][0] = nx.x; p->t[0][1] = nx.y; p->t[0][2] = nx.z;
  p->t[1][0] = ny.x; p->t[1][1] = ny.y; p->t[1][2] = ny.z;
  p->t[2][0] = nz.x; p->t[2][1] = nz.y; p->t[2][2] = nz.z;
# ifdef NP
  int n = NP*NP;
# else
  int n = 6;
# endif
#endif /* 3D */
  p->M = (double **) matrix_new (n, n, sizeof(double));  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{
#if dimension == 2
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#else /* 3D */
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y, z1 = m.z - p->o.z;
  double x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  double y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  double z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
# ifdef NP
  for (int i = 0; i < NP; i++)
    for (int j = 0; j < NP; j++) {
      for (int k = 0; k < NP; k++)
	for (int l = 0; l < NP; l++)
	  p->M[i*NP + j][k*NP + l] += w*pow(x, i + k)*pow(y, j + l);
      p->rhs[i*NP + j] += w*z*pow(x, i)*pow(y, j);
    }
# else // !NP 
  double x2 = x*x, x3 = x2*x, x4 = x3*x;
  double y2 = y*y, y3 = y2*y, y4 = y3*y;
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2; 
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;
# endif // !NP
#endif /* 3D */
}

static double parabola_fit_solve (ParabolaFit * p)
{
#if dimension == 2
  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  double pivmin = matrix_inverse (p->M, 3, 1e-10);
  if (pivmin) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else /* this may be a degenerate/isolated interface fragment */
    p->a[0] = p->a[1] = 0.;
#else /* 3D */
# ifdef NP
  double pivmin = matrix_inverse (p->M, NP*NP, 1e-10);
  if (pivmin)
    for (int i = 0; i < NP*NP; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < NP*NP; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else /* this may be a degenerate/isolated interface fragment */
    for (int i = 0; i < NP*NP; i++)
      p->a[i] = 0.;
# else // !NP
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  for (int i = 1; i < 6; i++)
    for (int j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  double pivmin = matrix_inverse (p->M, 6, 1e-10);
  if (pivmin)
    for (int i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < 6; j++)
	p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else /* this may be a degenerate/isolated interface fragment */
    for (int i = 0; i < 6; i++)
      p->a[i] = 0.;
# endif // !NP
#endif /* 3D */  
  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
				      double kappamax, double * kmax)
{
  double kappa;
#if dimension == 2
  double dnm = 1.[0] + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#else /* 3D */
# ifdef NP
  double hxx = 2.*p->a[2*NP], hyy = 2.*p->a[2], hxy = p->a[NP + 1];
  double hx = p->a[NP], hy = p->a[1];
# else
  double hxx = 2.*p->a[0], hyy = 2.*p->a[1], hxy = p->a[2];
  double hx = p->a[3], hy = p->a[4];
# endif
  double dnm = 1. + sq(hx) + sq(hy);
  kappa = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
    /sqrt (dnm*dnm*dnm);
  if (kmax) {
    double kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    double a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }
#endif /* 3D */
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}

#if AXI
static void parabola_fit_axi_curvature (const ParabolaFit * p,
					double r, double h,
					double * kappa, double * kmax)
{
  double nr = (p->m.x*p->a[1] + p->m.y)/sqrt (1. + sq(p->a[1]));
  /* limit the minimum radius to half the grid size */
  double kaxi = nr/max(r, h/2.);
  *kappa += kaxi;
  if (kmax)
    *kmax = max (*kmax, fabs (kaxi));
}
#endif /* 2D */
