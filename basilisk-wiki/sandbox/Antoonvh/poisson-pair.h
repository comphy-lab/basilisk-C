/**
# A pairwise Poisson-problem solver

An adaptation of [poisson.h](/src/poisson.h) with a special treatment
of the embedded boundary where fluxes are paritioned. Please mind the
different weighing of the implied weighing of `alpha` with the
embedded boundary face fractions.
 */
#include "embed.h"
#include "poisson.h"

double kappav1 = 1, kappav2 = 1; // Diffusivity near the interface

#ifndef FLUXFUN
#define FLUXFUN (0) // net flux
#endif
typedef double * doublep;

void embed_flux_pair (Point point, scalar * s,
		      double * val, bool homogenious);

struct Poisson_pair {
  scalar *a, *b;
  (const) face vector alpha1; 
  (const) face vector alpha2;
  (const) scalar lambda1;
  (const) scalar lambda2;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  void (* embed_flux_pair) (Point, scalar *, double *, bool);
#endif
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void relax_pair (scalar * al, scalar * bl, int l, void * data)
{
  scalar a1 = al[0], b1 = bl[0];
  scalar a2 = al[1], b2 = bl[1];
  struct Poisson_pair * p = (struct Poisson_pair *) data;
  (const) face vector alpha1 = p->alpha1;
  (const) scalar lambda1 = p->lambda1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda2 = p->lambda2;
  
  /**
  We use either Jacobi (under)relaxation or we directly reuse values
  as soon as they are updated. For Jacobi, we need to allocate space

  for the new field *c*. Jacobi is useful mostly as it gives results
  which are independent of the order in which the cells are
  traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or
  a multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  
#if JACOBI
  scalar c1[], c2[];
#else
  scalar c1 = a1;
  scalar c2 = a2;
#endif
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. An
  important change is that the embedded boundary is weighted here*/
  
  foreach_level_or_leaf (l) {
    double n1 = - sq(Delta)*b1[], d1 = - lambda1[]*sq(Delta);
    double n2 = - sq(Delta)*b2[], d2 = - lambda2[]*sq(Delta);
    foreach_dimension() {
      n1 += alpha1.x[1]*a1[1]*fs.x[1] +
	alpha1.x[]*a1[-1]*fs.x[];
      d1 += alpha1.x[1]*fs.x[1] + alpha1.x[]*fs.x[];
      n2 += alpha2.x[1]*a2[1]*(1. - fs.x[1]) +
	alpha2.x[]*a2[-1]*(1. - fs.x[]);
      d2 += alpha2.x[1]*(1. - fs.x[1]) + alpha2.x[]*(1. - fs.x[]);
    }
#if EMBED
    if (p->embed_flux_pair) {
      doublep val = malloc(2*sizeof(double));
      p->embed_flux_pair (point, al, val, true);
      n1 -= val[0]*sq(Delta);
      n2 -= val[1]*sq(Delta);
      free (val);
    }
    if (!d1) {
      c1[] = b1[] = 0.;
    }
    else
#endif // EMBED
      c1[] = n1/d1;
#if EMBED
    if (!d2) {
      c2[] = b2[] = 0.;
    }
    else
#endif //EMBED
      c2[] = n2/d2;
  }
  /**
     For weighted Jacobi we under-relax with a weight of 2/3. */
#if JACOBI
  foreach_level_or_leaf (l) {
    a1[] = (a1[] + 2.*c1[])/3.;
    a2[] = (a2[] + 2.*c2[])/3.;
  }
#endif
  
#if TRASH
  scalar a11[], a12[];
  foreach_level_or_leaf (l) {
    a11[] = a1[];
    a12[] = a2[];
  }
  trash ({a1, a2});
  foreach_level_or_leaf (l) {
    a1[] = a11[];
    a2[] = a12[];
  }
#endif
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual_pair (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a1 = al[0], b1 = bl[0], res1 = resl[0];
  scalar a2 = al[1], b2 = bl[1], res2 = resl[1];
  
  struct Poisson_pair * p = (struct Poisson_pair *) data;
  (const) face vector alpha1 = p->alpha1;
  (const) scalar lambda1 = p->lambda1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda2 = p->lambda2;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g1[], g2[];
  foreach_face() {
    g1.x[] = alpha1.x[]*face_gradient_x (a1, 0)*fs.x[];
    g2.x[] = alpha2.x[]*face_gradient_x (a2, 0)*(1 - fs.x[]);
  }
  foreach (reduction(max:maxres), nowarning) {
    res1[] = b1[] - lambda1[]*a1[];
    res2[] = b2[] - lambda2[]*a2[];
    
    foreach_dimension() {
      res1[] -= (g1.x[1] - g1.x[])/Delta;
      res2[] -= (g2.x[1] - g2.x[])/Delta;
    }
#if EMBED
    if (p->embed_flux_pair) {
      doublep val = malloc(2*sizeof(double));
      val[0] = nodata;
      p->embed_flux_pair (point, al, val, false);
      res1[] += val[0]; //- e1*a1[];
      res2[] += val[1]; // - e2*a2[];
      free (val);
    }
#endif // EMBED    
    if (fabs (res1[]) > maxres)
      maxres = fabs (res1[]);
    if (fabs (res2[]) > maxres)
      maxres = fabs (res2[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    res1[] = b1[] - lambda1[]*a1[];
    res2[] = b2[] - lambda2[]*a2[];
    foreach_dimension() {
      res1[] += (alpha1.x[0]*face_gradient_x (a1, 0)*fs.x[] -
		 alpha1.x[1]*face_gradient_x (a1, 1)*fs.x[1])/Delta;
      res2[] += (alpha2.x[0]*face_gradient_x (a2, 0)*(1 - fs.x[]) -
		 alpha2.x[1]*face_gradient_x (a2, 1)*(1 - fs.x[1]))/Delta;
    }
#if EMBED
    if (p->embed_flux_pair) {
      doublep val = malloc(2*sizeof(double));
      p->embed_flux_pair (point, al, val, false);
      res1[] += val[0]; //- e1*a1[];
      res2[] += val[1]; // - e2*a2[];
      free (val);
    }
#endif // EMBED
    if (fabs (res1[]) > maxres)
      maxres = fabs (res1[]);
    if (fabs (res2[]) > maxres)
      maxres = fabs (res2[]);
  }
#endif // !TREE
  return maxres;
}

/**
## User interface
$$ */

mgstats poisson_pair (struct Poisson_pair p)
{

  /**
  If $\alpha$s or $\lambda$s are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha1.x.i)
    p.alpha1 = unityf;
  if (!p.alpha2.x.i)
    p.alpha2 = unityf;
  if (!p.lambda1.i)
    p.lambda1 = zeroc;
  if (!p.lambda2.i)
    p.lambda2 = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha1 = p.alpha1;
  scalar lambda1 = p.lambda1;
  face vector alpha2 = p.alpha2;
  scalar lambda2 = p.lambda2;
  restriction ({alpha1, lambda1, alpha2, lambda2});
  
  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  
  //scalar a = p.a, b = p.b;
#if EMBED
  if (!p.embed_flux_pair)
    p.embed_flux_pair = embed_flux_pair;
#endif // EMBED
  mgstats s = mg_solve (p.a, p.b, residual_pair, relax_pair,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;
  
  return s;
}

/**

   ## Flux partitioning

*/


foreach_dimension()
  void nearby_vali_x(Point point, scalar s, scalar cs, face vector fs,
		  coord n, coord p, double * v, double * d) {
  bool defined = true;
  foreach_dimension()
    if (defined && !(1 - fs.x[(n.x > 0.)]))
      defined = false;
  if (defined) {
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (1 - fs.x[i + (i < 0),j] && 1 - fs.y[i,j] && 1 - fs.y[i,j+1] &&
	  1 - cs[i,j-1] && 1 - cs[i,j] && 1 - cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
    double z = p.z + d*n.z;
    int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
    z -= k;
    bool defined = 1 - fs.x[i + (i < 0),j,k];
    for (int m = -1; m <= 1 && defined; m++)
      if (!(1 - fs.y[i,j,k+m]) || !(1 - fs.y[i,j+1,k+m]) ||
	  (1 - !fs.z[i,j+m,k]) || !(1 - fs.z[i,j+m,k+1]) ||
	  !( 1 - cs[i,j+m,k-1]) || !(1 - cs[i,j+m,k]) || !(1 - cs[i,j+m,k+1]))
	defined = false;
    if (defined)
      // bi-quadratic interpolation
      v[l] =
	quadratic (z,
		   quadratic (y1,
			      (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		   quadratic (y1,
			      (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		   quadratic (y1,
			      (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
    else
      break;
    }
  }
  if (v[0] == nodata) {
    d[0] = max(1e-3, fabs(p.x/n.x));
    v[0] = s[];
    v[1] = nodata;
  }
}

foreach_dimension()
  void nearby_val_x(Point point, scalar s, scalar cs, face vector fs,
		  coord n, coord p, double * v, double * d) {
  bool defined = true;
   
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined) {
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1]) {
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
      }
      
#else // dimension == 3
    double z = p.z + d[l]*n.z;
    int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
    z -= k;
    bool defined = fs.x[i + (i < 0),j,k];
    for (int m = -1; m <= 1 && defined; m++)
      if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	  !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	  !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	defined = false;
    if (defined)
      // bi-quadratic interpolation
      v[l] =
	quadratic (z,
		   quadratic (y1,
			      (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		   quadratic (y1,
			      (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		   quadratic (y1,
			      (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
    }
  }
  if (v[0] == nodata) {
    printf ("yo\n");
    d[0] = max(1e-3, fabs(p.x/n.x));
    v[0] = s[];
  }
}

foreach_dimension()
  void pair_gradient_x  (Point point, scalar * sl, scalar cs, face vector fs,
			 coord n, coord p, double g[2], bool homogenious) {
  foreach_dimension()
    n.x = - n.x;
  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }
  nearby_val_x(point, sl[0], cs, fs, n, p, v1, d1);
  foreach_dimension()
    n.x = - n.x;
  nearby_vali_x(point, sl[1], cs, fs, n, p, v2, d2);
  double FLUXv = homogenious ? 0 :FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    double Ts0 = (v1[0]*(kappav1)/(Delta*d1[0]) + v2[0]*(kappav2)/(Delta*d2[0]) + FLUXv)/
      ((kappav1)/(Delta*d1[0]) + (kappav2)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {
    double wbc1 = (-d1[0] - d1[1])/(d1[0]*d1[1]);
    double wbc2 = (-d2[0] - d2[1])/(d2[0]*d2[1]);
    double wv01 = (-d1[1])/(sq(d1[0]) - d1[0]*d1[1]);
    double wv02 = (-d2[1])/(sq(d2[0]) - d2[0]*d2[1]);
    double wv11 = (d1[0])/(d1[0]*d1[1] - sq(d1[1]));
    double wv12 = (d2[0])/(d2[0]*d2[1] - sq(d2[1]));
    double Ts0 = -(wv01*v1[0]*kappav1 + wv02*v2[0]*kappav2 + (FLUXv*Delta) +
		  kappav1*wv11*v1[1] + kappav2*wv12*v2[1])/(kappav1*wbc1 + kappav2*wbc2);
    
    g[0]  = -(wv11*v1[1] + wv01*v1[0] + wbc1*Ts0)/(Delta);
    g[1]  = -(wv12*v2[1] + wv02*v2[0] + wbc2*Ts0)/(Delta);
  }
}

double total_flux_embed (Point point, coord p, coord n);
#if EMBED
void embed_flux_pair (Point point, scalar * s, double * val, bool homogenious)
{
  val[0] = 0.;
  val[1] = 0.;
  if (cs[] >= 1. || cs[] <= 0.)
    return;
  double g12[2];
  g12[0] = nodata; g12[1] = nodata;
  
  coord n = facet_normal (point, cs, fs), p;
  double alpha = plane_alpha (cs[], n);
  double area = plane_area_center (n, alpha, &p);
  normalize (&n);
    
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y)) {
      pair_gradient_x (point, s, cs, fs, n, p, g12, homogenious);
    }

#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z)) {
      pair_graident_x (point, s, cs, fs, n, p, g12, homogenious);
    }
  }
  else if (fabs(n.y) >= fabs(n.z)) {
    pair_gradient_y (point, s, cs, fs, n, p, g12, homogenious);
  } else {
    pair_gradient_z (point, s, cs, fs, n, p, g12, homogenious);
  }
#endif // dimension == 3
  double mua1 = kappav1, mua2 = kappav2;
  assert (g12[0] != nodata);
  assert (g12[1] != nodata);
  val[0] = -mua1*g12[0]*area/Delta;
  val[1] = -mua2*g12[1]*area/Delta;
}
#endif

/**
## Tests
* [A simple Poisson problem split on to two domains](test-pair.c)

## Usage
* [An implicit conjugate heatflux diffusion solver](diffusion-pair.h)
 */
