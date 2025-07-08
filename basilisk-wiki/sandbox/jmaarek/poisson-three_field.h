/**
# A pairwise Poisson-problem solver

An adaptation of [poisson.h](/src/poisson.h) with a special treatment
of the embedded boundary where fluxes are paritioned. Please mind the
different weighing of the implied weighing of `alpha` with the
embedded boundary face fractions.
 */
//#include "poisson.h"
#include "cross_sectional_area.h"

double kappav1 = 1, kappav2 = 1, kappav3 = 1;// Diffusivity near the interface

#ifndef FLUXFUN
#define FLUXFUN (0.0) // net flux
#endif
typedef double * doublep;

void embed_flux_three_field (Point point, scalar * s,
		      double * val, bool homogenious, scalar, scalar, scalar, face vector);

struct Poisson_three_field {
  scalar *a, *b;
  scalar cs2, cs3;
  face vector fs2;
  scalar f;
  (const) face vector alpha1; 
  (const) face vector alpha2;
  (const) face vector alpha3;
  (const) scalar lambda1;
  (const) scalar lambda2;
  (const) scalar lambda3;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  void (* embed_flux_three_field) (Point, scalar *, double *, bool, scalar, scalar, scalar, face vector);
#endif
};


/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void relax_three_field (scalar * al, scalar * bl, int l, void * data){
  scalar a1 = al[0], b1 = bl[0];
  scalar a2 = al[1], b2 = bl[1];
  scalar a3 = al[2], b3 = bl[2];
  struct Poisson_three_field * p = (struct Poisson_three_field *) data;

  scalar f = p->f;
  scalar cs2 = p->cs2;
  scalar cs3 = p->cs3;
  face vector fs2 = p->fs2;

  (const) face vector alpha1 = p->alpha1;
  (const) scalar lambda1 = p->lambda1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda2 = p->lambda2;
  (const) face vector alpha3 = p->alpha3;
  (const) scalar lambda3 = p->lambda3;
  
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
  scalar c1[], c2[], c3[];
#else
  scalar c1 = a1;
  scalar c2 = a2;
  scalar c3 = a3;
#endif
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. An
  important change is that the embedded boundary is weighted here*/
  
  foreach_level_or_leaf (l) {
    double n1 = - sq(Delta)*b1[], d1 = - lambda1[]*sq(Delta);
    double n2 = - sq(Delta)*b2[], d2 = - lambda2[]*sq(Delta);
    double n3 = - sq(Delta)*b3[], d3 = - lambda3[]*sq(Delta);
    foreach_dimension() {
      n1 += alpha1.x[1]*a1[1]*(1-fs.x[1]) + alpha1.x[]*a1[-1]*(1-fs.x[]);
      d1 += alpha1.x[1]*(1-fs.x[1]) + alpha1.x[]*(1-fs.x[]);
      n2 += alpha2.x[1]*a2[1]*(fs2.x[1]) + alpha2.x[]*a2[-1]*(fs2.x[]);
      d2 += alpha2.x[1]*(fs2.x[1]) + alpha2.x[]*(fs2.x[]);
      n3 += alpha3.x[1]*a3[1]*(fs.x[1]-fs2.x[1]) + alpha3.x[]*a3[-1]*(fs.x[] - fs2.x[]);
      d3 += alpha3.x[1]*(fs.x[1] - fs2.x[1]) + alpha3.x[]*(fs.x[] - fs2.x[]);
    }
#if EMBED
    if (p->embed_flux_three_field) {
      doublep val = malloc(3*sizeof(double));
      p->embed_flux_three_field (point, al, val, true, f, cs2, cs3, fs2);
      // it is observed to be required to modify the relaxation function by weighting the correction by the corresponding volume fraction to avoid problems in triple point cells
      n1 -= val[0]*sq(Delta)*(1-cs[]);
      n2 -= val[1]*sq(Delta)*(cs2[]);
      n3 -= val[2]*sq(Delta)*(cs3[]);
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
#if EMBED
    if (!d3) {
      c3[] = b3[] = 0.;
    }
    else
#endif //EMBED
      c3[] = n3/d3;
  }
  /**
     For weighted Jacobi we under-relax with a weight of 2/3. */
#if JACOBI
  foreach_level_or_leaf (l) {
    a1[] = (a1[] + 2.*c1[])/3.;
    a2[] = (a2[] + 2.*c2[])/3.;
    a3[] = (a3[] + 2.*c3[])/3.;
  }
#endif
  
#if TRASH
  scalar a11[], a12[], a13[];
  foreach_level_or_leaf (l) {
    a11[] = a1[];
    a12[] = a2[];
    a13[] = a3[];
  }
  trash ({a1, a2, a3});
  foreach_level_or_leaf (l) {
    a1[] = a11[];
    a2[] = a12[];
    a3[] = a13[];
  }
#endif
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */


static double residual_three_field (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a1 = al[0], b1 = bl[0], res1 = resl[0];
  scalar a2 = al[1], b2 = bl[1], res2 = resl[1];
  scalar a3 = al[2], b3 = bl[2], res3 = resl[2];
  
  struct Poisson_three_field * p = (struct Poisson_three_field *) data;

  scalar f = p->f;
  scalar cs2 = p->cs2;
  scalar cs3 = p->cs3;
  face vector fs2 = p->fs2;

  (const) face vector alpha1 = p->alpha1;
  (const) scalar lambda1 = p->lambda1;
  (const) face vector alpha2 = p->alpha2;
  (const) scalar lambda2 = p->lambda2;
  (const) face vector alpha3 = p->alpha3;
  (const) scalar lambda3 = p->lambda3;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g1[], g2[], g3[];
  foreach_face() {
    g1.x[] = alpha1.x[]*face_gradient_x (a1, 0)*(1-fs.x[]);
    g2.x[] = alpha2.x[]*face_gradient_x (a2, 0)*(fs2.x[]);
    g3.x[] = alpha3.x[]*face_gradient_x (a3, 0)*(fs.x[]-fs2.x[]);
  }
  foreach (reduction(max:maxres), nowarning) {
    res1[] = b1[] - lambda1[]*a1[];
    res2[] = b2[] - lambda2[]*a2[];
    res3[] = b3[] - lambda3[]*a3[];
    
    foreach_dimension() {
      res1[] -= (g1.x[1] - g1.x[])/Delta;
      res2[] -= (g2.x[1] - g2.x[])/Delta;
      res3[] -= (g3.x[1] - g3.x[])/Delta;
    }
#if EMBED
    if (p->embed_flux_three_field) {
      doublep val = malloc(3*sizeof(double));
      val[0] = nodata;
      val[1] = nodata;
      val[2] = nodata;
      p->embed_flux_three_field (point, al, val, false, f, cs2, cs3, fs2);
      assert (val[0] != nodata);
      assert (val[1] != nodata);
      assert (val[2] != nodata);
      res1[] += val[0]; //- e1*a1[];
      res2[] += val[1]; // - e2*a2[];
      res3[] += val[2]; // - e2*a2[];
      free (val);
    }
#endif // EMBED    
    if (fabs (res1[]) > maxres)
      maxres = fabs (res1[]);
    if (fabs (res2[]) > maxres)
      maxres = fabs (res2[]);
    if (fabs (res3[]) > maxres)
      maxres = fabs (res3[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres), nowarning) {
    res1[] = b1[] - lambda1[]*a1[];
    res2[] = b2[] - lambda2[]*a2[];
    res3[] = b3[] - lambda3[]*a3[];
    foreach_dimension() {
      res1[] += (alpha1.x[0]*face_gradient_x (a1, 0)*(1-fs.x[]) - alpha1.x[1]*face_gradient_x (a1, 1)*(1-fs.x[1]))/Delta;
      res2[] += (alpha2.x[0]*face_gradient_x (a2, 0)*(fs2.x[]) - alpha2.x[1]*face_gradient_x (a2, 1)*(fs2.x[1]))/Delta;
      res3[] += (alpha3.x[0]*face_gradient_x (a3, 0)*(fs.x[]-fs2.x[]) - alpha3.x[1]*face_gradient_x (a3, 1)*(fs.x[1] - fs2.x[1]))/Delta;
    }
#if EMBED
    if (p->embed_flux_three_field) {
      doublep val = malloc(3*sizeof(double));
      p->embed_flux_three_field (point, al, val, false, f, cs2, cs3, fs2);
      res1[] += val[0]; //- e1*a1[];
      res2[] += val[1]; // - e2*a2[];
      res3[] += val[2]; // - e2*a2[];
      free (val);
    }
#endif // EMBED
    if (fabs (res1[]) > maxres)
      maxres = fabs (res1[]);
    if (fabs (res2[]) > maxres)
      maxres = fabs (res2[]);
    if (fabs (res3[]) > maxres)
      maxres = fabs (res3[]);
  }
#endif // !TREE
  return maxres;
}


/**
## User interface
$$ */



mgstats poisson_three_field (struct Poisson_three_field p)
{

  face vector alpha1 = p.alpha1;
  scalar lambda1 = p.lambda1;
  face vector alpha2 = p.alpha2;
  scalar lambda2 = p.lambda2;
  face vector alpha3 = p.alpha3;
  scalar lambda3 = p.lambda3;

  scalar cs2 = p.cs2;
  scalar cs3 = p.cs3;
  scalar f = p.f;
  face vector fs2 = p.fs2;


  /**
  If $\alpha$s or $\lambda$s are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!alpha1.x.i)
    alpha1 = unityf;
  if (!alpha2.x.i)
    alpha2 = unityf;
  if (!alpha3.x.i)
    alpha3 = unityf;
  if (!lambda1.i)
    lambda1 = zeroc;
  if (!lambda2.i)
    lambda2 = zeroc;
  if (!lambda3.i)
    lambda3 = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */


  restriction ({alpha1, lambda1, alpha2, lambda2, alpha3, lambda3});
  restriction ({cs, cs2, cs3, fs2, f});
  
  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  
  //scalar a = p.a, b = p.b;
#if EMBED
  if (!p.embed_flux_three_field)
    p.embed_flux_three_field = embed_flux_three_field;
#endif // EMBED
  mgstats s = mg_solve (p.a, p.b, residual_three_field, relax_three_field,
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
  void nearby_val1_x(Point point, scalar s, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2, coord n, double alpha, coord p, coord p2, double * v, double * d) {
    foreach_dimension()
      n.x = - n.x;
    bool defined = true;
    //in triple point cells we use the degenerate case
    if ((cs[] > 0) && (cs[] < 1) && (cs2[] > 0) && (cs2[] < 1) && (cs3[] > 0) && (cs3[] < 1))
      defined = false;
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
        if ((1 - fs.x[i + (i < 0),j]) && (1 - fs.y[i,j]) && (1 - fs.y[i,j+1]) &&
      (1 - cs[i,j-1]) && (1 - cs[i,j]) && (1 - cs[i,j+1]))
    v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
  #else // dimension == 3
      double z = p.z + d*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = 1 - fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
        if (!(1 - fs.y[i,j,k+m]) || !(1 - fs.y[i,j+1,k+m]) ||
      !(1 - fs.z[i,j+m,k]) || !(1 - fs.z[i,j+m,k+1]) ||
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
      d[0] =  max(fabs(n.x*p2.x + n.y*p2.y + n.z*p2.z + alpha), 1e-3);
      v[0] = s[];
      v[1] = nodata;
    }
  }


foreach_dimension()
  void nearby_val2_x(Point point, scalar s, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2, coord n, double alpha, coord p, coord p2, double * v, double * d) {
    foreach_dimension()
      n.x = - n.x;
    bool defined = true;

    if ((cs[] > 0) && (cs[] < 1) && (cs2[] > 0) && (cs2[] < 1) && (cs3[] > 0) && (cs3[] < 1))
      defined = false;
     
    foreach_dimension()
      if (defined && !fs2.x[(n.x > 0.)])
        defined = false;
    if (defined) {
      for (int l = 0; l <= 1; l++) {
        int i = (l + 1)*sign(n.x);
        d[l] = (i - p.x)/n.x;
        double y1 = p.y + d[l]*n.y;
        int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
        y1 -= j;
  #if dimension == 2
        if (fs2.x[i + (i < 0),j] && fs2.y[i,j] && fs2.y[i,j+1] &&
      cs2[i,j-1] && cs2[i,j] && cs2[i,j+1]) {
    v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
        }
        
  #else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs2.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
        if (!fs2.y[i,j,k+m] || !fs2.y[i,j+1,k+m] ||
      !fs2.z[i,j+m,k] || !fs2.z[i,j+m,k+1] ||
      !cs2[i,j+m,k-1] || !cs2[i,j+m,k] || !cs2[i,j+m,k+1])
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
      d[0] =  max(fabs(n.x*p2.x + n.y*p2.y + n.z*p2.z + alpha), 1e-3);
      v[0] = s[];
      v[1] = nodata;
    }
  }


foreach_dimension()
  void nearby_val3_x(Point point, scalar s, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2, coord n, double alpha, coord p, coord p2, double * v, double * d) {
    foreach_dimension()
      n.x = - n.x;
    bool defined = true;

    if ((cs[] > 0) && (cs[] < 1) && (cs2[] > 0) && (cs2[] < 1) && (cs3[] > 0) && (cs3[] < 1))
      defined = false;
    foreach_dimension()
      if (defined && !(fs.x[(n.x > 0.)] - fs2.x[(n.x > 0.)]))
        defined = false;
    if (defined) {
      for (int l = 0; l <= 1; l++) {
        int i = (l + 1)*sign(n.x);
        d[l] = (i - p.x)/n.x;
        double y1 = p.y + d[l]*n.y;
        int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
        y1 -= j;
  #if dimension == 2
        if ((fs.x[i + (i < 0),j] - fs2.x[i + (i < 0),j]) && (fs.y[i,j] - fs2.y[i,j]) && (fs.y[i,j+1] - fs2.y[i,j+1]) &&
      cs3[i,j-1] && cs3[i,j] && cs3[i,j+1])
    v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
  #else // dimension == 3
      double z = p.z + d*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k] - fs2.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
        if (!(fs.y[i,j,k+m] - fs2.y[i,j,k+m]) || !(fs.y[i,j+1,k+m] - fs2.y[i,j+1,k+m]) ||
          !(fs.z[i,j+m,k] - fs2.z[i,j+m,k]) || !(fs.z[i,j+m,k+1] - fs2.z[i,j+m,k+1]) ||
          !( cs3[i,j+m,k-1]) || !(cs3[i,j+m,k]) || !(cs3[i,j+m,k+1]))
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
      d[0] = max(fabs(n.x*p2.x + n.y*p2.y + n.z*p2.z + alpha), 1e-3);
      v[0] = s[];
      v[1] = nodata;
    }
  }


/* This formulation uses quadratic interpolation and a four point stencil. This might lead to errors when combined with an advection scheme that is inaccurate in interfacial cells   */

/*
foreach_dimension()
  void pair_gradient_12_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }
  nearby_val1_x(point, sl[0], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val2_x(point, sl[1], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);

  //fprintf (fout, "temp %g %g %g %g %g %g %g %g\n", d1[0], d1[1], d2[0], d2[1], v1[0], v1[1], v2[0], v2[1]);

  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 1\n");
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

foreach_dimension()
  void pair_gradient_13_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }
  nearby_val1_x(point, sl[0], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val3_x(point, sl[2], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);
  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 2\n");
    double Ts0 = (v1[0]*(kappav1)/(Delta*d1[0]) + v2[0]*(kappav3)/(Delta*d2[0]) + FLUXv)/
      ((kappav1)/(Delta*d1[0]) + (kappav3)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {
    double wbc1 = (-d1[0] - d1[1])/(d1[0]*d1[1]);
    double wbc2 = (-d2[0] - d2[1])/(d2[0]*d2[1]);
    double wv01 = (-d1[1])/(sq(d1[0]) - d1[0]*d1[1]);
    double wv02 = (-d2[1])/(sq(d2[0]) - d2[0]*d2[1]);
    double wv11 = (d1[0])/(d1[0]*d1[1] - sq(d1[1]));
    double wv12 = (d2[0])/(d2[0]*d2[1] - sq(d2[1]));
    double Ts0 = -(wv01*v1[0]*kappav1 + wv02*v2[0]*kappav3 + (FLUXv*Delta) +
      kappav1*wv11*v1[1] + kappav3*wv12*v2[1])/(kappav1*wbc1 + kappav3*wbc2);
    
    g[0]  = -(wv11*v1[1] + wv01*v1[0] + wbc1*Ts0)/(Delta);
    g[1]  = -(wv12*v2[1] + wv02*v2[0] + wbc2*Ts0)/(Delta);
  }
}

foreach_dimension()
  void pair_gradient_23_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }

  nearby_val2_x(point, sl[1], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val3_x(point, sl[2], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);

  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 3\n");
    double Ts0 = (v1[0]*(kappav2)/(Delta*d1[0]) + v2[0]*(kappav3)/(Delta*d2[0]) + FLUXv)/
      ((kappav2)/(Delta*d1[0]) + (kappav3)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {
    double wbc1 = (-d1[0] - d1[1])/(d1[0]*d1[1]);
    double wbc2 = (-d2[0] - d2[1])/(d2[0]*d2[1]);
    double wv01 = (-d1[1])/(sq(d1[0]) - d1[0]*d1[1]);
    double wv02 = (-d2[1])/(sq(d2[0]) - d2[0]*d2[1]);
    double wv11 = (d1[0])/(d1[0]*d1[1] - sq(d1[1]));
    double wv12 = (d2[0])/(d2[0]*d2[1] - sq(d2[1]));
    double Ts0 = -(wv01*v1[0]*kappav2 + wv02*v2[0]*kappav3 + (FLUXv*Delta) +
      kappav2*wv11*v1[1] + kappav3*wv12*v2[1])/(kappav2*wbc1 + kappav3*wbc2);
    
    g[0]  = -(wv11*v1[1] + wv01*v1[0] + wbc1*Ts0)/(Delta);
    g[1]  = -(wv12*v2[1] + wv02*v2[0] + wbc2*Ts0)/(Delta);
  }
}*/


/* This formulation uses volume fraction weighted interpolation and a four point stencil. This should be more robust when combined with an advection scheme that is inaccurate in interfacial cells   */


foreach_dimension()
  void pair_gradient_12_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }
  nearby_val1_x(point, sl[0], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val2_x(point, sl[1], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);

  //fprintf (fout, "temp %g %g %g %g %g %g %g %g\n", d1[0], d1[1], d2[0], d2[1], v1[0], v1[1], v2[0], v2[1]);

  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 1\n");
    double Ts0 = (v1[0]*(kappav1)/(Delta*d1[0]) + v2[0]*(kappav2)/(Delta*d2[0]) + FLUXv)/
      ((kappav1)/(Delta*d1[0]) + (kappav2)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {

    double f_local = (1-cs[]); // local fraction value at the interface
    double denom1 = f_local/d1[0] + (1. - f_local)/d1[1];
    double denom2 = (1. - f_local)/d2[0] + f_local/d2[1];

    double num1 = f_local * v1[0]/d1[0] + (1. - f_local) * v1[1]/d1[1];
    double num2 = (1. - f_local) * v2[0]/d2[0] + f_local * v2[1]/d2[1];

    double Ts0 = (kappav1 * num1 + kappav2 * num2 + FLUXv) / (kappav1 * denom1 + kappav2 * denom2);

    
    g[0] = (f_local*(Ts0 - v1[0])/d1[0] + (1. - f_local)*(Ts0 - v1[1])/d1[1]) / Delta;
    g[1] = ((1. - f_local)*(Ts0 - v2[0])/d2[0] + f_local*(Ts0 - v2[1])/d2[1]) / Delta;
  }
}

foreach_dimension()
  void pair_gradient_13_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }
  nearby_val1_x(point, sl[0], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val3_x(point, sl[2], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);
  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 2\n");
    double Ts0 = (v1[0]*(kappav1)/(Delta*d1[0]) + v2[0]*(kappav3)/(Delta*d2[0]) + FLUXv)/
      ((kappav1)/(Delta*d1[0]) + (kappav3)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {

    double f_local = (1-cs[]); // local fraction value at the interface
    double denom1 = f_local/d1[0] + (1. - f_local)/d1[1];
    double denom2 = (1. - f_local)/d2[0] + f_local/d2[1];

    double num1 = f_local * v1[0]/d1[0] + (1. - f_local) * v1[1]/d1[1];
    double num2 = (1. - f_local) * v2[0]/d2[0] + f_local * v2[1]/d2[1];

    double Ts0 = (kappav1 * num1 + kappav3 * num2 + FLUXv) / (kappav1 * denom1 + kappav3 * denom2);

    
    g[0] = (f_local*(Ts0 - v1[0])/d1[0] + (1. - f_local)*(Ts0 - v1[1])/d1[1]) / Delta;
    g[1] = ((1. - f_local)*(Ts0 - v2[0])/d2[0] + f_local*(Ts0 - v2[1])/d2[1]) / Delta;
  }
}

foreach_dimension()
  void pair_gradient_23_x  (Point point, scalar * sl, scalar cs, scalar cs2, scalar cs3, face vector fs, face vector fs2,
       coord n, double alpha, coord p_1, coord p2_1, coord n2, double alpha2, coord p_2, coord p2_2, double g[2], bool homogenious) {


  double d1[2], v1[2] = {nodata, nodata};
  double d2[2], v2[2] = {nodata, nodata};
  for (int i = 0; i < 2; i++) {
    g[i] = nodata;
  }

  nearby_val2_x(point, sl[1], cs, cs2, cs3, fs, fs2, n, alpha, p_1, p2_1, v1, d1);

  nearby_val3_x(point, sl[2], cs, cs2, cs3, fs, fs2, n2, alpha2, p_2, p2_2, v2, d2);

  double FLUXv = homogenious ? 0 : FLUXFUN;
  if (v1[1] == nodata || v2[1] == nodata ) {
    //fprintf (fout, "test 3\n");
    double Ts0 = (v1[0]*(kappav2)/(Delta*d1[0]) + v2[0]*(kappav3)/(Delta*d2[0]) + FLUXv)/
      ((kappav2)/(Delta*d1[0]) + (kappav3)/(Delta*d2[0]));
    g[0]  = -(v1[0] - Ts0)/(d1[0]*Delta);
    g[1]  = -(v2[0] - Ts0)/(d2[0]*Delta);
  } else {

    double f_local = cs2[]; // local fraction value at the interface
    double denom1 = f_local/d1[0] + (1. - f_local)/d1[1];
    double denom2 = (1. - f_local)/d2[0] + f_local/d2[1];

    double num1 = f_local * v1[0]/d1[0] + (1. - f_local) * v1[1]/d1[1];
    double num2 = (1. - f_local) * v2[0]/d2[0] + f_local * v2[1]/d2[1];

    double Ts0 = (kappav2 * num1 + kappav3 * num2 + FLUXv) / (kappav2 * denom1 + kappav3 * denom2);

    
    g[0] = (f_local*(Ts0 - v1[0])/d1[0] + (1. - f_local)*(Ts0 - v1[1])/d1[1]) / Delta;
    g[1] = ((1. - f_local)*(Ts0 - v2[0])/d2[0] + f_local*(Ts0 - v2[1])/d2[1]) / Delta;
  }
}

double total_flux_embed (Point point, coord p, coord n);
#if EMBED
void embed_flux_three_field (Point point, scalar * s, double * val, bool homogenious, scalar f, scalar cs2, scalar cs3, face vector fs2)
{
  val[0] = 0.;
  val[1] = 0.;
  val[2] = 0.;

  if ( (cs[] >= 1. || cs[] <= 0.) && (cs2[] >= 1. || cs2[] <= 0.) && (cs3[] >= 1. || cs3[] <= 0.) )
    return;
  double g12[2], g13[2], g23[2];
  g12[0] = 0; g12[1] = 0;
  g13[0] = 0; g13[1] = 0;
  g23[0] = 0; g23[1] = 0;

  //pair 12
  if ((cs[] > 0) && (cs[] < 1) && (cs2[] > 0) && (cs2[] < 1)){
    if ((cs3[] >= 1. || cs3[] <= 0.)){ // non triple point

      coord n = facet_normal (point, cs, fs), p;
      foreach_dimension()
        n.x = -n.x;
      double alpha = plane_alpha (1-cs[], n);
      double area = plane_area_center (n, alpha, &p);

      coord p2;
      plane_center (n, alpha, 1-cs[], &p2);

      coord n2;
      foreach_dimension()
        n2.x = -n.x;


      coord p2_2;
      plane_center (n2, alpha, cs[], &p2_2);

      alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
      normalize (&n);
      normalize (&n2);


      #if dimension == 2
        foreach_dimension()
          if (fabs(n.x) >= fabs(n.y)) {
            pair_gradient_12_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious); 
          }

      #else // dimension == 3
        if (fabs(n.x) >= fabs(n.y)) {
          if (fabs(n.x) >= fabs(n.z)) {
            pair_gradient_12_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
          }
        }
        else if (fabs(n.y) >= fabs(n.z)) {
          pair_gradient_12_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
        } else {
          pair_gradient_12_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
        }
      #endif // dimension == 3

        double mua1 = kappav1, mua2 = kappav2;
        assert (g12[0] != nodata);
        assert (g12[1] != nodata);
        val[0] += -mua1*g12[0]*area/Delta;
        val[1] += -mua2*g12[1]*area/Delta;
    }
    else{

        coord n = facet_normal (point, cs, fs), p;
        foreach_dimension()
          n.x = -n.x;
        double alpha = plane_alpha (1-cs[], n);

        double area = plane_area_center (n, alpha, &p);

        coord p2;
        plane_center (n, alpha, 1-cs[], &p2);

        coord n2;
        foreach_dimension()
          n2.x = -n.x;

        //to approximate the position of the centroid of cs2 we perform a PLIC reconstruction
        coord ntemp = interface_normal (point, cs2);
        double alpha_temp = plane_alpha (cs2[], ntemp);

        coord p2_2;
        plane_center (ntemp, alpha_temp, cs2[], &p2_2);

        alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
        normalize (&n);
        normalize (&n2);

        coord n_temp = (coord){0,0,0};
        coord nn = interface_normal (point, f);

        #if dimension == 2
          if (fabs(nn.x) >= fabs(nn.y)){
            n_temp.x = -sign(nn.x);
          }
          else{
            n_temp.y = -sign(nn.y);
          }
        #else // dimension == 3
          if (fabs(nn.x) >= fabs(nn.y)) {
            if (fabs(nn.x) >= fabs(nn.z))
              n_temp.x = -sign(nn.x);
          }
          else{
            if (fabs(nn.y) >= fabs(nn.z)){
              n_temp.y = -sign(nn.y);
            }
            else{
              n_temp.z = -sign(nn.z);
            }
          }
        #endif // dimension == 3

        coord aa = (coord){-0.5,-0.5,-0.5};
        coord bb = (coord){0.5,0.5,0.5};

        foreach_dimension(){
          if (n_temp.x != 0){
            if (n_temp.x > 0){
              aa.x += (1-f[]);
            }
            else{
              bb.x -= (1-f[]);
            }
          }
        }

        area = compute_cross_section_area(n, alpha, aa, bb);

        #if dimension == 2
          foreach_dimension()
            if (fabs(n.x) >= fabs(n.y)) {
              pair_gradient_12_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious); 
            }

        #else // dimension == 3
          if (fabs(n.x) >= fabs(n.y)) {
            if (fabs(n.x) >= fabs(n.z)) {
              pair_gradient_12_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
            }
          }
          else if (fabs(n.y) >= fabs(n.z)) {
            pair_gradient_12_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
          } else {
            pair_gradient_12_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g12, homogenious);
          }
        #endif // dimension == 3

        ////area = 0.0;

        double mua1 = kappav1, mua2 = kappav2;
        assert (g12[0] != nodata);
        assert (g12[1] != nodata);
        val[0] += -mua1*g12[0]*area/Delta;
        val[1] += -mua2*g12[1]*area/Delta;
    }
  }

   //pair 13
  if ((cs[] > 0) && (cs[] < 1) && (cs3[] > 0) && (cs3[] < 1)){
    if ((cs2[] >= 1. || cs2[] <= 0.)){ // non triple point

      coord n = facet_normal (point, cs, fs), p;
      foreach_dimension()
        n.x = -n.x;
      double alpha = plane_alpha (1-cs[], n);
      double area = plane_area_center (n, alpha, &p);

      coord p2;
      plane_center (n, alpha, 1-cs[], &p2);

      coord n2;
      foreach_dimension()
        n2.x = -n.x;


      coord p2_2;
      plane_center (n2, alpha, cs[], &p2_2);

      alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
      normalize (&n);
      normalize (&n2);


      #if dimension == 2
        foreach_dimension()
          if (fabs(n.x) >= fabs(n.y)) {
            pair_gradient_13_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious); 
          }

      #else // dimension == 3
        if (fabs(n.x) >= fabs(n.y)) {
          if (fabs(n.x) >= fabs(n.z)) {
            pair_gradient_13_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
          }
        }
        else if (fabs(n.y) >= fabs(n.z)) {
          pair_gradient_13_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
        } else {
          pair_gradient_13_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
        }
      #endif // dimension == 3

        double mua1 = kappav1, mua2 = kappav3;
        assert (g13[0] != nodata);
        assert (g13[1] != nodata);
        val[0] += -mua1*g13[0]*area/Delta;
        val[2] += -mua2*g13[1]*area/Delta;
    }
    else{

        coord n = facet_normal (point, cs, fs), p;
        foreach_dimension()
          n.x = -n.x;
        double alpha = plane_alpha (1-cs[], n);

        double area = plane_area_center (n, alpha, &p);

        coord p2;
        plane_center (n, alpha, 1-cs[], &p2);

        coord n2;
        foreach_dimension()
          n2.x = -n.x;

        //to approximate the position of the centroid of cs2 we perform a PLIC reconstruction
        coord ntemp = interface_normal (point, cs3);
        double alpha_temp = plane_alpha (cs3[], ntemp);

        coord p2_2;
        plane_center (ntemp, alpha_temp, cs3[], &p2_2);

        alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
        normalize (&n);
        normalize (&n2);

        coord n_temp = (coord){0,0,0};
        coord nn = interface_normal (point, f);
        foreach_dimension()
          nn.x = -nn.x;

        #if dimension == 2
          if (fabs(nn.x) >= fabs(nn.y)){
            n_temp.x = -sign(nn.x);
          }
          else{
            n_temp.y = -sign(nn.y);
          }
        #else // dimension == 3
          if (fabs(nn.x) >= fabs(nn.y)) {
            if (fabs(nn.x) >= fabs(nn.z))
              n_temp.x = -sign(nn.x);
          }
          else{
            if (fabs(nn.y) >= fabs(nn.z)){
              n_temp.y = -sign(nn.y);
            }
            else{
              n_temp.z = -sign(nn.z);
            }
          }
        #endif // dimension == 3

        coord aa = (coord){-0.5,-0.5,-0.5};
        coord bb = (coord){0.5,0.5,0.5};

        foreach_dimension(){
          if (n_temp.x != 0){
            if (n_temp.x > 0){
              aa.x += (f[]);
            }
            else{
              bb.x -= (f[]);
            }
          }
        }

        area = compute_cross_section_area(n, alpha, aa, bb);

        #if dimension == 2
          foreach_dimension()
            if (fabs(n.x) >= fabs(n.y)) {
              pair_gradient_13_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious); 
            }

          #else // dimension == 3
            if (fabs(n.x) >= fabs(n.y)) {
              if (fabs(n.x) >= fabs(n.z)) {
                pair_gradient_13_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
              }
            }
            else if (fabs(n.y) >= fabs(n.z)) {
              pair_gradient_13_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
            } else {
              pair_gradient_13_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g13, homogenious);
            }
        #endif // dimension == 3

        ////area = 0.0;

        double mua1 = kappav1, mua2 = kappav3;
        assert (g13[0] != nodata);
        assert (g13[1] != nodata);
        val[0] += -mua1*g13[0]*area/Delta;
        val[2] += -mua2*g13[1]*area/Delta;
    }
  }

   //pair 23
  if ((cs2[] > 0) && (cs2[] < 1) && (cs3[] > 0) && (cs3[] < 1)){
    if ((cs[] >= 1. || cs[] <= 0.)){ // non triple point

      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      double area = plane_area_center (n, alpha, &p);

      coord p2;
      plane_center (n, alpha, f[], &p2);

      coord n2;
      foreach_dimension()
        n2.x = -n.x;

      coord p2_2;
      plane_center (n2, alpha, 1-f[], &p2_2);

      alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
      normalize (&n);
      normalize (&n2);


      #if dimension == 2
        foreach_dimension()
          if (fabs(n.x) >= fabs(n.y)) {
            pair_gradient_23_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious); 
          }

      #else // dimension == 3
        if (fabs(n.x) >= fabs(n.y)) {
          if (fabs(n.x) >= fabs(n.z)) {
            pair_gradient_23_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
          }
        }
        else if (fabs(n.y) >= fabs(n.z)) {
          pair_gradient_23_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
        } else {
          pair_gradient_23_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
        }
      #endif // dimension == 3

        double mua1 = kappav2, mua2 = kappav3;
        assert (g12[0] != nodata);
        assert (g12[1] != nodata);
        val[1] += -mua1*g23[0]*area/Delta;
        val[2] += -mua2*g23[1]*area/Delta;
    }
    else{

        coord n = interface_normal (point, f), p;
        double alpha = plane_alpha (f[], n);
        double area = plane_area_center (n, alpha, &p);

        coord ntemp = interface_normal (point, cs2);
        double alpha_temp = plane_alpha (cs2[], ntemp);

        coord p2;
        plane_center (ntemp, alpha_temp, cs2[], &p2);

        ntemp = interface_normal (point, cs3);
        alpha_temp = plane_alpha (cs3[], ntemp);

        coord p2_2;
        plane_center (ntemp, alpha_temp, cs3[], &p2_2);

        coord n2;
        foreach_dimension()
          n2.x = -n.x;

        alpha /= sqrt(sq(n.x) + sq(n.y) + sq(n.z));
        normalize (&n);
        normalize (&n2);

        coord n_temp = (coord){0,0,0};
        coord nn = facet_normal (point, cs, fs);

        #if dimension == 2
          if (fabs(nn.x) >= fabs(nn.y)){
            n_temp.x = -sign(nn.x);
          }
          else{
            n_temp.y = -sign(nn.y);
          }
        #else // dimension == 3
          if (fabs(nn.x) >= fabs(nn.y)) {
            if (fabs(nn.x) >= fabs(nn.z))
              n_temp.x = -sign(nn.x);
          }
          else{
            if (fabs(nn.y) >= fabs(nn.z)){
              n_temp.y = -sign(nn.y);
            }
            else{
              n_temp.z = -sign(nn.z);
            }
          }
        #endif // dimension == 3

        coord aa = (coord){-0.5,-0.5,-0.5};
        coord bb = (coord){0.5,0.5,0.5};

        foreach_dimension(){
          if (n_temp.x != 0){
            if (n_temp.x > 0){
              aa.x += (1-cs[]);
            }
            else{
              bb.x -= (1-cs[]);
            }
          }
        }

        area = compute_cross_section_area(n, alpha, aa, bb);

        #if dimension == 2
          foreach_dimension()
            if (fabs(n.x) >= fabs(n.y)) {
              pair_gradient_23_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious); 
            }

        #else // dimension == 3
          if (fabs(n.x) >= fabs(n.y)) {
            if (fabs(n.x) >= fabs(n.z)) {
              pair_gradient_23_x  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
            }
          }
          else if (fabs(n.y) >= fabs(n.z)) {
            pair_gradient_23_y  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
          } else {
            pair_gradient_23_z  (point, s, cs, cs2, cs3, fs, fs2, n, alpha, p, p2, n2, alpha, p, p2_2, g23, homogenious);
          }
        #endif // dimension == 3

        //area = 0.0;


        double mua1 = kappav2, mua2 = kappav3;
        assert (g23[0] != nodata);
        assert (g23[1] != nodata);
        val[1] += -mua1*g23[0]*area/Delta;
        val[2] += -mua2*g23[1]*area/Delta;
    }
  }

}
#endif

/**
## Tests
* [A simple Poisson problem split on to two domains](test-pair.c)

## Usage
* [An implicit conjugate heatflux diffusion solver](diffusion-pair.h)
 */

