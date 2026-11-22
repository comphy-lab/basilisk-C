/**
# Integral formulation for surface tension

See [Al Saud et al., 2018](#alsaud2018), [Popinet & Zaleski,
1999](#popinet1999), and [Popinet, 2018](#popinet2018) for details.

This header extends the 2D integral formulation to 1D, 2D, axisymmetric
(`AXI`), and 3D meshes (uniform multigrid, adaptive quadtrees/octrees,
multi-root bricks, and MPI) using $3^d$ (`BGHOSTS = 1`) stencils:

* **Constant $\sigma$ (`is_constant (d.sigmaf)` with
  `INTEGRAL_TENSOR == 0`)**: Uses level-set curvature projected onto
  $d = 0$ via the parallel-surface Weingarten relation with the
  well-balanced CSF force $\alpha_f \sigma \kappa_f \nabla_h f$.
* **Variable $\sigma(\mathbf{x})$ (or `#define INTEGRAL_TENSOR 1`)**:
  Uses the conservative capillary stress tensor
  $\mathbf{S} = \sigma (\mathbf{I} - \mathbf{n}\otimes\mathbf{n})
  \delta_\Gamma$ across staggered momentum control volumes (exact 2D
  edge cuts and 3D 4-triangle face cuts).

The surface tension field $\sigma$ will be associated to each levelset
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */

#ifndef INTEGRAL_TENSOR
# define INTEGRAL_TENSOR 0
#endif

extern scalar * tracers;

attribute {
  scalar sigmaf;
}

/**
Surface tension is a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {
  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave:
$$
T = \sqrt{\frac{\rho_m \Delta_{\min}^3}{\pi \sigma}}
$$
with $\rho_m = (\rho_1 + \rho_2)/2$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++)
{
  double sigma = 0.;
  bool has_var = false;
  for (scalar d in tracers)
    if (d.sigmaf.i) {
      if (is_constant (d.sigmaf))
        sigma += constant (d.sigmaf);
      else
        has_var = true;
    }
  double sigmamax = sigma;

  if (sigmamax || has_var) {

    /**
    We first compute the minimum and maximum values of $\alpha/f_m =
    1/\rho$, as well as $\Delta_{\min}$. */

    double amin = HUGE, amax = -HUGE, dmin = HUGE;
    foreach_face (reduction(min:amin) reduction(max:amax)
                  reduction(min:dmin) reduction(max:sigmamax))
      if (fm.x[] > 0.) {
        if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
        if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
        if (Delta < dmin) dmin = Delta;

        /**
        The maximum timestep is set using the sum of surface tension
        coefficients. */

        if (has_var) {
          double sigmai = sigma;
          for (scalar d in tracers)
            if (d.sigmaf.i && !is_constant (d.sigmaf) &&
                fabs (d[]) < 2.5*Delta) {
              scalar sigmaf = d.sigmaf;
              sigmai += sigmaf[];
            }
          if (sigmai > sigmamax)
            sigmamax = sigmai;
        }
      }
    double rhom = (1./amin + 1./amax)/2.;

    if (sigmamax > 0.) {
      double dt = sqrt (rhom*cube(dmin)/(pi*sigmamax));
      if (dt < dtmax)
        dtmax = dt;
    }
  }
}

/**
## Tree prolongation and restriction for `nodata`-masked curvature */

#if TREE
static void distance_curvature_restriction (Point point, scalar kappa)
{
  double k = 0., s = 0.;
  foreach_child()
    if (kappa[] != nodata)
      k += kappa[], s++;
  kappa[] = s ? k/s : nodata;
}

static void distance_curvature_prolongation (Point point, scalar kappa)
{
  foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)
# if dimension > 1
      for (int j = 0; j <= 1; j++)
# endif
# if dimension > 2
        for (int k = 0; k <= 1; k++)
# endif
          if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
            sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    kappa[] = s ? sk/s : nodata;
  }
}
#endif // TREE

/**
## Curvature

We compute the total curvature $\kappa = -\nabla \cdot (\nabla d /
|\nabla d|)$ from the levelset function *d* and project it onto the
zero level set $d = 0$ via the parallel-surface Weingarten relation so
that $p_{d > 0} - p_{d < 0} = +\sigma \kappa$. */

#ifndef PARALLEL_SURFACE
# define PARALLEL_SURFACE 1
#endif

#if dimension == 2
static inline double distance_curvature_plane (Point point, scalar d)
{
  double dx = (d[1] - d[-1])/(2.*Delta);
  double dy = (d[0,1] - d[0,-1])/(2.*Delta);
  double dn2 = sq(dx) + sq(dy);
  if (dn2 <= 1e-12)
    return nodata;
  double dn = sqrt (dn2);

  double dxx = (d[1] - 2.*d[] + d[-1])/Delta;
  double dyy = (d[0,1] - 2.*d[] + d[0,-1])/Delta;
  double dxy = (d[1,1] - d[-1,1] - d[1,-1] + d[-1,-1])/(4.*Delta);
  double kappa = - (sq(dx)*dyy - 2.*dx*dy*dxy + sq(dy)*dxx)/cube(dn)/Delta;

# if PARALLEL_SURFACE
  double dperp = d[]/dn;
  if (fabs (dperp) < 2.5*Delta) {
    double denom = 1. + kappa*dperp;
    if (denom > 0.5 && denom < 2.0)
      kappa /= denom;
  }
# endif

  if (fabs (kappa) > 1./Delta)
    kappa = (kappa > 0. ? 1. : -1.)/Delta;
  return kappa;
}
#endif // dimension == 2

static inline double distance_curvature (Point point, scalar d)
{
#if dimension == 1
  return 0. [-1];

#elif dimension == 2
  double kappa = distance_curvature_plane (point, d);
  if (kappa == nodata)
    return nodata;
# if AXI
  double dx = (d[1] - d[-1])/(2.*Delta);
  double dy = (d[0,1] - d[0,-1])/(2.*Delta);
  double dn = sqrt (sq(dx) + sq(dy));
  double r = max (y, Delta/2.);
  double ktheta = - dy/(dn*r);
#  if PARALLEL_SURFACE
  double dperp = d[]/dn;
  if (fabs (dperp) < 2.5*Delta) {
    double denom_t = 1. + ktheta*dperp;
    if (denom_t > 0.5 && denom_t < 2.0)
      ktheta /= denom_t;
  }
#  endif
  kappa += ktheta;
  if (fabs (kappa) > 2./Delta)
    kappa = (kappa > 0. ? 2. : -2.)/Delta;
# endif
  return kappa;

#elif dimension == 3
  double dx = (d[1] - d[-1])/(2.*Delta);
  double dy = (d[0,1] - d[0,-1])/(2.*Delta);
  double dz = (d[0,0,1] - d[0,0,-1])/(2.*Delta);
  double dn2 = sq(dx) + sq(dy) + sq(dz);
  if (dn2 <= 1e-12)
    return nodata;
  double dn = sqrt (dn2);
  double nx = dx/dn, ny = dy/dn, nz = dz/dn;

  double hxx = (d[1] - 2.*d[] + d[-1])/(dn*sq(Delta));
  double hyy = (d[0,1] - 2.*d[] + d[0,-1])/(dn*sq(Delta));
  double hzz = (d[0,0,1] - 2.*d[] + d[0,0,-1])/(dn*sq(Delta));
  double hxy = (d[1,1,0] - d[-1,1,0] -
                d[1,-1,0] + d[-1,-1,0])/(4.*dn*sq(Delta));
  double hxz = (d[1,0,1] - d[-1,0,1] -
                d[1,0,-1] + d[-1,0,-1])/(4.*dn*sq(Delta));
  double hyz = (d[0,1,1] - d[0,-1,1] -
                d[0,1,-1] + d[0,-1,-1])/(4.*dn*sq(Delta));

  double vx = hxx*nx + hxy*ny + hxz*nz;
  double vy = hxy*nx + hyy*ny + hyz*nz;
  double vz = hxz*nx + hyz*ny + hzz*nz;
  double nhn = nx*vx + ny*vy + nz*vz;

  double kappa = - (hxx + hyy + hzz - nhn);

# if PARALLEL_SURFACE
  double dperp = d[]/dn;
  if (fabs (dperp) < 2.5*Delta) {
    double hnorm2 = sq(hxx) + sq(hyy) + sq(hzz) +
      2.*(sq(hxy) + sq(hxz) + sq(hyz));
    double trW2 = hnorm2 - 2.*(sq(vx) + sq(vy) + sq(vz)) + sq(nhn);
    double gaussK = 0.5*(sq(kappa) - trW2);
    double denom = 1. + kappa*dperp + gaussK*sq(dperp);
    if (denom > 0.5 && denom < 2.0)
      kappa = (kappa + 2.*gaussK*dperp)/denom;
  }
# endif

  if (fabs (kappa) > 2./Delta)
    kappa = (kappa > 0. ? 2. : -2.)/Delta;
  return kappa;
#endif
}

static void compute_curvature (scalar d, scalar kappa, bool plane)
{
#if TREE
  kappa.refine = kappa.prolongation = distance_curvature_prolongation;
  kappa.restriction = distance_curvature_restriction;
#endif
  foreach() {
    if (fabs (d[]) <= 2.5*Delta) {
#if dimension == 2
      kappa[] = plane ? distance_curvature_plane (point, d) :
        distance_curvature (point, d);
#else
      (void) plane;
      kappa[] = distance_curvature (point, d);
#endif
    }
    else
      kappa[] = nodata;
  }
}

/**
## Staggered edge shear (2D) and face cut geometry (3D)

In 2D, the off-diagonal capillary shear on each staggered control-volume
edge is evaluated inline from the $3^2$ stencil. In 3D, each square face
$[0, 1]^2$ bounding a staggered momentum control volume is split
symmetrically into 4 triangles around `(0.5, 0.5)` so interface loops
close identically across all 6 faces, giving machine-zero force on
tilted planes and exact discrete momentum conservation. */

#if dimension == 2
static inline double edge_shear_2d (double dm0, double dm1,
                                    double dp0, double dp1,
                                    double sm0, double sm1,
                                    double sp0, double sp1,
                                    double Delta)
{
  double dm = dm1 + dm0, dp = dp1 + dp0;
  if (dp*dm > 0. || dm == dp)
    return 0.*sm0/Delta;
  double xi = dm/(dm - dp);
  double gy = (xi*(dp1 - dm1 + dm0 - dp0) + dm1 - dm0)/Delta;
  double gx = (dp - dm)/(2.*Delta);
  double dn = sqrt (sq(gx) + sq(gy));
  double ny = dn > 0. ? gy/dn : 0.;
  double sigmai = (sm1 + sm0 + xi*(sp1 + sp0 - sm1 - sm0))/2.;
  double sgn = dp > 0. ? 1. : (dp < 0. ? -1. : (dm < 0. ? 1. : -1.));
  return - sigmai*sgn*ny/Delta;
}
#elif dimension == 3
static inline int tri_iso_seg (double u0, double v0, double f0, double g0,
                               double u1, double v1, double f1, double g1,
                               double u2, double v2, double f2, double g2,
                               double * uLo, double * vLo,
                               double * uPAo, double * vPAo, double * gPAo,
                               double * uPBo, double * vPBo, double * gPBo,
                               double * tAo, double * tBo, bool * pLo)
{
  bool p0 = (f0 >= 0.), p1 = (f1 >= 0.), p2 = (f2 >= 0.);
  int npos = (int) p0 + (int) p1 + (int) p2;
  if (npos == 3)
    return 1;
  if (npos == 0)
    return 0;

  double uL, vL, fL, gL, uA, vA, fA, gA, uB, vB, fB, gB;
  bool pL;
  if (p0 != p1 && p1 == p2) {
    uL = u0; vL = v0; fL = f0; gL = g0; pL = p0;
    uA = u1; vA = v1; fA = f1; gA = g1;
    uB = u2; vB = v2; fB = f2; gB = g2;
  }
  else if (p1 != p2 && p2 == p0) {
    uL = u1; vL = v1; fL = f1; gL = g1; pL = p1;
    uA = u2; vA = v2; fA = f2; gA = g2;
    uB = u0; vB = v0; fB = f0; gB = g0;
  }
  else {
    uL = u2; vL = v2; fL = f2; gL = g2; pL = p2;
    uA = u0; vA = v0; fA = f0; gA = g0;
    uB = u1; vB = v1; fB = f1; gB = g1;
  }

  double tA = fL/(fL - fA), tB = fL/(fL - fB);
  *uLo  = uL; *vLo  = vL; *pLo  = pL;
  *tAo  = tA; *tBo  = tB;
  *uPAo = uL + tA*(uA - uL);
  *vPAo = vL + tA*(vA - vL);
  *gPAo = gL + tA*(gA - gL);
  *uPBo = uL + tB*(uB - uL);
  *vPBo = vL + tB*(vB - vL);
  *gPBo = gL + tB*(gB - gL);
  return 2;
}

static inline void tri_cut_2d (double u0, double v0, double f0, double g0,
                               double u1, double v1, double f1, double g1,
                               double u2, double v2, double f2, double g2,
                               double * L_tot, double * su_tot,
                               double * sv_tot, double * ssig_tot,
                               double * area_pos)
{
  double uL = 0., vL = 0., uPA = 0., vPA = 0., gPA = 0.*g0;
  double uPB = 0., vPB = 0., gPB = 0.*g0, tA = 0., tB = 0.;
  bool pL = false;
  int st = tri_iso_seg (u0, v0, f0, g0, u1, v1, f1, g1, u2, v2, f2, g2,
                        &uL, &vL, &uPA, &vPA, &gPA, &uPB, &vPB, &gPB,
                        &tA, &tB, &pL);
  if (st == 1) {
    *area_pos += 0.25;
    return;
  }
  if (st == 0)
    return;

  double seg = sqrt (sq(uPA - uPB) + sq(vPA - vPB));
  *L_tot    += seg;
  *su_tot   += seg*0.5*(uPA + uPB);
  *sv_tot   += seg*0.5*(vPA + vPB);
  *ssig_tot += seg*0.5*(gPA + gPB);

  double tri_frac = 0.25*tA*tB;
  *area_pos += pL ? tri_frac : (0.25 - tri_frac);
}

static inline double face_cut_2d (double d00, double d10,
                                  double d01, double d11,
                                  double s00, double s10,
                                  double s01, double s11,
                                  double * um, double * vm,
                                  double * sig_m, double * area_pos)
{
  double dc = 0.25*(d00 + d10 + d01 + d11);
  double sc = 0.25*(s00 + s10 + s01 + s11);
  double L_tot = 0., su_tot = 0., sv_tot = 0., ssig_tot = 0.*sc, apos = 0.;

  tri_cut_2d (0., 0., d00, s00, 1., 0., d10, s10, 0.5, 0.5, dc, sc,
              &L_tot, &su_tot, &sv_tot, &ssig_tot, &apos);
  tri_cut_2d (1., 0., d10, s10, 1., 1., d11, s11, 0.5, 0.5, dc, sc,
              &L_tot, &su_tot, &sv_tot, &ssig_tot, &apos);
  tri_cut_2d (1., 1., d11, s11, 0., 1., d01, s01, 0.5, 0.5, dc, sc,
              &L_tot, &su_tot, &sv_tot, &ssig_tot, &apos);
  tri_cut_2d (0., 1., d01, s01, 0., 0., d00, s00, 0.5, 0.5, dc, sc,
              &L_tot, &su_tot, &sv_tot, &ssig_tot, &apos);

  if (L_tot > 1e-14) {
    *um    = su_tot/L_tot;
    *vm    = sv_tot/L_tot;
    *sig_m = ssig_tot/L_tot;
  }
  else {
    *um    = 0.5;
    *vm    = 0.5;
    *sig_m = sc;
  }
  *area_pos = apos;
  return L_tot;
}

static inline double tri_shear_3d (double u0, double v0,
                                   double f0, double g0,
                                   double u1, double v1,
                                   double f1, double g1,
                                   double u2, double v2,
                                   double f2, double g2,
                                   double d00, double d10,
                                   double d01, double d11,
                                   double gb00, double gb10,
                                   double gb01, double gb11,
                                   double Delta)
{
  double uL = 0., vL = 0., uPA = 0., vPA = 0., gPA = 0.*g0;
  double uPB = 0., vPB = 0., gPB = 0.*g0, tA = 0., tB = 0.;
  bool pL = false;
  if (tri_iso_seg (u0, v0, f0, g0, u1, v1, f1, g1, u2, v2, f2, g2,
                   &uL, &vL, &uPA, &vPA, &gPA, &uPB, &vPB, &gPB,
                   &tA, &tB, &pL) != 2)
    return 0.*g0/Delta;

  double wa = vPB - vPA;
  double wc = -(uPB - uPA);
  if (wa == 0. && wc == 0.)
    return 0.*g0/Delta;
  double dot = wa*(uL - uPA) + wc*(vL - vPA);
  if (pL ? (dot < 0.) : (dot > 0.))
    wa = - wa;

  double um = 0.5*(uPA + uPB), vm = 0.5*(vPA + vPB);
  double sig_m = 0.5*(gPA + gPB);
  double ga_m = ((1. - vm)*(d10 - d00) + vm*(d11 - d01))/Delta;
  double gc_m = ((1. - um)*(d01 - d00) + um*(d11 - d10))/Delta;
  double gb_m = (1. - um)*(1. - vm)*gb00 + um*(1. - vm)*gb10 +
    (1. - um)*vm*gb01 + um*vm*gb11;
  double dn = sqrt (sq(ga_m) + sq(gc_m) + sq(gb_m));
  return dn > 0. ? - sig_m*(gb_m/dn)*wa/Delta : 0.*g0/Delta;
}

static inline double edge_shear_3d (double dm00, double dm10,
                                    double dm0m, double dm1m,
                                    double dm0p, double dm1p,
                                    double dp00, double dp10,
                                    double dp0m, double dp1m,
                                    double dp0p, double dp1p,
                                    double sm00, double sm10,
                                    double sm0m, double sm1m,
                                    double sm0p, double sm1p,
                                    double sp00, double sp10,
                                    double sp0m, double sp1m,
                                    double sp0p, double sp1p,
                                    double Delta)
{
  double d00 = 0.25*(dm00 + dm10 + dm0m + dm1m);
  double d10 = 0.25*(dp00 + dp10 + dp0m + dp1m);
  double d01 = 0.25*(dm00 + dm10 + dm0p + dm1p);
  double d11 = 0.25*(dp00 + dp10 + dp0p + dp1p);
  bool p0 = (d00 >= 0.);
  if (p0 == (d10 >= 0.) && p0 == (d01 >= 0.) && p0 == (d11 >= 0.))
    return 0.*sm00/Delta;

  double s00 = 0.25*(sm00 + sm10 + sm0m + sm1m);
  double s10 = 0.25*(sp00 + sp10 + sp0m + sp1m);
  double s01 = 0.25*(sm00 + sm10 + sm0p + sm1p);
  double s11 = 0.25*(sp00 + sp10 + sp0p + sp1p);
  double dc = 0.25*(d00 + d10 + d01 + d11);
  double sc = 0.25*(s00 + s10 + s01 + s11);

  double gb00 = 0.5*(dm10 - dm00 + dm1m - dm0m)/Delta;
  double gb10 = 0.5*(dp10 - dp00 + dp1m - dp0m)/Delta;
  double gb01 = 0.5*(dm10 - dm00 + dm1p - dm0p)/Delta;
  double gb11 = 0.5*(dp10 - dp00 + dp1p - dp0p)/Delta;

  return tri_shear_3d (0., 0., d00, s00, 1., 0., d10, s10, 0.5, 0.5, dc, sc,
                       d00, d10, d01, d11, gb00, gb10, gb01, gb11, Delta) +
    tri_shear_3d (1., 0., d10, s10, 1., 1., d11, s11, 0.5, 0.5, dc, sc,
                  d00, d10, d01, d11, gb00, gb10, gb01, gb11, Delta) +
    tri_shear_3d (1., 1., d11, s11, 0., 1., d01, s01, 0.5, 0.5, dc, sc,
                  d00, d10, d01, d11, gb00, gb10, gb01, gb11, Delta) +
    tri_shear_3d (0., 1., d01, s01, 0., 0., d00, s00, 0.5, 0.5, dc, sc,
                  d00, d10, d01, d11, gb00, gb10, gb01, gb11, Delta);
}
#endif // dimension == 3

/**
## Surface tension term

The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

#if !INTEGRAL_TENSOR || AXI
# include "fractions.h"
#endif

event acceleration (i++)
{

  /**
  We check whether surface tension is associated with any levelset
  function *d*. */

  for (scalar d in tracers)
    if (d.sigmaf.i) {
      (const) scalar sigma = d.sigmaf;

#if !INTEGRAL_TENSOR
      if (is_constant (sigma)) {

        /**
        For constant surface tension we use the well-balanced CSF
        formulation with parallel-surface level-set curvature. */

        extern scalar f;

        scalar kappa[];
        compute_curvature (d, kappa, false);
        foreach()
          f[] = clamp (f[], 0., 1.);

# if TREE
        set_prolongation (f, p.prolongation);
# endif

        face vector av = a;
        foreach_face()
          if (f[] != f[-1] && fm.x[] > 0.) {
            double kappaf =
              (kappa[] < nodata && kappa[-1] < nodata) ?
              (kappa[] + kappa[-1])/2. :
              kappa[] < nodata ? kappa[] :
              kappa[-1] < nodata ? kappa[-1] :
              0.;
            av.x[] += alpha.x[]/(fm.x[] + SEPS)*sigma[]*kappaf*
              (f[] - f[-1])/Delta;
          }

# if TREE
        set_prolongation (f, fraction_refine);
# endif
        continue;
      }
#endif // !INTEGRAL_TENSOR

      /**
      For variable surface tension (or when `INTEGRAL_TENSOR` is set)
      we compute the divergence of the capillary stress tensor. We
      allocate the diagonal stress components manually as plain scalars
      because the symmetries of default vector/tensor allocations in
      Basilisk are antisymmetric on normal boundaries. */

#if dimension == 1
      (void) sigma;

#elif dimension == 2
      scalar kappa[];
      compute_curvature (d, kappa, true);

      scalar Sxx[], Syy[];
      vector Sdiag; Sdiag.x = Sxx, Sdiag.y = Syy;
      foreach()
        foreach_dimension() {
          Sdiag.y[] = 0.*sigma[]/Delta;
          for (int i = -1; i <= 1; i += 2)
            if (d[]*(d[] + d[i]) < 0.) {
              double xi = d[]/(d[] - d[i]);
              double gx = ((d[1] - d[-1])/2. +
                           xi*i*(d[-1] - 2.*d[] + d[1]))/Delta;
              double gy0 = (d[0,1] - d[0,-1])/(2.*Delta);
              double gyi = (d[i,1] - d[i,-1])/(2.*Delta);
              double gy = gy0 + xi*(gyi - gy0);
              double dn = sqrt (sq(gx) + sq(gy));
              double nx = dn > 0. ? fabs (gx)/dn : 0.;
              double sigmai = sigma[] + xi*(sigma[i] - sigma[]);
              double ki = kappa[] < nodata ? kappa[] : 0.;
              Sdiag.y[] += sigmai*(nx/Delta +
                                   (d[] >= 0. ? 1. : -1.)*ki*(0.5 - xi));
            }
        }

      face vector av = a;
      foreach_face() {
        double sxy0 = edge_shear_2d (d[-1,-1], d[-1,0], d[0,-1], d[0,0],
                                     sigma[-1,-1], sigma[-1,0],
                                     sigma[0,-1], sigma[0,0], Delta);
        double sxy1 = edge_shear_2d (d[-1,0], d[-1,1], d[0,0], d[0,1],
                                     sigma[-1,0], sigma[-1,1],
                                     sigma[0,0], sigma[0,1], Delta);
        av.x[] += alpha.x[]/(fm.x[] + SEPS)*
          (Sdiag.x[] - Sdiag.x[-1] + sxy1 - sxy0)/Delta;

        /**
        ### Axisymmetric term

        The axisymmetric acceleration is computed using the volumetric
        surface tension formulation as $\frac{1}{\rho}\sigma
        \kappa_\theta \mathbf{n}\delta_s$ with $\kappa_\theta = n^r/r$
        and $\mathbf{n}\delta_s \approx \nabla f$. */

# if AXI
        coord grad_d = {
          (d[] - d[-1])/Delta,
          (d[0,1] + d[-1,1] - d[0,-1] - d[-1,-1])/(4.*Delta)
        };
        double dn = sqrt (sq(grad_d.x) + sq(grad_d.y));
        extern scalar f;
        if (dn > 0.)
          av.x[] -= alpha.x[]/sq(fm.x[] + SEPS)*
            (sigma[] + sigma[-1])/2.*(grad_d.y/dn)*(f[] - f[-1])/Delta;
# endif // AXI
      }

#elif dimension == 3
      scalar kappa[];
      compute_curvature (d, kappa, false);

      scalar Sxx[], Syy[], Szz[];
      vector Sdiag; Sdiag.x = Sxx, Sdiag.y = Syy, Sdiag.z = Szz;
      foreach()
        foreach_dimension() {
          Sdiag.x[] = 0.*sigma[]/Delta;
          double d00 = 0.25*(d[0,0,0] + d[0,-1,0] + d[0,0,-1] + d[0,-1,-1]);
          double d10 = 0.25*(d[0,0,0] + d[0,1,0] + d[0,0,-1] + d[0,1,-1]);
          double d01 = 0.25*(d[0,0,0] + d[0,-1,0] + d[0,0,1] + d[0,-1,1]);
          double d11 = 0.25*(d[0,0,0] + d[0,1,0] + d[0,0,1] + d[0,1,1]);
          bool p0 = (d00 >= 0.);
          if (p0 != (d10 >= 0.) || p0 != (d01 >= 0.) || p0 != (d11 >= 0.)) {
            double s00 = 0.25*(sigma[0,0,0] + sigma[0,-1,0] +
                               sigma[0,0,-1] + sigma[0,-1,-1]);
            double s10 = 0.25*(sigma[0,0,0] + sigma[0,1,0] +
                               sigma[0,0,-1] + sigma[0,1,-1]);
            double s01 = 0.25*(sigma[0,0,0] + sigma[0,-1,0] +
                               sigma[0,0,1] + sigma[0,-1,1]);
            double s11 = 0.25*(sigma[0,0,0] + sigma[0,1,0] +
                               sigma[0,0,1] + sigma[0,1,1]);
            double um, vm, sig_m, area_pos;
            double L_seg = face_cut_2d (d00, d10, d01, d11,
                                        s00, s10, s01, s11,
                                        &um, &vm, &sig_m, &area_pos);
            if (L_seg > 0.) {
              double gb_m = ((1. - vm)*(d10 - d00) + vm*(d11 - d01))/Delta;
              double gc_m = ((1. - um)*(d01 - d00) + um*(d11 - d10))/Delta;
              double ga00 = 0.125*(d[1,0,0] - d[-1,0,0] +
                                   d[1,-1,0] - d[-1,-1,0] +
                                   d[1,0,-1] - d[-1,0,-1] +
                                   d[1,-1,-1] - d[-1,-1,-1])/Delta;
              double ga10 = 0.125*(d[1,0,0] - d[-1,0,0] +
                                   d[1,1,0] - d[-1,1,0] +
                                   d[1,0,-1] - d[-1,0,-1] +
                                   d[1,1,-1] - d[-1,1,-1])/Delta;
              double ga01 = 0.125*(d[1,0,0] - d[-1,0,0] +
                                   d[1,-1,0] - d[-1,-1,0] +
                                   d[1,0,1] - d[-1,0,1] +
                                   d[1,-1,1] - d[-1,-1,1])/Delta;
              double ga11 = 0.125*(d[1,0,0] - d[-1,0,0] +
                                   d[1,1,0] - d[-1,1,0] +
                                   d[1,0,1] - d[-1,0,1] +
                                   d[1,1,1] - d[-1,1,1])/Delta;
              double ga_m = (1. - um)*(1. - vm)*ga00 +
                um*(1. - vm)*ga10 + (1. - um)*vm*ga01 + um*vm*ga11;
              double n_bc = sqrt (sq(gb_m) + sq(gc_m));
              double dn = sqrt (sq(ga_m) + sq(n_bc));
              double conorm_a = dn > 0. ? n_bc/dn : 0.;

              double frac_other = d[] >= 0. ? 1. - area_pos : area_pos;
              double s_sign = d[] >= 0. ? 1. : -1.;
              double ki = kappa[] < nodata ? kappa[] : 0.;
              Sdiag.x[] = sig_m*conorm_a*L_seg/Delta +
                s_sign*sig_m*ki*frac_other;
            }
          }
        }

      face vector av = a;
      foreach_face()
        if (fm.x[] > 0.) {
          double sxy_m = edge_shear_3d (
            d[-1,-1,0], d[-1,0,0], d[-1,-1,-1], d[-1,0,-1],
            d[-1,-1,1], d[-1,0,1],
            d[0,-1,0], d[0,0,0], d[0,-1,-1], d[0,0,-1],
            d[0,-1,1], d[0,0,1],
            sigma[-1,-1,0], sigma[-1,0,0], sigma[-1,-1,-1], sigma[-1,0,-1],
            sigma[-1,-1,1], sigma[-1,0,1],
            sigma[0,-1,0], sigma[0,0,0], sigma[0,-1,-1], sigma[0,0,-1],
            sigma[0,-1,1], sigma[0,0,1], Delta);
          double sxy_p = edge_shear_3d (
            d[-1,0,0], d[-1,1,0], d[-1,0,-1], d[-1,1,-1],
            d[-1,0,1], d[-1,1,1],
            d[0,0,0], d[0,1,0], d[0,0,-1], d[0,1,-1],
            d[0,0,1], d[0,1,1],
            sigma[-1,0,0], sigma[-1,1,0], sigma[-1,0,-1], sigma[-1,1,-1],
            sigma[-1,0,1], sigma[-1,1,1],
            sigma[0,0,0], sigma[0,1,0], sigma[0,0,-1], sigma[0,1,-1],
            sigma[0,0,1], sigma[0,1,1], Delta);
          double sxz_m = edge_shear_3d (
            d[-1,0,-1], d[-1,0,0], d[-1,-1,-1], d[-1,-1,0],
            d[-1,1,-1], d[-1,1,0],
            d[0,0,-1], d[0,0,0], d[0,-1,-1], d[0,-1,0],
            d[0,1,-1], d[0,1,0],
            sigma[-1,0,-1], sigma[-1,0,0], sigma[-1,-1,-1], sigma[-1,-1,0],
            sigma[-1,1,-1], sigma[-1,1,0],
            sigma[0,0,-1], sigma[0,0,0], sigma[0,-1,-1], sigma[0,-1,0],
            sigma[0,1,-1], sigma[0,1,0], Delta);
          double sxz_p = edge_shear_3d (
            d[-1,0,0], d[-1,0,1], d[-1,-1,0], d[-1,-1,1],
            d[-1,1,0], d[-1,1,1],
            d[0,0,0], d[0,0,1], d[0,-1,0], d[0,-1,1],
            d[0,1,0], d[0,1,1],
            sigma[-1,0,0], sigma[-1,0,1], sigma[-1,-1,0], sigma[-1,-1,1],
            sigma[-1,1,0], sigma[-1,1,1],
            sigma[0,0,0], sigma[0,0,1], sigma[0,-1,0], sigma[0,-1,1],
            sigma[0,1,0], sigma[0,1,1], Delta);
          av.x[] += alpha.x[]/(fm.x[] + SEPS)*
            (Sdiag.x[] - Sdiag.x[-1] + sxy_p - sxy_m + sxz_p - sxz_m)/Delta;
        }
#endif // dimension == 3
    }
}

/**
## References

~~~bib
@hal{alsaud2018, hal-01706565}

@hal{popinet2018, hal-01528255}

@article{popinet1999,
  title={A front-tracking algorithm for accurate representation of surface tension},
  author={S. Popinet and S. Zaleski},
  journal={International Journal for Numerical Methods in Fluids},
  volume={30},
  number={6},
  pages={775--793},
  year={1999},
  publisher={Wiley Online Library}
}
~~~
*/
