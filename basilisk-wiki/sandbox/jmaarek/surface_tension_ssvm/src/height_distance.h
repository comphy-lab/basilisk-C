#include "heights.h"
#include "redistance.h"

/**
## Determinism helper

`slp_barrier()` routes a value through a `volatile` temporary so the
compiler cannot fold, reassociate, or SLP-vectorize it together with
neighbouring operations (the standard requires volatile accesses to
occur exactly as written). Several short, fixed-length accumulations
below (squared-normal norms, weighted averages of `qx`/`qy`/`qz`)
feed directly into the returned closest-point vector `p`, which
becomes the signed distance `dist[]` -- and `dist[]` in turn gates a
hard threshold comparison in `viscosity_st.h`. Reordering by SLP
auto-vectorization can therefore flip that branch between otherwise
equivalent compiler versions, producing large result differences
rather than small rounding noise. Wrapping the sensitive
accumulations fixes their evaluation order regardless of compiler or
optimisation level. Guarded so this header and `viscosity_st.h` can
both be included without a redefinition error, whichever is included
first.

## Compiler-dependence: what was actually happening, and how it was found

Test cases in this project produced substantially different results
when built with GCC 9.4 versus GCC 15.2 on otherwise-matched hardware
and flags. `-march`, `-ffp-contract=off`, and the multigrid sweep
order (Gauss-Seidel vs. red/black) were all ruled out in turn; the
divergence tracked exactly with `-ftree-slp-vectorize` (SLP
"superword-level parallelism" auto-vectorization), enabled by
default at `-O2` in both compiler versions but evidently exercising a
different, version-specific cost model. Disabling it alone
(`-fno-tree-slp-vectorize`) fixed the discrepancy.

SLP vectorization does not require a loop: it packs any two (or more)
*isomorphic* statements sitting in the same basic block -- same
shape, same operators, differing only in which operand/field they
touch -- into a single vector instruction, whether or not those
statements originated from a `foreach_dimension()` expansion. Two
consequences made this dangerous here rather than merely a rounding
curiosity:

1. Several such statement pairs feed a `sqrt`/`pow` call (e.g. the
   `a`/`b` transverse-slope weights in `height_closest`, 2D and 3D).
   A compiler is free to satisfy these via a vectorized math-library
   routine instead of two scalar `libm` calls, and vector math
   implementations are not guaranteed bit-identical to their scalar
   counterparts -- so packing can change the *value*, not just the
   arithmetic order, of intermediate quantities.
2. Every one of these statement pairs ultimately contributes to the
   returned point `p`, which becomes `dist[]`, which is compared
   against a hard threshold
   (`fabs(dist[]) < 2.*L0/((double)(1 << grid->maxdepth))`) in
   `viscosity_st.h` to decide whether a cell receives surface-tension
   coefficients at all. A last-bit difference here does not stay a
   last-bit difference: it can flip that branch at a handful of
   cells, handing the multigrid solve a structurally different linear
   system rather than the same system solved to marginally different
   precision -- which is why the resulting discrepancy was large
   rather than small.

Bisection (function-level, via
`__attribute__((optimize("no-tree-slp-vectorize")))` on candidate
functions, confirmed by inspecting the qcc-generated intermediate
source to ensure the attribute actually survived macro expansion)
localised the effect specifically to `height_closest()` -- not
`height_closest_x`/`_y`/`_z`, not `residual_viscosity_st`, not
`relax_viscosity_st` (which turned out not to be SLP-vectorized by
either compiler at all), and not `redistance.h`. Within
`height_closest()`, every isomorphic statement pair that feeds the
returned `p` -- the transverse-slope weight blocks, the
normal-normalisation sum, and the final weighted-average / fallback
assignments, in both the 2D and 3D branches -- has been rewritten
below using `slp_barrier()` to fix its evaluation to a specific
scalar order, which is functionally equivalent to
`-fno-tree-slp-vectorize` but scoped only to the handful of
expressions that actually matter, leaving the rest of the file free
for the compiler to vectorize normally. With these fixes in place,
results are compiler-independent under default flags (no
`-fno-tree-slp-vectorize` required). If retested against a different
compiler/version and a divergence reappears, treat any adjacent pair
of similarly-shaped statements -- reduction or plain assignment,
looped or hand-written -- as a suspect, not just explicit
accumulation loops. */

#ifndef SLP_BARRIER_DEFINED
#define SLP_BARRIER_DEFINED
static inline double slp_barrier (double x)
{
  volatile double tmp = x;
  return tmp;
}
#endif

double minimum (double a, double b, double c, double x)
{
  double f;
  int nmax = 10;
  do {
    f = (x + (2.*a*x + b)*(x*(a*x + b) + c))/
      (a*(6.*x*(a*x + b) + 2.*c) + b*b + 1.);
    x -= f;
  }
  while (nmax-- && fabs (f) > 1e-6);
  return x;
}

static inline bool interfacial2 (Point point, scalar c)
{
  if (c[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] <= 0.)
	  return true;
  }
  else if (c[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (c[i] >= 1.)
	  return true;
  }
  else // c[] > 0. && c[] < 1.
    return true;
  return false;
}

#if dimension == 2

#define XMIN 1.5

foreach_dimension()
coord height_closest_y (Point point, vector h)
{
  coord p = { nodata, nodata };
  if (h.y[-1] != nodata && h.y[] != nodata && h.y[1] != nodata) {
    double a = (h.y[1] + h.y[-1] - 2.*h.y[])/2.;
    double b = (h.y[1] - h.y[-1])/2.;
    double c = height (h.y[]);
    double xm = nodata, ym = nodata, dm = HUGE;
    for (double x = -1; x <= 1.; x += 1.) {
      double xp = minimum (a, b, c, x), yp = xp*(a*xp + b) + c;
      double dp = sq(xp) + sq(yp);
      if (dp < dm)
	xm = xp, ym = yp, dm = dp;
    }
    if (fabs (xm) < XMIN)
      p.x = xm*Delta, p.y = ym*Delta;
  }
  return p;
}

coord height_closest (Point point, vector h, scalar c, int * s)
{
  coord p;
  coord qx = height_closest_x (point, h);
  coord qy = height_closest_y (point, h);

  if (qx.x != nodata) {
    *s = (orientation(h.x[]) ? height(h.x[]) < 0. : height(h.x[]) > 0.) ? 1 : -1;
    if (qy.x != nodata) {
      /* Slope of each height function in its transverse direction.
         nodata neighbor → 1e11 (interface out of stencil ≡ very steep slope).
         The three pairs of statements below (a/b computed, normalised,
         then weighted) are structurally identical, differing only in
         which field (h.x vs h.y) they touch -- exactly the isomorphic-
         statement shape SLP auto-vectorization targets even without a
         loop. If the compiler packs a/b into a vectorized sqrt/pow
         call pair, it can use a different (less precise) vector math
         implementation than the scalar libm call, changing the actual
         VALUES of a and b, not just their combination order -- so
         barrier-fixing only the final weighted average (as done
         previously) was not sufficient. Each stage is barrier-wrapped
         here to force scalar evaluation throughout. */
      double a = slp_barrier ((h.x[0,1] != nodata && h.x[0,-1] != nodata) ?
                 fabs (h.x[0,1] - h.x[0,-1])/2. : 1e6);
      double b = slp_barrier ((h.y[1]   != nodata && h.y[-1]   != nodata) ?
                 fabs (h.y[1]   - h.y[-1])  /2. : 1e6);
      a = slp_barrier (a / slp_barrier (sqrt (1. + a*a)));    /* sin(theta_x) = |ny|     */
      b = slp_barrier (b / slp_barrier (sqrt (1. + b*b)));    /* sin(theta_y) = |nx|     */
      a = slp_barrier (pow (slp_barrier (1. - a), 4));     /* weight for QX           */
      b = slp_barrier (pow (slp_barrier (1. - b), 4));     /* weight for QY           */
      foreach_dimension()
        p.x = slp_barrier (a*qx.x + b*qy.x)/(a + b);   /* weighted average,
							   order-fixed */
    }
    else
      p = qx;
  }
  else {
    if (qy.x != nodata){
      *s = (orientation(h.y[]) ? height(h.y[]) < 0. : height(h.y[]) > 0.) ? 1 : -1;
    p = qy;}
    else{
     if (interfacial2(point, c)){
	coord nn = interface_normal (point, c);
	double alpha = plane_alpha (c[], nn);
        double nn2 = 0.;
    	foreach_dimension() nn2 = slp_barrier (nn2 + slp_barrier (sq(nn.x)));
	foreach_dimension()
		p.x = slp_barrier (alpha / nn2 * nn.x * Delta);
	*s = (c[] > 0.5) ? 1 : -1;

     }
    }
 
  }
  return p;
}


trace
void height_distance (scalar c, scalar d, double weight = 1., int imax = 16)
{
  vector h[];
  heights (c, h);
  
  foreach() {
    int s;
    coord p = height_closest (point, h, c, &s);
    if (p.x != nodata){
      d[] = (1. - weight)*d[] + weight*s*sqrt(sq(p.x) + sq(p.y));
      bool interfacial = false;
      /*foreach_neighbor(1){
       if (interfacial2 (point, c))
        interfacial = true;
      }
      if (!interfacial)
        d[] = (c[] > 0.5) ? 10*Delta : -10*Delta;*/
    }
    else{
      bool isrefined = false;
      foreach_neighbor(1)
	if (is_refined(cell))
	   isrefined = true;
      if (c[] > 0.5){
       d[] = isrefined ? 3*Delta : 6*Delta;
      }
      else{
         d[] = isrefined ? -3*Delta : -6*Delta;
      }
    }
  }

  redistance (d, imax = imax, cfl = 0.5, phixxmin = HUGE);

}

#undef XMIN

#endif

#if dimension == 3

/* -----------------------------------------------------------------------
   Analytical closest point on the tangent plane (strictly convex,
   always gives a good starting guess near the interface).
   Solves: min_{t1,t2}  t1^2 + t2^2 + (F + D*t1 + E*t2)^2
   ----------------------------------------------------------------------- */
static void tangent_plane_closest (double D, double E, double F,
                                   double * t1_0, double * t2_0)
{
  double denom = 1. + D*D + E*E;   /* always >= 1 */
  *t1_0 = -D*F / denom;
  *t2_0 = -E*F / denom;
}

/* -----------------------------------------------------------------------
   Levenberg-Marquardt Newton for closest point on paraboloid:
   col(t1,t2) = A*t1^2 + B*t2^2 + C*t1*t2 + D*t1 + E*t2 + F
   Minimises D^2 = t1^2 + t2^2 + col^2.
   LM regularisation guarantees a descent direction even when D^2
   is locally non-convex (cell inside osculating sphere).
   Step tolerance matches the 1D minimum() solver exactly.
   ----------------------------------------------------------------------- */
static void minimum3d (double A, double B, double C,
                       double D, double E, double F,
                       double * t1p, double * t2p)
{
  double t1 = *t1p, t2 = *t2p;
  const double eps  = 1e-6;   /* minimum eigenvalue margin    */
  const double stol = 1e-5;   /* step tolerance (matches 1D)  */
  double mu = 0.;
  int nmax = 40;

  do {
    double col  = A*t1*t1 + B*t2*t2 + C*t1*t2 + D*t1 + E*t2 + F;
    double ct1  = 2.*A*t1 + C*t2 + D;   /* d(col)/d(t1) */
    double ct2  = C*t1 + 2.*B*t2 + E;   /* d(col)/d(t2) */

    /* Gradient of D^2/2 */
    double ft1  = t1 + col*ct1;
    double ft2  = t2 + col*ct2;

    /* Hessian of D^2/2 */
    double Hxx  = 1. + ct1*ct1 + 2.*A*col;
    double Hxy  = ct1*ct2 + C*col;
    double Hyy  = 1. + ct2*ct2 + 2.*B*col;

    /* Minimum eigenvalue (analytical 2x2) */
    double tr   = Hxx + Hyy;
    double diff = Hxx - Hyy;
    double disc = sqrt (diff*diff + 4.*Hxy*Hxy);
    double lmin = 0.5*(tr - disc);
    mu = (lmin < eps) ? (eps - lmin) : 0.;

    /* LM-regularised Newton step */
    double rHxx = Hxx + mu;
    double rHyy = Hyy + mu;
    double det  = rHxx*rHyy - Hxy*Hxy;  /* > 0 guaranteed by LM */
    if (fabs (det) < 1e-14) break;
    double dt1  = -( rHyy*ft1 - Hxy*ft2) / det;
    double dt2  = -(-Hxy*ft1 + rHxx*ft2) / det;
    t1 += dt1;
    t2 += dt2;

    if (fabs (dt1) < stol && fabs (dt2) < stol) break;
  }
  while (nmax--);

  *t1p = t1;  *t2p = t2;
}

/* -----------------------------------------------------------------------
   Find the closest point on the interface to the cell centre, using
   the height function h.z (column direction = z, transverse = x,y).
   foreach_dimension() generates height_closest_x and height_closest_y.

   Variable naming convention (avoids shadowing Basilisk coordinates):
     t1  = first transverse direction  (x in height_closest_z template)
     t2  = second transverse direction (y in height_closest_z template)
     tc  = column direction            (z in height_closest_z template)
   ----------------------------------------------------------------------- */
#define XMIN 1.5

foreach_dimension()
coord height_closest_z (Point point, vector h, int * hierarchy)
{
  *hierarchy = 4;
  coord p = {nodata, nodata, nodata};
  if (h.z[0,0,0] == nodata) return p;
  double F = height (h.z[0,0,0]);

  /* Transverse stencil availability                                    */
  bool ht1m = (h.z[-1, 0, 0] != nodata), ht1p = (h.z[ 1, 0, 0] != nodata);
  bool ht2m = (h.z[ 0,-1, 0] != nodata), ht2p = (h.z[ 0, 1, 0] != nodata);
  bool have_cross = ht1m && ht1p && ht2m && ht2p;
  bool have_all_corners = have_cross &&
                          (h.z[-1,-1, 0] != nodata) &&
                          (h.z[ 1,-1, 0] != nodata) &&
                          (h.z[-1, 1, 0] != nodata) &&
                          (h.z[ 1, 1, 0] != nodata);

  /* ------------------------------------------------------------------ */
  /* Level 1/2: full 2D paraboloid fit                                   */
  /* Level 1 = all corners (includes twist C), O(Delta^2)               */
  /* Level 2 = cross only, C=0,               O(Delta^2)               */
  /* ------------------------------------------------------------------ */
  if (have_cross) {
    double A  = (h.z[ 1, 0, 0] + h.z[-1, 0, 0] - 2.*h.z[0,0,0]) / 2.;
    double B  = (h.z[ 0, 1, 0] + h.z[ 0,-1, 0] - 2.*h.z[0,0,0]) / 2.;
    double D  = (h.z[ 1, 0, 0] - h.z[-1, 0, 0]) / 2.;
    double E  = (h.z[ 0, 1, 0] - h.z[ 0,-1, 0]) / 2.;
    double C  = have_all_corners ?
                (h.z[ 1, 1, 0] - h.z[-1, 1, 0]
               - h.z[ 1,-1, 0] + h.z[-1,-1, 0]) / 4. : 0.;

    /* Tangent-plane starting guess (always near interface, D^2 convex) */
    double t1_0, t2_0;
    tangent_plane_closest (D, E, F, &t1_0, &t2_0);

    /* Search all seeds, keep global minimum                             */
    struct { double t1, t2, tc, d; } best = {nodata, nodata, nodata, HUGE};

    double s1[] = {t1_0,-1., 0., 1.,-1., 0., 1.,-1., 0., 1.};
    double s2[] = {t2_0,-1.,-1.,-1., 0., 0., 0., 1., 1., 1.};

    for (int k = 0; k < 10; k++) {
      double t1 = s1[k], t2 = s2[k];
      minimum3d (A, B, C, D, E, F, &t1, &t2);
      if (fabs (t1) < XMIN && fabs (t2) < XMIN) {
        double tc = A*t1*t1 + B*t2*t2 + C*t1*t2 + D*t1 + E*t2 + F;
        double dd = sq(t1) + sq(t2) + sq(tc);
        if (dd < best.d)
          best.t1 = t1, best.t2 = t2, best.tc = tc, best.d = dd;
      }
    }

    if (best.t1 != nodata) {
      p.x = best.t1 * Delta;   /* t1 = x in z-template */
      p.y = best.t2 * Delta;   /* t2 = y in z-template */
      p.z = best.tc * Delta;   /* col= z in z-template */
    }
    *hierarchy = have_all_corners ? 1 : 2;
    return p;
  }

/* ------------------------------------------------------------------ */
  /* Level 3: single-row fallback.                                       */
  /* Compute both available 1D parabolas and pick the more informative.  */
  /* A row is "informative" if it has non-negligible slope or curvature. */
  /* If the only available row is flat (a~0, b~0 -- i.e. the column     */
  /* direction is parallel to the interface), return nodata: the cell    */
  /* is too far from the interface for the height function to give a     */
  /* useful transverse estimate. Redistancing handles these cells.       */
  /* ------------------------------------------------------------------ */
  {
    /* t1-row coefficients (if available) */
    double a1 = 0., b1 = 0., info1 = 0.;
    if (ht1m && ht1p) {
      a1 = (h.z[ 1,0,0] + h.z[-1,0,0] - 2.*h.z[0,0,0]) / 2.;
      b1 = (h.z[ 1,0,0] - h.z[-1,0,0]) / 2.;
      info1 = sq(a1) + sq(b1);
    }

    /* t2-row coefficients (if available) */
    double a2 = 0., b2 = 0., info2 = 0.;
    if (ht2m && ht2p) {
      a2 = (h.z[0, 1,0] + h.z[0,-1,0] - 2.*h.z[0,0,0]) / 2.;
      b2 = (h.z[0, 1,0] - h.z[0,-1,0]) / 2.;
      info2 = sq(a2) + sq(b2);
    }

    /* Neither row has useful information: return nodata */
    if (info1 < 1e-10 && info2 < 1e-10) return p;

    /* Pick the more informative row */
    double a = (info1 >= info2) ? a1 : a2;
    double b = (info1 >= info2) ? b1 : b2;

    double t0  = -b*F / (1. + b*b);
    double tm  = nodata, dm = HUGE;
    double st[] = {t0, -1., 0., 1.};
    for (int k = 0; k < 4; k++) {
      double t1 = minimum (a, b, F, st[k]);
      if (fabs (t1) < XMIN) {
        double tc = a*t1*t1 + b*t1 + F;
        double dd = sq(t1) + sq(tc);
        if (dd < dm) tm = t1, dm = dd;
      }
    }

    if (tm != nodata) {
      if (info1 >= info2) {
        p.x = tm * Delta;   /* t1 direction */
        p.y = 0.;
        p.z = (a*tm*tm + b*tm + F) * Delta;
      } else {
        p.x = 0.;
        p.y = tm * Delta;   /* t2 direction */
        p.z = (a*tm*tm + b*tm + F) * Delta;
      }
    }
    /* If tm==nodata: closest point outside XMIN, return nodata.
       Redistancing will propagate from nearby well-covered cells.     */
    *hierarchy = 3;
    return p;
  }

  return p;  /* Level 4: nodata */
}

#undef XMIN

/* -----------------------------------------------------------------------
   Combine estimates from all three height-function directions.
   Sign is taken from whichever direction has a valid height function.
   ----------------------------------------------------------------------- */
/* -----------------------------------------------------------------------
   Combine estimates from all three height-function directions.
   Sign is taken from the VOF field c[] directly (c[]>0.5 = inside = negative SDF).
   The old height-function-based sign determination is left as a comment below;
   it fails for closed objects near axis-aligned poles where the column stencil
   captures only one interface crossing and the column sum exceeds 3.5 despite
   the cell being geometrically inside the object.
   ----------------------------------------------------------------------- */
coord height_closest (Point point, vector h, scalar c, int *s)
{
  coord p;
  int hierarchy_x = 0;
  int hierarchy_y = 0;
  int hierarchy_z = 0;
  coord qx = height_closest_x (point, h, &hierarchy_x);
  coord qy = height_closest_y (point, h, &hierarchy_y);
  coord qz = height_closest_z (point, h, &hierarchy_z);

#ifdef DIAGNOSE_HD
  if (qx.x != nodata)
    fprintf (ferr, "QX (%g %g %g): px=%g py=%g pz=%g dist=%g hierarchy=%d\n",
             x,y,z, qx.x,qx.y,qx.z, sqrt(sq(qx.x)+sq(qx.y)+sq(qx.z)), hierarchy_x);
  if (qy.x != nodata)
    fprintf (ferr, "QY (%g %g %g): px=%g py=%g pz=%g dist=%g hierarchy=%d\n",
             x,y,z, qy.x,qy.y,qy.z, sqrt(sq(qy.x)+sq(qy.y)+sq(qy.z)), hierarchy_y);
  if (qz.x != nodata)
    fprintf (ferr, "QZ (%g %g %g): px=%g py=%g pz=%g dist=%g hierarchy=%d\n",
             x,y,z, qz.x,qz.y,qz.z, sqrt(sq(qz.x)+sq(qz.y)+sq(qz.z)), hierarchy_z);
#endif

  /* Weights from h-function STENCIL SLOPES, not from the returned p-vector.
     Missing neighbor → 1e11 (nodata means interface moves sharply out of
     stencil range, implying the column is nearly tangent to the interface).

     Each block below computes two transverse slope terms via
     structurally identical ternary expressions (e.g. sy/sz), then
     combines them into a single weight. This is the same isomorphic-
     statement-pair shape that turned out to be the actual source of
     compiler-dependent divergence in the 2D height_closest() (the
     a/b slope block there) -- SLP auto-vectorization can pack two
     adjacent, similarly-shaped statements even with no loop present,
     potentially using a different (less precise) vectorized sqrt/pow
     implementation than the scalar path. Barrier-wrapped throughout
     for the same reason. */
  double wx = 0., wy = 0., wz = 0.;

  if (qx.x != nodata) {
    /* QX: col=x, t1=y, t2=z */
    double sy = slp_barrier ((h.x[0, 1,0] != nodata && h.x[0,-1,0] != nodata) ?
                fabs (h.x[0,1,0] - h.x[0,-1,0]) / 2. : 1e11);
    double sz = slp_barrier ((h.x[0,0, 1] != nodata && h.x[0,0,-1] != nodata) ?
                fabs (h.x[0,0,1] - h.x[0,0,-1]) / 2. : 1e11);
    double slope = slp_barrier (sqrt (slp_barrier (sq(sy) + sq(sz))));
    wx = slp_barrier (pow (slp_barrier (1. - slope/slp_barrier (sqrt (1. + sq(slope)))), 4));
    wx = max (wx, 1e-8);
  }

  if (qy.x != nodata) {
    /* QY: col=y, t1=z, t2=x */
    double sz = slp_barrier ((h.y[0,0, 1] != nodata && h.y[0,0,-1] != nodata) ?
                fabs (h.y[0,0,1] - h.y[0,0,-1]) / 2. : 1e11);
    double sx = slp_barrier ((h.y[ 1,0,0] != nodata && h.y[-1,0,0] != nodata) ?
                fabs (h.y[1,0,0] - h.y[-1,0,0]) / 2. : 1e11);
    double slope = slp_barrier (sqrt (slp_barrier (sq(sz) + sq(sx))));
    wy = slp_barrier (pow (slp_barrier (1. - slope/slp_barrier (sqrt (1. + sq(slope)))), 4));
    wy = max (wy, 1e-8);
  }

  if (qz.x != nodata) {
    /* QZ: col=z, t1=x, t2=y */
    double sx = slp_barrier ((h.z[ 1,0,0] != nodata && h.z[-1,0,0] != nodata) ?
                fabs (h.z[1,0,0] - h.z[-1,0,0]) / 2. : 1e11);
    double sy = slp_barrier ((h.z[0, 1,0] != nodata && h.z[0,-1,0] != nodata) ?
                fabs (h.z[0,1,0] - h.z[0,-1,0]) / 2. : 1e11);
    double slope = slp_barrier (sqrt (slp_barrier (sq(sx) + sq(sy))));
    wz = slp_barrier (pow (slp_barrier (1. - slope/slp_barrier (sqrt (1. + sq(slope)))), 4));
    wz = max (wz, 1e-8);
  }

  /* wtot and the weighted blend of qx/qy/qz below feed directly into
     the returned closest-point vector p, which becomes dist[] in
     height_distance() -- and dist[] gates the interfacial-band
     threshold comparison in viscosity_st.h. The three accumulation
     statements below (one per direction, each touching all of
     p.x/p.y/p.z via foreach_dimension()) are structurally identical
     and sit in the same basic block, which is exactly the pattern
     SLP auto-vectorization targets for reassociation. Each addition
     is barrier-fixed to a left-to-right order so the result is the
     same regardless of compiler/SLP heuristics. */
  double wtot = slp_barrier (slp_barrier (wx + wy) + wz);
  int n = 0;
  foreach_dimension() p.x = 0.;

  if (qx.x != nodata) {
    foreach_dimension() p.x = slp_barrier (p.x + slp_barrier (wx * qx.x));
    n++;
  }
  if (qy.x != nodata) {
    foreach_dimension() p.x = slp_barrier (p.x + slp_barrier (wy * qy.x));
    n++;
  }
  if (qz.x != nodata) {
    foreach_dimension() p.x = slp_barrier (p.x + slp_barrier (wz * qz.x));
    n++;
  }

  if (n > 0) {
    foreach_dimension() p.x /= wtot;

    /* Sign from VOF field: c[]>0.5 means inside the reference phase → negative SDF.
       This is unambiguous for any interface geometry.
       
       Old height-function-based sign (kept for reference):
       Sign was taken from the first available height function direction,
       using the column sum relative to 3.5 to determine inside/outside:
         if (qx.x != nodata)
           *s = (orientation(h.x[]) ? height(h.x[]) < 0.
                                     : height(h.x[]) > 0.) ? 1 : -1;
         else if (qy.x != nodata)
           *s = (orientation(h.y[]) ? height(h.y[]) < 0.
                                     : height(h.y[]) > 0.) ? 1 : -1;
         else
           *s = (orientation(h.z[]) ? height(h.z[]) < 0.
                                     : height(h.z[]) > 0.) ? 1 : -1;
       This fails for closed objects near axis-aligned poles: when the 7-cell
       column stencil captures only the near-side interface crossing of a sphere,
       the column sum exceeds 3.5 for inside cells, inverting the sign.
       Example: cell at (-0.684, -0.027, -0.020) inside sphere R=0.7, normal≈(-1,0,0):
         column sum = 5.50 cells → height > 0 → s=+1 (outside) → WRONG.             */
    *s = (c[] > 0.5) ? 1 : -1;
  }
  else {
    p.x = nodata;
    /* Sign still needed for fallback in height_distance */
    *s = (c[] > 0.5) ? 1 : -1;
    if (interfacial2(point, c)){
        coord nn = interface_normal (point, c);
        double alpha = plane_alpha (c[], nn);
        double nn2 = 0.;
        foreach_dimension() nn2 = slp_barrier (nn2 + slp_barrier (sq(nn.x)));
        foreach_dimension()
                p.x = slp_barrier (alpha / nn2 * nn.x * Delta);
     }
  }

  return p;
}


/* -----------------------------------------------------------------------
   Main entry point: compute signed distance field from VOF field c.
   ----------------------------------------------------------------------- */
trace
void height_distance (scalar c, scalar d, double weight = 1., int imax = 16)
{
  vector h[];
  heights (c, h);

  foreach() {
    int s;
    //coord p = height_closest (point, h, &s);
coord p = height_closest (point, h, c, &s);
    if (p.x != nodata){
      d[] = (1. - weight)*d[] +
            weight*s*sqrt (sq(p.x) + sq(p.y) + sq(p.z));
      /*bool interfacial = false;
      foreach_neighbor(1){
       if (interfacial2 (point, c))
    	interfacial = true;
      }
	foreach_neighbor(2){
	    if (c[] > 0 && c[] < 1)
        interfacial = true;
      }
	for ( int kk = -2 ; kk <= 2; kk++){
	for ( int ii = -2 ; ii <= 2; ii++){
		for (int jj = -2; jj <= 2; jj++){
			if (c[kk,ii,jj] > 0 && c[kk,ii,jj] < 1 && (fabs(ii) + fabs(jj) + fabs(kk)) <= 6)
 				interfacial = true;


		}
	}
}
      if (!interfacial)
      	d[] = (c[] > 0.5) ? 10*Delta : -10*Delta;*/
      }
      else{
      bool isrefined = false;
      foreach_neighbor(1)
        if (is_refined(cell))
           isrefined = true;
      if (c[] > 0.5){
       d[] = isrefined ? 3*Delta : 6*Delta;
      }
      else{
         d[] = isrefined ? -3*Delta : -6*Delta;
      }
    }
  }

#ifdef DIAGNOSE_HD
  foreach()
    if (fabs(d[]) != 6*Delta){
 	    fprintf (ferr, "BEFORE_REDISTANCE (%g %g %g): dist=%g\n",
             x,y,z,d[]);
  }
#endif

  redistance (d, imax = imax, cfl = 0.5, phixxmin = HUGE);
#ifdef DIAGNOSE_HD
  foreach()
    if (fabs(d[]) < 2*Delta){
            fprintf (ferr, "AFTER_REDISTANCE (%g %g %g): dist=%g\n",
             x,y,z,d[]);
  }
#endif

}

#endif /* dimension == 3 */