/**
# Common Evaporation

We define tolerances and functions useful
for evaporation models. */

#ifndef F_ERR
# define F_ERR 1.e-10
#endif
#ifndef P_ERR
# define P_ERR 1.e-10
#endif
#ifndef T_ERR
# define T_ERR 0.
#endif
#ifndef I_TOL
# define I_TOL 1.e-8
#endif

//#include "curvature.h"
#include "adapt_wavelet_leave_interface.h"
#include "mass_refine_prolongation.h"

macro foreach_interfacial (scalar f, double tol = 1e-10,
    Reduce reductions = None)
{
  foreach(0, reductions) {
    POINT_VARIABLES();
    if (f[] > tol && f[] < 1. - tol)
    {...}
  }
}

macro foreach_interfacial_plic (scalar f, double tol = 1e-10,
    Reduce reductions = None)
{
  foreach_interfacial (f, tol, reductions) {
    coord m = interface_normal (point, f);
    double alpha = plane_alpha (f[], m);
    coord prel;
    double area = plane_area_center (m, alpha, &prel);
#if AXI
    double dirac = area*(y + prel.y*Delta)/(Delta*y)*cm[];
#else
    double dirac = area/Delta*cm[];
#endif
    NOT_UNUSED (dirac);
    {...}
  }
}

/**
We define a macro describing a radial
profile. It can be useful for field initialization. */

macro number radialprofile (double r,
    double r1, double r2,
    double T1, double T2)
{
  return (r <= r1) ? T1
       : (r >= r2) ? T2
       : (-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2));
}

/**
## Compute the normal in an interfacial cell
*/

coord normal (Point point, scalar c)
{
  coord n = interface_normal (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
## *avg_neighbor()*: Compute the average value of a scalar field Y in a 3x3 stencil around the current cell
* *point*: current cell location
* *Y*: field to average
* *f*: vof volume fraction field
*/

double avg_neighbor (Point point, scalar Y, scalar f) {
  double fYnei = 0., fnei = 0.;
  foreach_neighbor(1) {
    double ff = Y.inverse ? 1. - f[] : f[];
    fYnei += ff*Y[];
    fnei  += ff;
  }
  return fYnei / fnei;
}

/**
## *avg_interface()*: Compute the average value of a scalar field Y in a region around the interface
* *Y*: fields to average
* *f*: vof volume fraction field
*/

double avg_interface (scalar Y, scalar f, double tol=F_ERR) {
  double Yavg = 0.;
  int counter = 0;
  foreach(reduction(+:Yavg) reduction(+:counter)) {
    if (f[] > tol && f[] < 1.-tol) {
      counter++;
      Yavg += Y[];
    }
  }
  return Yavg / counter;
}

/**
## *vof_reconstruction()*: VOF reconstruction step

The VOF reconstruction step must be frequently performed to compute
the Dirac delta which allows surface integrals to be transformed
into volume intergrals. Used for the evaporation source terms.

* *point*: current cell in a *foreach()* loop
* *f*: vof field
*/

typedef struct {
  coord m, prel;
  double alpha, area, dirac;
} vofrecon;

vofrecon vof_reconstruction (Point point, scalar f) {
  coord m = interface_normal (point, f);
  double alpha = plane_alpha (f[], m);
  coord prel;
  double area = plane_area_center (m, alpha, &prel);
#if AXI
  double dirac = area*(y + prel.y*Delta)/(Delta*y)*cm[];
#else
  double dirac = area/Delta*cm[];
#endif

  return (vofrecon){m, prel, alpha, area, dirac};
}

/**
## *shift_field()*: Shift a field localized at the interface toward the closest pure gas or liquid cells

* *fts*: field to shift
* *f*: vof volume fraction field
* *dir*: shifting direction: 1 liquid, gas otherwise
*/

scalar avg[];

trace
void shift_field (scalar fts, scalar f, int dir) {

  // scalar avg[]; // fixme: if avg is local and not global, there is a problem.
  // My impression is that `no_restriction` is not working - but it should be
  // investigated more deeply
  avg.c = f;
#if TREE
  avg.refine = refinement_avg;
  set_prolongation (avg, refinement_avg);
  set_restriction (avg, no_restriction);
#endif

  // Compute avg
  foreach() {
    avg[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      if (dir == 1) {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] > 1.-F_ERR && cm[]) // Number of pure-liquid cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
      else {
        int count = 0;
        foreach_neighbor (1) {
          if (f[] < F_ERR && cm[]) // Number of pure-gas cells close to the interfacial cell
            count ++;
        }
        avg[] = count;
      }
    }
  }

  scalar sf0[];
  foreach() {
    sf0[] = fts[];
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      fts[] = 0.;
  }

  // Compute m
  foreach() {
    if (dir == 1) {
      if (f[] > 1.-F_ERR && cm[]) { // Move toward pure-liquid
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] += val;
      }
    }
    else {
      if (f[] < F_ERR && cm[]) { // Move toward pure-gas
        double val = 0.;
        foreach_neighbor (1) {
          if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
            val += sf0[]/avg[];
          }
        }
        fts[] += val;
      }
    }
  }
}

#include "diffusion-flux.h"

trace
void shift_diffusion (scalar fts, scalar f,
    double cells = 0.1, double tol = 1e-3)
{
  foreach()
    fts[] = (cm[] > 0) ? fts[]/cm[] : 0.;
  double totold = statsf (fts).sum;

  double delta = L0/(1 << grid->maxdepth);
  double diff = cells*sq(delta);

  scalar theta[];
  foreach()
    theta[] = cm[];

  face vector D[];
  foreach_face() {
    D.x[] = diff*fm.x[];
    if (f[-1] > F_ERR && f[-1] < 1.-F_ERR && fts[] != 0.)
      D.x[] *= 0.;
    else if (f[] > F_ERR && f[] < 1.-F_ERR && fts[-1] != 0.)
      D.x[] *= 0.;
  }

  double TOLERANCE_BACKUP = TOLERANCE;
  TOLERANCE = tol;
  diffusion (fts, 1., D = D, theta = theta);
  TOLERANCE = TOLERANCE_BACKUP;

  double totnew = statsf (fts).sum;
  foreach()
    fts[] *= cm[]*totold/totnew;
}

/**
## *copy_bcs()*: Copy the boundary conditions from a target field to a list of fields

* *fts*: field to shift
* *f*: vof volume fraction field
* *dir*: shifting direction: 1 liquid, gas otherwise
*/

void copy_bcs (scalar * dest, scalar orig) {
  for (int bid = 0; bid < nboundary; bid++) {
    for (scalar s in dest) {
      s.boundary[bid] = orig.boundary[bid];
      s.boundary_homogeneous[bid] = orig.boundary_homogeneous[bid];
    }
  }
}

#if TREE
attribute {
  scalar rho;
}

void density_refine (Point point, scalar rhov) {
  refine_bilinear (point, rhov);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhov[];
  double drho = rhov[] - rhou/((1 << dimension)*(cm[] + SEPS));
  foreach_child()
    rhov[] += drho;
}

void density_restriction (Point point, scalar rhov) {
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhov[];
  rhov[] = rhou/((1 << dimension)*(cm[] + SEPS));
  //restriction_volume_average (point, rhov);
}

void restriction_mass_average (Point point, scalar s) {
  scalar rhov = s.rho;
  double sum = 0., mass = 0.;
  foreach_child() {
    sum += cm[]*rhov[]*s[];
    mass += cm[]*rhov[];
  }
  s[] = sum/(mass + 1e-30);
}
#endif


/**
## Utilities for evaporation with embedded boundaries

When the interface between the liquid and the gas phase meets a solid described
by an embedded boundary (see [embed.h](/src/embed.h)), the fragment of embedded
boundary contained in the cut cells at the contact line is wetted partly by the
liquid and partly by the gas. The following functions compute the portion of
the fragment wetted by each phase, which is used to split the fluxes through the
embedded boundary (e.g. a heat flux imposed on the wall) between the phases. They
work in two and three dimensions. 

todo: the following functions with geometric operations were written using
Claude and never reviewed
*/

#if EMBED

/**
The scalar fields of a phase (e.g. temperature and mass fractions) can be
associated with a field which corrects the flux through the embedded boundary
computed by `embed_flux()` (see *embed_fraction()* and
`phase_embed_flux()` in [phase.h](phase.h)). As for the other scalar attributes,
the default index 0 means that the field is not set. */

attribute {
  scalar wetting;
}

/**
### *embed_polygon_area_center()*: area and centroid of a convex polygon

The vertices *v* of the polygon are ordered, and the area is the sum of the areas
of the triangles which share the first vertex. */

#if dimension == 3
static double embed_polygon_area_center (const coord * v, int nv, coord * p)
{
  double area = 0.;
  coord c = {0., 0., 0.};
  for (int i = 1; i < nv - 1; i++) {
    coord a, b;
    foreach_dimension() {
      a.x = v[i].x - v[0].x;
      b.x = v[i+1].x - v[0].x;
    }
    coord n = {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
    double da = sqrt (sq(n.x) + sq(n.y) + sq(n.z))/2.;
    foreach_dimension()
      c.x += da*(v[0].x + v[i].x + v[i+1].x)/3.;
    area += da;
  }
  if (area > 0.)
    foreach_dimension()
      c.x /= area;
  *p = c;
  return area;
}
#endif

/**
### *embed_wetted_fraction()*: fraction of the embedded fragment wetted by a phase

This function returns the fraction of the fragment of embedded boundary
contained in the cell which lies in the phase with volume fraction *c*. The
fragment, reconstructed from `cs` and `fs`, is intersected with the interface
reconstructed from *c*. Consistently with [vof.h](/src/vof.h), which ignores the
volume occupied by the solid in the cut cells, the interface is reconstructed as
if the whole cell was occupied by the fluid. The phase lies on the side of the
interface where $\mathbf{m}\cdot\mathbf{x} < \alpha$, with $\mathbf{m}$ the
normal to the interface and $\alpha$ its intercept, therefore the wetted
fractions of the phases with volume fractions $c$ and $1 - c$ sum to one. If a
metric is defined (`metric_embed_factor()`, e.g. in axisymmetric domains), the
fraction refers to the area of the boundary including the metric factor. Cells
which do not contain a fragment of embedded boundary return zero.

* *c*: volume fraction of the phase
* *tol*: tolerance on the volume fraction
*/

double embed_wetted_fraction (Point point, scalar c, double tol = 1e-10)
{
  if (cs[] <= 0. || cs[] >= 1.)
    return 0.;
  double cc = clamp (c[], 0., 1.);
  if (cc <= tol)
    return 0.;
  if (cc >= 1. - tol)
    return 1.;

  /**
  We reconstruct the fragment of embedded boundary and the interface, in the
  coordinates of the cell. */

  coord n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (cs[], n);
  coord m = interface_normal (point, c);
  double alphac = plane_alpha (cc, m);

#if dimension == 2

  /**
  In two dimensions, the fragment is the segment from $\mathbf{a}$ to
  $\mathbf{b}$, and the portion with parameter $s_0 \le s \le s_1$ lies in the
  phase. */

  coord p[2];
  if (facets (n, alpha, p) < 2)
    return cc;
  coord a = p[0], b = p[1];
  double da = m.x*a.x + m.y*a.y - alphac;
  double db = m.x*b.x + m.y*b.y - alphac;
  double s0 = 0., s1 = 0.;
  if (da <= 0. && db <= 0.)
    s1 = 1.;
  else if (da < 0. || db < 0.) {
    double s = da/(da - db);
    if (da < 0.)
      s1 = s;
    else
      s0 = s, s1 = 1.;
  }
  double w = s1 - s0;

  /**
  The metric factor is linear along the segment, therefore its averages over
  the wetted portion and over the whole segment are the values at their
  midpoints. */

  if (metric_embed_factor && w > 0.) {
    coord pw, pc;
    foreach_dimension() {
      pw.x = a.x + (b.x - a.x)*(s0 + s1)/2.;
      pc.x = (a.x + b.x)/2.;
    }
    double mc = metric_embed_factor (point, pc);
    if (mc > 0.)
      w *= metric_embed_factor (point, pw)/mc;
  }
  return w;

#else // dimension == 3

  /**
  In three dimensions, the fragment is a convex polygon, which is clipped by the
  plane of the interface: we keep the vertices which lie in the phase, and we
  add the intersections of the edges with the interface. The clipped polygon has
  at most one vertex more than the fragment. */

  coord v[12], vw[13];
  int nv = facets (n, alpha, v, 1.), nw = 0;
  if (nv < 3)
    return cc;
  for (int i = 0; i < nv; i++) {
    coord a = v[i], b = v[(i + 1) % nv];
    double da = m.x*a.x + m.y*a.y + m.z*a.z - alphac;
    double db = m.x*b.x + m.y*b.y + m.z*b.z - alphac;
    if (da <= 0.)
      vw[nw++] = a;
    if ((da < 0. && db > 0.) || (da > 0. && db < 0.)) {
      double s = da/(da - db);
      foreach_dimension()
        vw[nw].x = a.x + s*(b.x - a.x);
      nw++;
    }
  }
  coord pc, pw;
  double area = embed_polygon_area_center (v, nv, &pc);
  if (area <= 0.)
    return cc;
  double w = nw > 2 ? embed_polygon_area_center (vw, nw, &pw)/area : 0.;
  if (metric_embed_factor && w > 0.) {
    double mc = metric_embed_factor (point, pc);
    if (mc > 0.)
      w *= metric_embed_factor (point, pw)/mc;
  }
  return w;

#endif // dimension == 3
}

/**
### *embed_fraction()*: wetting correction of the flux through the embedded boundary

The flux through the fragment of embedded boundary computed by `embed_flux()`
(see [embed.h](/src/embed.h)) uses the coefficient
$\sum_f \mu_f/\sum_f f_{m,f}$, i.e. the average of the coefficient $\mu$ over the
faces of the cell. When the coefficient of a phase includes the face fractions
$f_c$ of the phase, $\mu = \lambda f_m f_c$ (as in [phase.h](phase.h)), the flux
of the cut cells at the contact line is split between the phases according to
the average face fraction of each phase, $\sum f_m f_c/\sum f_m$, rather than to
the portion of the fragment wetted by each phase. The error is of $O(\Delta)$,
and it depends on the position of the interface within the cells at the contact
line. This function computes the factor
$$
  w = w_{wet}\dfrac{\sum_f f_{m,f}}{\sum_f f_{m,f} f_{c,f}}
$$
which converts the flux of `embed_flux()` into the flux through the portion
$w_{wet}$ of the fragment wetted by the phase (see *embed_wetted_fraction()*).
In the cells which do not contain a fragment of embedded boundary $w = 1$. If no
face of the cut cell is open to the phase, the phase is isolated in the cell and
it does not receive any flux ($w = 0$), as with `embed_flux()`. Without the face
fractions *fc*, the factor is the wetted fraction $w_{wet}$.

* *w*: wetting correction (filled by this function)
* *c*: volume fraction of the phase
* *fc*: face fractions included in the coefficient of the flux
* *tol*: tolerance on the volume fraction
*/

void embed_fraction (scalar w, scalar c,
    (const) face vector fc = {{-1}}, double tol = 1e-10)
{
  if (fc.x.i < 0)
    fc[] = {1.,1.,1.};

  foreach (nowarning) {
    w[] = 1.;
    if (cs[] > 0. && cs[] < 1.) {
      double sfm = 0., sfc = 0.;
      foreach_dimension() {
        sfm += fm.x[] + fm.x[1];
        sfc += fm.x[]*fc.x[] + fm.x[1]*fc.x[1];
      }
      w[] = (sfc > 0.) ? embed_wetted_fraction (point, c, tol)*sfm/sfc : 0.;
    }
  }

  /**
  The correction is needed on all the levels of the multigrid solver. */

  restriction ({w});
}

#endif // EMBED
