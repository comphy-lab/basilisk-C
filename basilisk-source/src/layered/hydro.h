/**
# The hydrostatic multilayer solver for free-surface flows

The theoretical basis and main algorithms for this solver are
described in [Popinet, 2020](/Bibliography#popinet2020). Note however
that this version of the solver is more recent and may not match the
details of Popinet, 2020.

The system of $n$ layers of incompressible fluid pictured below is
modelled as the set of (hydrostatic) equations:
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
\end{aligned}
$$
with $\mathbf{u}_k$ the velocity vector, $h_k$ the layer thickness,
$g$ the acceleration of gravity and $\eta$ the free-surface
height. The non-hydrostatic pressure $\phi_k$ and vertical velocity
$w_k$ can be added using [the non-hydrostatic extension](nh.h). See
also the [general introduction](README).

![Definition of the $n$ layers.](../figures/layers.svg){ width="60%" }

## Fields and parameters

The `zb` and `eta` fields define the bathymetry and free-surface
height. The `h` and `u` fields define the layer thicknesses and velocities.

The acceleration of gravity is `G` and `dry` controls the minimum
layer thickness. The hydrostatic CFL criterion is defined by
`CFL_H`. It is set to a very large value by default, but will be tuned
either by the user or by the default solver settings (typically
depending on whether time integration is explicit or implicit).

The `gradient` pointer gives the function used for limiting.

`tracers` is a list of tracers for each layer. By default it contains
only the components of velocity (unless `linearised` is set to `true`
in which case the tracers list is empty). */

#define BGHOSTS 2
#define LAYERS 1
#include "utils.h"

scalar zb[], eta, h;
vector u;
double G = 1., dry = 1e-12, CFL_H = 1e40;
double (* gradient) (double, double, double) = minmod2;

scalar * tracers = NULL, eta_r;
bool linearised = false;

/**
## Setup

We first define refinement and restriction functions for the
free-surface height `eta` which is just the sum of all layer
thicknesses and of bathymetry. */

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child() {
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
  }
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[];
  foreach_layer()
    eta[] += h[];
}
#endif // TREE

/**
The allocation of fields for each layer is split in two parts, because
we need to make sure that the layer thicknesses and $\eta$ are
allocated first (in the case where other layer fields are added, for
example for the [non-hydrostatic extension](nh.h)). */

event defaults0 (i = 0)
{
  assert (nl > 0);
  h = new scalar[nl];
  h.gradient = gradient;
#if TREE
  h.refine = h.prolongation = refine_linear;
  h.restriction = restriction_volume_average;
  h.dirty = true;
#endif
  eta_r = eta = new scalar;
  reset ({h, zb}, 0.);

  /**
  We set the proper gradient and refinement/restriction functions. */
  
  zb.gradient = gradient;
  eta.gradient = gradient;
#if TREE
  zb.refine = zb.prolongation = refine_linear;
  zb.restriction = restriction_volume_average;
  zb.dirty = true;
  eta.prolongation = refine_linear;
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
  eta.depends = list_copy ({zb,h});
  eta.dirty = true;
#endif // TREE
}

#include "run.h"

/**
Other fields, such as $\mathbf{u}_k$ here, are added by this
event. */

event defaults (i = 0)
{

  /**
  The (velocity) CFL is limited by the unsplit advection scheme, so is
  dependent on the dimension. The (gravity wave) CFL is set to 1/2 (if
  not already set by the user). */
  
  CFL = 1./(2.*dimension);
  if (CFL_H == 1e40)
    CFL_H = 0.5 [0];
  
  u = new vector[nl];
  reset ({u}, 0.);
    
  if (!linearised)
    foreach_dimension()
      tracers = list_append (tracers, u.x);
  
  /**
  The gradient and prolongation/restriction functions are set for all
  tracer fields. */

  for (scalar s in tracers) {
    s.gradient = gradient;
#if TREE
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
    s.dirty = true;
#endif
  }

  /**
  We setup the default display. */

  display ("squares (color = 'eta > zb ? eta : nodata', spread = -1);");
}

/**
After user initialisation, we define the free-surface height $\eta$. */

event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
    dimensional (h[] == Delta);
  }
}

/**
The maximum timestep `dtmax` can be used to impose additional stability
conditions. */

double dtmax;
event set_dtmax (i++,last) dtmax = DT;

/**
The macro below can be overloaded to define the barotropic
acceleration. By default it is just the slope of the free-surface
times gravity. */

#define gmetric(i) (2.*fm.x[i]/(cm[i] + cm[i-1]))
#ifndef a_baro
# define a_baro(eta, i) (G*gmetric(i)*(eta[i-1] - eta[i])/Delta)
#endif

/**
## Computation of face values

At each timestep, temporary face fields are defined for the fluxes $(h
\mathbf{u})^{n+1/2}$, face height $h_f^{n+1/2}$ and height-weighted
face accelerations $(ha)^{n+1/2}$. */


static bool hydrostatic = true;
face vector hu, hf, ha;

event face_fields (i++, last)
{
  hu = new face vector[nl];
  hf = new face vector[nl];
  ha = new face vector[nl];

  /**
  The (CFL-limited) timestep is also computed by this function. A
  difficulty is that the prediction step below also requires an
  estimated timestep (the `pdt` variable below). The timestep at the
  previous iteration is used as estimate. */

  static double pdt = 0.;
  foreach_face (reduction (min:dtmax)) {
    double ax = a_baro (eta_r, 0);
    double H = 0., um = 0.;
    double Hr = 0., Hl = 0.;
    foreach_layer() {
      Hr += h[], Hl += h[-1];

      /**
      The face velocity is computed as the height-weighted average of
      the cell velocities. */
      
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      hu.x[] = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;

      /**
      If the left or central cell are dry, we consider a "step-like"
      bathymetry and define the face height as the water level above
      the step. */

      double hff;
      if (Hl <= dry)
	hff = fmax (fmin (zb[] + Hr - zb[-1], h[]), 0.);
      else if (Hr <= dry)
	hff = fmax (fmin (zb[-1] + Hl - zb[], h[-1]), 0.);
      else {

	/**
	In the default case, the face height is computed using a variant of the
	[BCG](/src/bcg.h) scheme. */

	double un = pdt*hu.x[]/Delta, a = sign(un);
	int i = - (a + 1.)/2.;
	double g = h.gradient ? h.gradient (h[i-1], h[i], h[i+1])/Delta :
	  (h[i+1] - h[i-1])/(2.*Delta);
	hff = h[i] + a*(1. - a*un)*g*Delta/2.;
      }
      hf.x[] = fm.x[]*hff;

      /**
      The maximum velocity is stored and the flux and height-weighted
      accelerations are computed. */
      
      if (fabs(hu.x[]) > um)
	um = fabs(hu.x[]);
      
      hu.x[] *= hf.x[];
      ha.x[] = hf.x[]*ax;
 
      H += hff;
    }

    /**
    The maximum timestep is computed using the total depth `H` and the
    advection and gravity wave CFL criteria. The gravity wave speed
    takes dispersion into account in the non-hydrostatic case. */
    
    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
	double dt = min(cm[], cm[-1])*Delta/(c*fm.x[]);
	if (dt < dtmax)
	  dtmax = dt;
      }
    }
  }

  /**
  The timestep is computed, taking into account the timing of events,
  and also stored in `pdt` (see comment above). */
  
  pdt = dt = dtnext (dtmax);
}

/**
## Advection and diffusion

The function below approximates the advection terms using estimates of
the face fluxes $h\mathbf{u}$ and face heights $h_f$. */

void advect (scalar * tracers, face vector hu, face vector hf, double dt)
{

  /**
  The fluxes are first limited according to the CFL condition to
  ensure strict positivity of the layer heights. This step is
  necessary due to the approximate estimation of the CFL condition in
  the timestep calculation above. */
  
  foreach_face()
    foreach_layer() {
      double hul = hu.x[];
      if (hul*dt/(Delta*cm[-1]) > CFL*h[-1])
	hul = CFL*h[-1]*Delta*cm[-1]/dt;
      else if (- hul*dt/(Delta*cm[]) > CFL*h[])
	hul = - CFL*h[]*Delta*cm[]/dt;
      
      /**
      In the case where the flux is limited, and for multiple layers, an
      attempt is made to conserve the total barotropic flux by merging
      the flux difference with the flux in the layer just above. This
      allows to maintain an accurate evolution for the free-surface
      height $\eta$. */
      
      if (hul != hu.x[]) {
	if (point.l < nl - 1)
	  hu.x[0,0,1] += hu.x[] - hul;
#if !_GPU
	else if (nl > 1)
	  fprintf (stderr, "src/layered/hydro.h:%d: warning: could not conserve barotropic flux "
		   "at %g,%g,%d\n", LINENO, x, y, point.l);
#endif
	hu.x[] = hul;
      }
    }

  face vector flux[];
  foreach_layer() {
        
    /**
    We compute the flux $(shu)_{i+1/2,k}$ for each tracer $s$, using a
    variant of the BCG scheme. */
    
    for (scalar s in tracers) {
      foreach_face() {
	double un = dt*hu.x[]/((hf.x[] + dry)*Delta), a = sign(un);
	int i = -(a + 1.)/2.;
	double g = s.gradient ?
	  s.gradient (s[i-1], s[i], s[i+1])/Delta :
	  (s[i+1] - s[i-1])/(2.*Delta);
	double s2 = s[i] + a*(1. - a*un)*g*Delta/2.;

#if dimension > 1
	if (hf.y[i] + hf.y[i,1] > dry) {
	  double vn = (hu.y[i] + hu.y[i,1])/(hf.y[i] + hf.y[i,1]);
	  double syy = (s.gradient ? s.gradient (s[i,-1], s[i], s[i,1]) :
			vn < 0. ? s[i,1] - s[i] : s[i] - s[i,-1]);
	  s2 -= dt*vn*syy/(2.*Delta);
	}
#endif // dimension > 1
	
	flux.x[] = s2*hu.x[];
      }
      
      /**
      We compute $(hs)^\star_i = (hs)^n_i + \Delta t 
      [(shu)_{i+1/2} -(shu)_{i-1/2}]/\Delta$. */
      
      foreach() {
	s[] *= h[];
	foreach_dimension()
	  s[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
      }
    }

    /**
    We then obtain $h^{n+1}$ and $s^{n+1}$ using
    $$
    \begin{aligned}
    h_i^{n+1} & = h_i^n + \Delta t \frac{(hu)_{i+1/2} - (hu)_{i-1/2}}{\Delta},\\
    s_i^{n+1} & = \frac{(hs)^\star_i}{h_i^{n+1}}
    \end{aligned}
    $$
    */

    foreach() {
      double h1 = h[];
      foreach_dimension()
	h1 += dt*(hu.x[] - hu.x[1])/(Delta*cm[]);
#if !_GPU
      if (h1 < - dry)
	fprintf (stderr, "src/layered/hydro.h:%d: warning: h1 = %g < - 1e-12 at %g,%g,%d,%g\n",
		 LINENO, h1, x, y, _layer, t);
#endif
      h[] = max(h1, 0.);
      if (h1 < dry) {
	for (scalar f in tracers)
	  f[] = 0.;
      }
      else
	for (scalar f in tracers)
	  f[] /= h1;      
    }
  }
}

/**
This is where the 'two-step advection' of the [implicit
scheme](implicit.h) plugs itself (nothing is done for the explicit
scheme). */

event half_advection (i++, last);

/**
Vertical diffusion (including viscosity) is added by this code. */

#include "diffusion.h"

/**
## Accelerations 

Acceleration terms are added here. In the simplest case, this is only
the pressure gradient due to the free-surface slope, as computed in
[face_fields](#face_fields). */

event acceleration (i++, last);

event pressure (i++, last)
{

  /**
  The acceleration is applied to the face fluxes... */
  
  foreach_face()
    foreach_layer() 
      hu.x[] += dt*ha.x[];

  /**
  ... and to the centered velocity field, using height-weighting. */
  
  foreach()
    foreach_layer() {
      foreach_dimension()
	u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
#if dimension == 2
      // metric terms
      double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
      double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
      double ux = u.x[], uy = u.y[];
      double fG = uy*dmdl - ux*dmdt;
      u.x[] += dt*fG*uy;
      u.y[] -= dt*fG*ux;
#endif // dimension == 2
    }
  delete ((scalar *){ha});

  /**
  The resulting fluxes are used to advect both tracers and layer
  heights. */
  
  advect (tracers, hu, hf, dt);
}
  
/**
Finally the free-surface height $\eta$ is updated. */

event update_eta (i++, last)
{
  delete ((scalar *){hu, hf});
  foreach() {
    double etap = zb[];
    foreach_layer()
      etap += h[];
    eta[] = etap;
  }
}

/**
## Vertical remapping and mesh adaptation

[Vertical remapping](remap.h) is applied here if necessary. */

event remap (i++, last);
 
#if TREE
event adapt (i++,last);
#endif

/**
## Cleanup

The fields and lists allocated in [`defaults()`](#defaults0) above
must be freed at the end of the run. */
   
event cleanup (t = end, last)
{
  delete ({eta, eta_r, h, u});
  free (tracers), tracers = NULL;
}

/**
# Horizontal pressure gradient

The macro below computes the horizontal pressure gradient
$$
pg_k = - \mathbf{{\nabla}} (h \phi)_k + \left[ \phi 
 \mathbf{{\nabla}} z \right]_k
$$
on the faces of each layer. The slope of the layer interfaces
$\mathbf{{\nabla}} z_{k+1/2}$ in the second-term is bounded by
`max_slope` (by default 30 degrees). */

double max_slope = 0.577350269189626 [0]; // = tan(30.*pi/180.)
#define slope_limited(dz) (fabs(dz) < max_slope ? (dz) :	\
			   ((dz) > 0. ? max_slope : - max_slope))

#define hpg(pg,phi,i,code) do {						\
  double dz = zb[i] - zb[i-1];						\
  foreach_layer() {							\
    double pg = 0.;							\
    if (h[i] + h[i-1] > dry) {						\
      double s = Delta*slope_limited(dz/Delta);				\
      pg = (h[i-1] - s)*phi[i-1] - (h[i] + s)*phi[i];                   \
      if (point.l < nl - 1) {						\
	double s = Delta*slope_limited((dz + h[i] - h[i-1])/Delta);     \
	pg += (h[i-1] + s)*phi[i-1,0,1] - (h[i] - s)*phi[i,0,1];        \
      }									\
      pg *= gmetric(i)*hf.x[i]/(Delta*(h[i] + h[i-1]));			\
    } code;								\
    dz += h[i] - h[i-1];						\
  }									\
} while (0)

/**
# Hydrostatic vertical velocity

For the hydrostatic solver, the vertical velocity is not defined by
default since it is usually not required. The function below can be
applied to compute it using the (diagnostic) incompressibility condition
$$
\mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k + \left[ w -
\mathbf{u} \cdot \mathbf{{\nabla}} (z) \right]_k = 0
$$
*/

void vertical_velocity (scalar w, face vector hu, face vector hf)
{
  foreach() {
    double dz = zb[1] - zb[-1];
    double wm = 0.;
    foreach_layer() {
      w[] = wm + (hu.x[] + hu.x[1])/(hf.x[] + hf.x[1] + dry)*
	(dz + h[1] - h[-1])/(2.*Delta);
      if (point.l > 0)
	foreach_dimension()
	  w[] -= (hu.x[0,0,-1] + hu.x[1,0,-1])
	  /(hf.x[0,0,-1] + hf.x[1,0,-1] + dry)*dz/(2.*Delta);
      foreach_dimension()
	w[] -= (hu.x[1] - hu.x[])/Delta;
      dz += h[1] - h[-1], wm = w[];
    }
  }
}

/**
# "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (*ref*). */

double _radiation (Point point, double ref, scalar s)
{
  double H = 0.;
  foreach_layer()
    H += h[];
#if 0  
  return H > dry ? sqrt(G/H)*(zb[] + H - ref) : 0.;
#else
  return sqrt (G*max(H,0.)) - sqrt(G*max(ref - zb[], 0.));
#endif
}
  
#define radiation(ref) _radiation(point, ref, _s)

/**
# Conservation of water surface elevation

We re-use some generic functions. */

#include "elevation.h"

/**
But we need to re-define the water depth refinement function. */

#if TREE
static void refine_layered_elevation (Point point, scalar h)
{

  /**
  We first check whether we are dealing with "shallow cells". */
  
  bool shallow = zb[] > default_sea_level;
  foreach_child()
    if (zb[] > default_sea_level) {
      shallow = true; break;
    }

  /**
  If we do, refined cells are just set to the default sea level. */
  
  if (shallow)
    foreach_child()
      h[] = max(0., default_sea_level - zb[]);

  /**
  Otherwise, we use the surface elevation of the parent cells to
  reconstruct the water depth of the children cells. */
  
  else {
    double eta = zb[] + h[];   // water surface elevation  
    coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
	g.x = gradient (zb[-1] + h[-1], eta, zb[1] + h[1])/4.;
    else
      foreach_dimension()
	g.x = (zb[1] + h[1] - zb[-1] - h[-1])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta;
      foreach_dimension()
	etac += g.x*child.x;
      h[] = max(0., etac - zb[]);
    }  
  }
}

/**
We overload the `conserve_elevation()` function. */

void conserve_layered_elevation (void)
{
  h.refine  = refine_layered_elevation;
  h.prolongation = prolongation_elevation;
  h.restriction = restriction_elevation;
  h.dirty = true; // boundary conditions need to be updated
}

#define conserve_elevation() conserve_layered_elevation()

#endif // TREE

#include "gauges.h"

#if dimension == 2

/**
# Fluxes through sections

These functions are typically used to compute fluxes (i.e. flow rates)
through cross-sections defined by two endpoints (i.e. segments). Note
that the orientation of the segment is taken into account when
computing the flux i.e the positive normal direction to the segment is
to the "left" when looking from the start to the end. 

This can be expressed mathematically as:
$$
\text{flux}[k] = \int_A^B h_k\mathbf{u}_k\cdot\mathbf{n}dl
$$
with $A$ and $B$ the endpoints of the segment, $k$ the layer index,
$\mathbf{n}$ the oriented segment unit normal and $dl$ the elementary
length. The function returns the sum (over $k$) of all the fluxes. */

double segment_flux (coord segment[2], double * flux, scalar h, vector u)
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  for (int l = 0; l < nl; l++)
    flux[l] = 0.;
  foreach_segment (segment, p, reduction(+:flux[:nl])) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      foreach_layer()
	flux[point.l] += dl/2.*
	interpolate_linear (point, h, a.x, a.y, 0.)*
	(m.x*interpolate_linear (point, u.x, a.x, a.y, 0.) +
	 m.y*interpolate_linear (point, u.y, a.x, a.y, 0.));
    }
  }
  double tot = 0.;
  for (int l = 0; l < nl; l++)
    tot += flux[l];
  return tot;
}

/**
A NULL-terminated array of *Flux* structures passed to
*output_fluxes()* will create a file (called *name*) for each
flux. Each time *output_fluxes()* is called a line will be appended to
the file. The line contains the time, the total flux and the value of
the flux for each $h$, $u$ pair in the layer. The *desc* field can be
filled with a longer description of the flux. */

typedef struct {
  char * name;
  coord s[2];
  char * desc;
  FILE * fp;
} Flux;

void output_fluxes (Flux * fluxes, scalar h, vector u)
{
  for (Flux * f = fluxes; f->name; f++) {
    double flux[nl];
    double tot = segment_flux (f->s, flux, h, u);
    if (pid() == 0) {
      if (!f->fp) {
	f->fp = fopen (f->name, "w");
	if (f->desc)
	  fprintf (f->fp, "%s\n", f->desc);
      }
      fprintf (f->fp, "%g %g", t, tot);
      for (int i = 0; i < nl; i++)
	fprintf (f->fp, " %g", flux[i]);
      fputc ('\n', f->fp);
      fflush (f->fp);
    }
  }
}

#endif // dimension == 2
