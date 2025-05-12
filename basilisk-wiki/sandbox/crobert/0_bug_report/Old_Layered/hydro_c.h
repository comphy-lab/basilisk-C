/**
# A vertically-Lagrangian multilayer solver for free-surface flows

This is the implementation of the solver described in [Popinet,
2019](/Bibliography#popinet2019). The system of $n$ layers of
incompressible fluid pictured below is modelled as the set of
(hydrostatic) equations:
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
$w_k$ can be added using [the non-hydrostatic extension](nh.h).

![Definition of the $n$ layers.](../figures/layers.svg){ width="60%" }

## Fields and parameters

The number of layers is stored in `nl`. Each scalar field has a new
attribute `l` which can be used to access the index of the
corresponding layer.

The `zb` and `eta` fields define the bathymetry and free-surface
height. The layer thicknesses and velocities for each layer are
defined in the `hl` and `ul` scalar and vector lists.

`tracers` is an array of lists of tracers for each layer. By default
it contains only the components of velocity.

The acceleration of gravity is `G` and `dry` controls the minimum
layer thickness.

The `gradient` pointer gives the function used for limiting.

The `ufl` and `al` lists contain the face velocity and face
acceleration fields for each layer. */

#include "run.h"
#include "bcg.h"

int nl = 1;
attribute { int l; }
scalar zb[], eta;
scalar * hl = NULL, ** tracers = NULL;
vector * ul = NULL;
scalar * cl = NULL;
double G = 1., dry = 1e-6;
double (* gradient) (double, double, double) = minmod2;

vector * ufl = NULL, * al = NULL;

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
    for (scalar h in hl)
      eta[] += h[];
  }
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[];
  for (scalar h in hl)
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
  assert (hl == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar h = new scalar;
    h.gradient = gradient;
#if TREE    
    h.refine = h.prolongation = refine_linear;
    h.restriction = restriction_volume_average;
#endif
    hl = list_append (hl, h);
    h.l = l;
  }
  eta = new scalar;
  reset (hl, 0.);
  reset ({zb}, 0.);

  /**
  We set the proper gradient and refinement/restriction functions. */
  
  zb.gradient = gradient;
  eta.gradient = gradient;
#if TREE
  zb.refine = zb.prolongation = refine_linear;
  zb.restriction = restriction_volume_average;
  eta.prolongation = refine_linear;
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
#endif // TREE

  assert (!tracers);
  tracers = calloc (nl, sizeof(scalar *));
}

/**
Other fields, such as $\mathbf{u}_k$ here, are added by this
event. */

event defaults (i = 0)
{
  assert (ul == NULL && ufl == NULL && al == NULL);
  for (int l = 0; l < nl; l++) {
    vector u = new vector;
    vector uf = new face vector;
    vector a = new face vector;
    ul = vectors_append (ul, u);
    ufl = vectors_append (ufl, uf);
    al = vectors_append (al, a);
    foreach_dimension() {
      u.x.l = uf.x.l = a.x.l = l;
      tracers[l] = list_append (tracers[l], u.x);
    }
  }
  reset (ul, 0.);
  reset (ufl, 0.);
  reset (al, 0.);

  /**
  The gradient and prolongation/restriction functions are set for all
  tracer fields. */

  for (int l = 0; l < nl; l++)
    for (scalar s in tracers[l]) {
      s.gradient = gradient;
#if TREE
      s.refine = s.prolongation = refine_linear;
      s.restriction = restriction_volume_average;
#endif
    }
}

/**
After user initialisation, we define the initial face velocities `uf`
and free-surface height $\eta$. */

double dtmax;

event init (i = 0)
{
  trash ((scalar *)ufl);
  foreach_face() {
    vector uf, u;
    scalar h;
    for (h,u,uf in hl,ul,ufl)
      uf.x[] = fm.x[]*(h[]*u.x[] + h[-1]*u.x[-1])/(h[] + h[-1] + dry);
  }
  boundary ((scalar *)ufl);
 
  foreach() {
    eta[] = zb[];
    for (scalar h in hl)
      eta[] += h[];
  }
  boundary (all);

  dtmax = DT;
  event ("stability");
}

/**
## Stability condition

The maximum timestep is set using the standard (hydrostatic) CFL condition
$$
\Delta t < \text{CFL}\frac{\Delta}{|\mathbf{u}|+\sqrt{gH}}
$$
or the modified (non-hydrostatic) condition ([Popinet,
2019](/Bibliography#popinet2019))
$$
\Delta t < 
\text{CFL}\frac{\Delta}{|\mathbf{u}|+\sqrt{g\Delta H\text{tanh}(H/\Delta)}}
$$
*/

#define face_is_wet() ((H > dry && Hm > dry) ||				\
		       (H > dry && eta[] >= zb[-1]) ||			\
		       (Hm > dry && eta[-1] >= zb[]))

event set_dtmax (i++,last) dtmax = DT;

static bool hydrostatic = true;

event stability (i++,last)
{
  foreach_face (reduction (min:dtmax)) {
    double H = 0., Hm = 0.;
#if !TREE
    for (scalar h in hl)
      H += h[], Hm += h[-1];
#else // TREE
    scalar h;
    vector uf;
    for (h,uf in hl,ufl) {
      H += h[], Hm += h[-1];
    
      /**
      We check and enforce that fluxes between dry and wet cells are
      zero. This is necessary only when using adaptive refinement,
      otherwise this is guaranteed by the way face velocities are
      computed in the [acceleration section](#acceleration). */
      
      if ((h[] < dry && uf.x[] < 0.) ||
	  (h[-1] < dry && uf.x[] > 0.))
	uf.x[] = 0.;
    }
    if (!face_is_wet()) {
      for (vector uf in ufl)
	uf.x[] = 0.;
    }
    else
#endif // TREE

    if (H + Hm > 0.) {
      H = (H + Hm)/2.;
      double cp = hydrostatic ? sqrt(G*H) : sqrt(G*Delta*tanh(H/Delta));
      for (vector uf in ufl) {
	double c = fm.x[]*cp + fabs(uf.x[]);
	if (c > 0.) {
	  double dt = cm[]*Delta/c;
	  if (dt < dtmax)
	    dtmax = dt;
	}
      }
    }
  }
  dt = dtnext (CFL*dtmax);
}

/**
## Advection (and diffusion)

This section first solves
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = 0, \\
  \partial_t \left( h s \right)_k + \mathbf{{\nabla}} \cdot \left(
  h s  \mathbf{u} \right)_k & = 0,
\end{aligned}
$$
where $s_k$ is any other tracer field added to the `tracers` list. */

event advection_term (i++,last)
{
  face vector F[], flux[];
  
  vector uf;
  scalar h;
  int l = 0;
  for (uf,h in ufl,hl) {

    /**
    We first compute the "thickness flux" $F_{i+1/2,k}=(hu)_{i+1/2,k}$
    using the [Bell--Collela--Glaz](/src/bcg.h) scheme. */

    tracer_fluxes (h, uf, F, dt, zeroc);

    /**
    We then compute the flux $(sF)_{i+1/2,k}$ for each tracer $s$, also
    using a variant of the BCG scheme. */
    
    for (scalar s in tracers[l]) {
      foreach_face() {
	double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), a = sign(un);
	int i = -(a + 1.)/2.;
	double g = s.gradient ?
	  s.gradient (s[i-1], s[i], s[i+1])/Delta :
	  (s[i+1] - s[i-1])/(2.*Delta);
	double s2 = s[i] + a*(1. - a*un)*g*Delta/2.;
#if dimension > 1
	if (fm.y[i] && fm.y[i,1]) {
	  double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
	  double syy = (s.gradient ? s.gradient (s[i,-1], s[i], s[i,1]) :
			vn < 0. ? s[i,1] - s[i] : s[i] - s[i,-1]);
	  s2 -= dt*vn*syy/(2.*Delta);
	}
#endif
	flux.x[] = s2*F.x[];
      }
      boundary_flux ({flux});

      /**
      We compute $(hs)^\star_i = (hs)^n_i + \Delta t 
      [(sF)_{i+1/2} -(sF)_{i-1/2}]/\Delta$. */
      
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
    h_i^{n+1} & = h_i^n + \Delta t \frac{F_{i+1/2} - F_{i-1/2}}{\Delta}, \\
    s_i^{n+1} & = \frac{(hs)^\star_i}{h_i^{n+1}}
    \end{aligned}
    $$
    */

    foreach() {
      foreach_dimension()
	h[] += dt*(F.x[] - F.x[1])/(Delta*cm[]);
      if (h[] < dry) {
	for (scalar f in tracers[l])
	  f[] = 0.;
      }
      else
	for (scalar f in tracers[l])
	  f[] /= h[];
    }
    l++;
  }

  /**
  Finally the free-surface height $\eta$ is updated and the boundary
  conditions are applied. */
  
  foreach() {
    double H = 0.;
    for (scalar h in hl)
      H += h[];
    eta[] = zb[] + H;
  }

  scalar * list = list_copy (hl);
  for (int l = 0; l < nl; l++)
    for (scalar s in tracers[l])
      list = list_append (list, s);
  list = list_append (list, eta);
  boundary (list);
  free (list);
}
 
/**
Vertical diffusion (including viscosity) is added by this code. */

#include "diffusion_visco.h"
#include "diffusion_solute_vert.h"
#include "diffusion_solute_hor.h"
/**
## Acceleration

We compute the acceleration as
$$
a_{i + 1 / 2, k} = - g \frac{\eta_{i + 1, k} - \eta_{i, k}}{\Delta}
$$
taking into account [metric
terms](/src/README#general-orthogonal-coordinates), and the face
velocity field as
$$
  u_{i + 1 / 2, k} = \frac{(hu)_{i + 1, k} + 
  (hu)_{i, k}}{h_{i + 1, k} + h_{i, k}} + \Delta ta_{i + 1 / 2, k}
$$
provided the "wet condition"
$$
\begin{aligned}
  H_i > \epsilon \quad &\text{and} \quad H_{i + 1} > \epsilon 
  \quad \text{or} \\
  H_i > \epsilon \quad &\text{and} \quad \eta_i \leq z_{i + 1, b} 
  \quad \text{or} \\
  H_{i + 1} > \epsilon \quad &\text{and} \quad \eta_{i + 1} \leq z_{i, b}
\end{aligned}
$$
is verified. */

event acceleration (i++,last)
{
  trash (ufl);
  foreach_face() {
    double H = 0., Hm = 0.;
    for (scalar h in hl)
      H += h[], Hm += h[-1];
    scalar h;
    face vector a, uf;
    vector u;
    if (face_is_wet()) {
      for (h,a,uf,u in hl,al,ufl,ul) {
        a.x[] = - sq(fm.x[])*G*(eta[] - eta[-1])/((cm[] + cm[-1])*Delta/2.);
        uf.x[] = fm.x[]*(h[]*u.x[] + h[-1]*u.x[-1])/(h[] + h[-1] + dry) +
	  dt*a.x[];
      }
    }
    else {
      for (a,uf in al,ufl)
	a.x[] = uf.x[] = 0.;
    }
  }
  boundary ((scalar *)ufl);
  boundary ((scalar *)al);
}

/**
Finally the acceleration is used to obtain the velocity field at time
$n+1$ as
$$
u^{n + 1}_{i, k} \leftarrow u^{n + 1}_{i, k} + \Delta t \frac{a_{i + 1 / 2,
k} + a_{i - 1 / 2, k}}{2}
$$
The missing shallow-water metric terms $(f_G v, -f_G u)$, with
$$
f_G \equiv \frac{v \partial_{\lambda} m_{\theta} - u \partial_{\theta}
m_{\lambda}}{m_{\lambda} m_{\theta}}
$$
are also added here ([Popinet, 2011](/src/references.bib#popinet2011)). */

event pressure (i++,last)
{
  foreach() {
    vector a, u;
    for (a,u in al,ul) {
      foreach_dimension()
        u.x[] += dt*(a.x[] + a.x[1])/(fm.x[] + fm.x[1]);
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
  }
  boundary ((scalar *) ul);
}

/**
## Mesh adaptation and cleanup

Adaptation plugs itself here. The fields and lists allocated in
[`defaults()`](#defaults0) above must be freed at the end of the run. */

#if TREE
event adapt (i++,last);
#endif

event cleanup (i = end, last)
{
  for (int l = 0; l < nl; l++)
    free (tracers[l]);
  free (tracers), tracers = NULL;
  delete ((scalar *) ul), free (ul), ul = NULL;
  delete ((scalar *) ufl), free (ufl), ufl = NULL;
  delete ((scalar *) al), free (al), al = NULL;
  delete (hl), free (hl), hl = NULL;
  delete ({eta});
}

/**
# Hydrostatic vertical velocity

For the hydrostatic solver, the vertical velocity is not defined by
default since it is usually not required. The function below can be
applied to compute it using the (diagnostic) incompressibility condition
$$
\mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k + \left[ w -
\mathbf{u} \cdot \mathbf{{\nabla}} (z) \right]_k = 0
$$
Note that the `wl` list must first be allocated separately. */

void vertical_velocity (scalar * wl)
{
  foreach() {
    scalar h, w;
    vector uf;
    int l = 0;
    double dz = zb[1] - zb[-1];
    double wm = 0.;
    for (h,w,uf in hl,wl,ufl) {
      w[] = wm + (uf.x[] + uf.x[1])*(dz + h[1] - h[-1])/(4.*Delta);
      if (l > 0) {
	vector ufm = ufl[l-1];
	w[] -= (ufm.x[] + ufm.x[1])*dz/(4.*Delta);
      }
      foreach_dimension()
	w[] -= ((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      l++, dz += h[1] - h[-1], wm = w[];
    }
  }
  boundary (wl);
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
  for (scalar h in hl)
    H += h[];
  return H > dry ? sqrt(G/H)*(zb[] + H - ref) : 0.;
}
  
#define radiation(ref) _radiation(point, ref, _s)

#include "elevation.h"
#include "gauges.h"
