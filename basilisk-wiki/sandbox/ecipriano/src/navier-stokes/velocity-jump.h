/**
# Incompressible Navier--Stokes solver with Phase Change Jump Condition

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations with phase change
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity.

We extend this scheme by including the velocity jump due to the phase
change:
$$
  \left(\mathbf{u}_g - \mathbf{u}_l\right)\cdot\mathbf{n}_\Gamma
  = \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)
$$
which is introduced into the velocity field directly, according to the
ghost velocity approach proposed by [Nguyen et al. 2001](#nguyen2001boundary),
that was modified by [Tanguy et al. 2007](#tanguy2007level) for droplet
evaporation problems, and by [Tanguy et al. 2014](#tanguy2014benchmarks)
for boiling simulations.

The idea is to set the velocity jump in the "ghost regions" by exploiting the
jump condition:
$$
  \mathbf{u}^{ghost}_l = \mathbf{u}_g -
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}_\Gamma
$$

$$
  \mathbf{u}^{ghost}_g = \mathbf{u}_l +
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}_\Gamma
$$
Two different momentum equations for the gas
and liquid phase velocities are solved, but a single projection step is
used to update the velocities at the new time step. Additional
projection steps are used to obtain a divergence-free velocity for
the advection of the volume fraction field.

This is a modified version of the method discussed by [Long et al. 2024](#long2024edge)
for phase change using the EBIT approach, and implemented in basilisk by [Tian Long](/sandbox/tianlong/src/semushin/double-evaporation.h).
Unlike the previous version, this approach has mainly been tested on evaporation
problems rather than on boiling simulations, and it was designed for VOF simulations.

The main advantage with respect to the [navier-stokes/centered-evaporation.h](/sandbox/ecipriano/src/navier-stokes/centered-evaporation.h)
approach is that we avoid the oscillations that arise due to the introduction of
the localized expansion term as a volumetric source in the projection step.
*/

#define VELOCITY_JUMP

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
# include "viscosity-embed.h"
#else
# include "viscosity.h"
#endif
#include "aslam.h"

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

/**
Other variables specific to this algorithm. */

vector u1[], u2[], uext[], u1g[], u2g[];
face vector uf1[], uf2[], ufext[], uf1g[], uf2g[];
scalar mEvapTot1[], mEvapTot2[];
scalar ls1[], ls2[];
scalar pg[];
vector n[];
scalar ps[];
face vector nf[];

extern scalar f;
extern scalar mEvapTot;
extern scalar rho1v, rho2v;

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;
scalar drhodt[], drhodtext[];

/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

trace
mgstats project_sf_ghost (face vector uf, face vector ufg, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  foreach()
    div[] += drhodtext[]/dt;

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    ufg.x[] = uf.x[] - dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

trace
mgstats project_extended (face vector uf, face vector ufs, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += ufs.x[1] - ufs.x[];
    div[] /= dt*Delta;
  }

  foreach()
    div[] += drhodtext[]/dt;

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] = ufs.x[] - dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

trace
mgstats project_sf_twofield (face vector uf1, face vector uf2, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    double div1 = 0., div2 = 0.;
    foreach_dimension() {
      div1 += uf1.x[1] - uf1.x[];
      div2 += uf2.x[1] - uf2.x[];
    }
    div[] = (ls1[] < 0.) ? div1 : div2;
    div[] /= dt*Delta;
  }

  /**
  We add the density lagrangian derivative. */

  foreach()
    div[] += drhodt[]/dt;

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face() {
#ifndef DECOUPLED
    uf1.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
#endif
    uf2.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
  }

  return mgp;
}

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :  \
            a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));
ps[right] = neumann (neumann_pressure(ghost));
ps[left]  = neumann (- neumann_pressure(0));
pg[right] = neumann (neumann_pressure(ghost));
pg[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
uf1.n[bottom] = 0.;
uf1.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
uf2.n[bottom] = 0.;
uf2.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
ps[top]    = neumann (neumann_pressure(ghost));
pg[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
ps[top]    = neumann (neumann_pressure(ghost));
ps[bottom] = neumann (- neumann_pressure(0));
pg[top]    = neumann (neumann_pressure(ghost));
pg[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
ps[front]  = neumann (neumann_pressure(ghost));
ps[back]   = neumann (- neumann_pressure(0));
pg[front]  = neumann (neumann_pressure(ghost));
pg[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

#define NOSHIFTING
#ifdef BOILING_SETUP
# define BYRHOGAS
# define CONSISTENTPHASE2
#endif

event defaults (i = 0)
{

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  ps.nodump = pg.nodump = true;

  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;
  uf1.x.refine = refine_face_solenoidal;
  uf2.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  uf1.x.refine = refine_face;
  foreach_dimension() {
    uf.x.prolongation = refine_embed_face_x;
    uf1.x.prolongation = refine_embed_face_x;
    uf2.x.prolongation = refine_embed_face_x;
  }
  for (scalar s in {p, pf, u, g, u1, u2}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf, ps, pg})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE
}

/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  trash ({uf1, uf2, uf});
  foreach_face() {
    uf.x[] = fm.x[]*face_value (u.x, 0);
    uf1.x[] = fm.x[]*face_value (u1.x, 0);
    uf2.x[] = fm.x[]*face_value (u2.x, 0);
  }

  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");

  /**
  We set the default divergence source term to zero (for the liquid phase) */

  foreach() {
    drhodt[] = 0.;
    drhodtext[] = 0.;
  }
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : min (timestep (uf1, dtmax), timestep
        (uf2, dtmax)));
}

/**
## Extrapolations

We use PDE-based Aslam's extrapolations, to extrapolate the
vaporization rate $\dot{m}$, from the interface to the liquid and gas
phases. */

void extrapolations (void)
{
  /**
  First, we store the liquid and gas phase volume fractions. */

  scalar faslam1[], faslam2[];
  foreach() {
    faslam1[] = f[];
    faslam2[] = 1. - f[];
  }

  /**
  We convert the volume fraction in level set. This step can be
  expensive, and it should be skipped if the CLSVOF method is used. */

#ifdef CLSVOF
  extern scalar d;
  foreach() {
    ls1[] = d[];
    ls2[] = -d[];
  }
#else
  vof_to_ls (f, ls1, imax = 5);
  foreach()
    ls2[] = -ls1[];
#endif

  /**
  The interface normals are computed using the level set. By doing so
  the normals are naturally defined over the whole domain. */

  gradients ({ls1}, {n});
  foreach() {
    double maggf = 0.;
    foreach_dimension()
      maggf += sq (n.x[]);
    maggf = sqrt (maggf);
    foreach_dimension()
      n.x[] /= (maggf + 1.e-10);
  }

  /**
  Update the interface normal vectors on the faces. */

  foreach_face()
    nf.x[] = 0.5*(n.x[] + n.x[-1]);

  /**
  We store the extrapolated vaporization rate $\hat{m}$ on the fields
  `mEvapTot1` and `mEvapTot2`. Using `constant_extrapolations` the
  liquid and gas phase cells are populated with the vaporization rate
  without changing its values. */

#ifndef DIFFUSIVE
  foreach() {
    double rhojump = (rho1v[] > 0. && rho2v[] > 0.) ?
      (1./rho2v[] - 1./rho1v[]) : 0.;

    mEvapTot1[] = mEvapTot[]*rhojump;
    mEvapTot2[] = mEvapTot[]*rhojump;
  }

  constant_extrapolation (mEvapTot1, ls1, 0.5, 20, c=faslam1, nl=0,
      nointerface=true);
  constant_extrapolation (mEvapTot1, ls2, 0.5, 20, c=faslam2, nl=0,
      nointerface=true);

  foreach()
    mEvapTot2[] = mEvapTot1[];
#endif
}

/**
We perform the extrapolations after the `phasechange` event in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

event phasechange (i++,last) {
  extrapolations();
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

/**
## Initialize Velocities

We set the initial ghost velocities for the first iteration. For
droplet evaporation problems we use the jump condition to update the
gas phase ghost velocity:

$$
  \mathbf{u}^{ghost}_g = \mathbf{u}_l -
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}_\Gamma
$$

while the initial liquid ghost velocity is set to zero. For boiling
problems we set to zero the initial gas ghost velocity, while we use
the jump condition to set the value of the liquid phase ghost
velocity:

$$
  \mathbf{u}^{ghost}_l = \mathbf{u}_g +
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}_\Gamma
$$
*/

void update_ghost_velocities (void) {
  foreach() {
    foreach_dimension() {
#ifdef BOILING_SETUP
      u2g.x[] = u2g.x[];
      u1g.x[] = u2.x[] + mEvapTot1[]*n.x[];
#else
      u1g.x[] = u1g.x[];
      u2g.x[] = u1.x[] - mEvapTot2[]*n.x[];
#endif
      u1.x[] = (ls1[] < 0.) ? u1.x[] : u1g.x[];
      u2.x[] = (ls1[] > 0.) ? u2.x[] : u2g.x[];
    }
  }

  foreach_face() {
#ifdef BOILING_SETUP
    double mEvapTot1f = 0.5*(mEvapTot1[] + mEvapTot1[-1]);
    uf2g.x[] = uf2g.x[];
    uf1g.x[] = uf2.x[] + mEvapTot1f*nf.x[]*fm.x[];
#else
    double mEvapTot2f = 0.5*(mEvapTot2[] + mEvapTot2[-1]);
    uf1g.x[] = uf1g.x[];
    uf2g.x[] = uf1.x[] - mEvapTot2f*nf.x[]*fm.x[];
#endif
    double lsf = 0.5*(ls1[] + ls1[-1]);
    uf1.x[] = (lsf < 0.) ? uf1.x[] : uf1g.x[];
    uf2.x[] = (lsf > 0.) ? uf2.x[] : uf2g.x[];
  }
}

event init_ghost (i = 0, last)
{
  /**
  We compute the value of the centered ghost velocities. */

  foreach() {
    foreach_dimension() {
#ifdef BOILING_SETUP
      u2g.x[] = 0.;
      u1g.x[] = u2.x[] + mEvapTot1[]*n.x[];
#else
      u1g.x[] = 0.;
      u2g.x[] = u1.x[] - mEvapTot2[]*n.x[];
#endif
    }
  }

  /**
  We compute the value of the face ghost velocities. The interface
  normal on the face is computed using the level set. */

  foreach_face() {
#ifdef BOILING_SETUP
    double mEvapTot1f = 0.5*(mEvapTot1[] + mEvapTot1[-1]);
    uf2g.x[] = 0.;
    uf1g.x[] = uf2.x[] + mEvapTot1f*nf.x[]*fm.x[];
#else
    double mEvapTot2f = 0.5*(mEvapTot2[] + mEvapTot2[-1]);
    uf1g.x[] = 0.;
    uf2g.x[] = uf1.x[] - mEvapTot2f*nf.x[]*fm.x[];
#endif
  }
}

/**
At the beginning of every iteration we initialize the liquid
and gas phase velocities using the ghost values:

$$
  \mathbf{u}_l =
  \begin{cases}
    \mathbf{u}_l  & \text{if } \phi < 0 \\
    \mathbf{u}_l^{ghost} & \text{if } \phi > 0
  \end{cases}
$$

$$
  \mathbf{u}_g =
  \begin{cases}
    \mathbf{u}_g  & \text{if } \phi > 0 \\
    \mathbf{u}_g^{ghost} & \text{if } \phi < 0
  \end{cases}
$$

Where the level set $\phi$ is negative in the liquid phase and
positive in the gas phase.
*/

event init_velocities (i++, last)
{
  foreach() {
    foreach_dimension() {
      u1.x[] = (ls1[] < 0.) ? u1.x[] : u1g.x[];
      u2.x[] = (ls1[] > 0.) ? u2.x[] : u2g.x[];
    }
  }

  foreach_face() {
    double lsf = 0.5*(ls1[] + ls1[-1]);
    uf1.x[] = (lsf < 0.) ? uf1.x[] : uf1g.x[];
    uf2.x[] = (lsf > 0.) ? uf2.x[] : uf2g.x[];
  }
}

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction(vector u, face vector uf)
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
    du.x[] = 0.;
  else
#endif
    du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
    du.x[] = 0.;
  else
#endif
    du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar *){du});
}

/**
We solve the advection equations for $\mathbf{u}_l$ and
$\mathbf{u}_g$. */

void advection_div (scalar * tracers, face vector u, double dt,
		scalar * src = NULL)
{
  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    tracer_fluxes (f, u, flux, dt, source);
#if !EMBED
    foreach() {
#if NO_ADVECTION_DIV
      double fold = f[];
#endif
      foreach_dimension()
#if NO_ADVECTION_DIV
        f[] += dt*(flux.x[] - flux.x[1] + fold*(u.x[1] - u.x[]))/(Delta*cm[]);
#else
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#endif
    }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}

event advection_term (i++, last)
{
  if (!stokes) {
    prediction (u1, uf1);
    prediction (u2, uf2);
    mgpf = project_sf_twofield (uf1, uf2, pf, alpha, dt/2., mgpf.nrelax);
    advection_div ((scalar *){u1}, uf1, dt, (scalar *){g});
    advection_div ((scalar *){u2}, uf2, dt, (scalar *){g});
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (vector u, double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
}

/**
Solving the viscous term we obtain the temporary velocities
$\mathbf{u}_l^*$ and $\mathbf{u}_g^*$. */

event viscous_term (i++, last)
{
  if (constant(mu.x) != 0.) {
    correction (u1, dt);
    correction (u2, dt);
    mgu = viscosity (u1, mu, rho, dt, mgu.nrelax);
    mgu = viscosity (u2, mu, rho, dt, mgu.nrelax);
    correction (u1, -dt);
    correction (u2, -dt);
  }

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  trash ({uf1,uf2});
  foreach_face() {
    uf1.x[] = fm.x[]*(face_value (u1.x, 0) + dt*a.x[]);
    uf2.x[] = fm.x[]*(face_value (u2.x, 0) + dt*a.x[]);
  }
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

/**
## Projection

To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  /**
  For boiling problems, we obtain a ghost pressure $p^{ghost}$ using
  the following projection step: 
  $$
    \nabla\cdot\left(\dfrac{1}{\rho}\nabla p^{ghost}\right) =
    \dfrac{\nabla\cdot\mathbf{u}_g^*}{\Delta t}
  $$
  from which we obtain the updated gas ghost velocity:
  $$
    \mathbf{u}_g^{ghost} = \mathbf{u}_g^*
    - \dfrac{\Delta t}{\rho}\nabla p^{ghost}
  $$

  For evaporation problems, we solve the ghost pressure using the
  divergence of the liquid velocity:
  $$
    \nabla\cdot\left(\dfrac{1}{\rho}\nabla p^{ghost}\right) =
    \dfrac{\nabla\cdot\mathbf{u}_l^*}{\Delta t}
  $$
  and we use this pressure to update the liquid ghost velocity:
  $$
    \mathbf{u}_l^{ghost} = \mathbf{u}_l^*
    - \dfrac{\Delta t}{\rho}\nabla p^{ghost}
  $$
  */

#ifdef BOILING_SETUP
  project_sf_ghost (uf2, uf2g, pg, alpha, dt, mgpf.nrelax);

  vector gpg[];
  centered_gradient (pg, gpg);
  foreach()
    foreach_dimension()
      u2g.x[] = u2.x[] + dt*gpg.x[];
#else
  project_sf_ghost (uf1, uf1g, pg, alpha, dt, mgpf.nrelax);

  vector gpg[];
  centered_gradient (pg, gpg);
  foreach()
    foreach_dimension()
      u1g.x[] = u1.x[] + dt*gpg.x[];
#endif

  /**
  The one-field pressure of the system is obtained from the following
  projection step:
  $$
    \nabla\cdot\left(\dfrac{1}{\rho}\nabla p^{ghost}\right) =
    \dfrac{\nabla\cdot\mathbf{u}^*}{\Delta t}
  $$

  where the divergence is computed based on the level set field:
  $$
    \nabla\cdot\mathbf{u}^* =
    \begin{cases}
      \nabla\cdot\mathbf{u}_l^* & \text{if } \phi < 0 \\
      \nabla\cdot\mathbf{u}_g^* & \text{if } \phi > 0
    \end{cases}
  $$

  and the two velocities are computed from the one-field pressure and
  update at the new time step as:

  $$
    \mathbf{u}_l^{n+1} = \mathbf{u}_l^* - \dfrac{\Delta t}{\rho}
    \nabla p
  $$
  $$
    \mathbf{u}_g^{n+1} = \mathbf{u}_g^* - \dfrac{\Delta t}{\rho}
    \nabla p
  $$
  */

  mgp = project_sf_twofield (uf1, uf2, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  correction (u1, dt);
  correction (u2, dt);
}

/**
We update the centered and face velocity values. */

event update_ghost (i++, last) {
  foreach() {
    foreach_dimension() {
#ifdef BOILING_SETUP
      u2g.x[] = u2g.x[];
      u1g.x[] = u2.x[] + mEvapTot1[]*n.x[];
#else
      u1g.x[] = u1g.x[];
      u2g.x[] = u1.x[] - mEvapTot2[]*n.x[];
#endif
    }
  }

  foreach_face() {
#ifdef BOILING_SETUP
    double mEvapTot1f = 0.5*(mEvapTot1[] + mEvapTot1[-1]);
    uf2g.x[] = uf2g.x[];
    uf1g.x[] = uf2.x[] + mEvapTot1f*nf.x[]*fm.x[];
#else
    double mEvapTot2f = 0.5*(mEvapTot2[] + mEvapTot2[-1]);
    uf1g.x[] = uf1g.x[];
    uf2g.x[] = uf1.x[] - mEvapTot2f*nf.x[]*fm.x[];
#endif
  }
}

/**
## Extended Velocity

According to [Tanguy et al. 2007](#tanguy2007level), an additional
projection step is required in evaporation problems to ensure that
the liquid velocity is conservative. This step could be skipped for
boiling problems. We define a velocity $\mathbf{u}^S$ as:
$$
  \mathbf{u}^S =
  \begin{cases}
    \mathbf{u}_l^{n+1} & \text{if } \phi < 0 \\
    \mathbf{u}_l^{ghost} & \text{if } \phi > 0
  \end{cases}
$$

And we solve an additional Poisson equation:
$$
  \nabla\cdot\left(\dfrac{1}{\rho}\nabla\psi\right) =
  \dfrac{\nabla\cdot\mathbf{u}^S}{\Delta t}
$$
from which we correct the liquid velocity ensuring that it respects
the divergence-free constraint:
$$
  \mathbf{u}_l^{n+1} = \mathbf{u}^S - \dfrac{\Delta t}{\rho}\nabla \psi
$$

If we perform this step also for boiling problems we need to do the
same operation replacing $\mathbf{u}_l$ with $\mathbf{u}_g$, assuming
that the VOF field is advected using the gas phase velocity. */

event extended_velocity (i++, last) {
#ifdef BOILING_SETUP
  face vector ufs[];
  foreach_face() {
    double lsf = 0.5*(ls1[] + ls1[-1]);
    ufs.x[] = (lsf > 0.) ? uf2.x[] : uf2g.x[];
  }
  project_extended (uf2, ufs, ps, alpha, dt, mgp.nrelax);

  foreach_face()
    uf2g.x[] = uf2.x[];
#else
  face vector ufs[];
  foreach_face() {
    double lsf = 0.5*(ls1[] + ls1[-1]);
    ufs.x[] = (lsf < 0.) ? uf1.x[] : uf1g.x[];
  }
  project_extended (uf1, ufs, ps, alpha, dt, mgp.nrelax);

  foreach_face()
    uf1g.x[] = uf1.x[];
#endif
}

/**
## Reconstructions

We reconstruct $\mathbf{u}_l^{n+1}$ and $\mathbf{u}_g^{n+1}$ using
the ghost liquid and gas phase velocities, for the next time step.
The one-field velocity $\mathbf{u}$ is simply reconstructed from a
volume average, just for visualization purposes. 

We also set the value of the extended velocity $\mathbf{u}^E$, used
for the VOF advection equation. */

event reconstructions (i++, last) {
  foreach() {
    foreach_dimension() {
      u1.x[] = (ls1[] < 0.) ? u1.x[] : u1g.x[];
      u2.x[] = (ls1[] > 0.) ? u2.x[] : u2g.x[];
    }
  }

  foreach_face() {
    double lsf = 0.5*(ls1[] + ls1[-1]);
    uf1.x[] = (lsf < 0.) ? uf1.x[] : uf1g.x[];
    uf2.x[] = (lsf > 0.) ? uf2.x[] : uf2g.x[];
#ifdef BOILING_SETUP
    ufext.x[] = uf2.x[];
#else
    ufext.x[] = uf1.x[];
#endif
  }

  foreach()
    foreach_dimension()
      u.x[] = u1.x[]*f[] + u2.x[]*(1. - f[]);

  foreach_face() {
    double ff = 0.5*(f[] + f[-1]);
    uf.x[] = uf1.x[]*ff + uf2.x[]*(1. - ff);
  }
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
#endif
  event ("properties");
}
#endif

/**
## References

~~~bib
@article{tanguy2007level,
  title={A level set method for vaporizing two-phase flows},
  author={Tanguy, S{\'e}bastien and M{\'e}nard, Thibaut and Berlemont, Alain},
  journal={Journal of Computational Physics},
  volume={221},
  number={2},
  pages={837--853},
  year={2007},
  publisher={Elsevier}
}

@article{tanguy2014benchmarks,
  title={Benchmarks and numerical methods for the simulation of boiling flows},
  author={Tanguy, S{\'e}bastien and Sagan, Micha{\"e}l and Lalanne, Benjamin and Couderc, Fr{\'e}d{\'e}ric and Colin, Catherine},
  journal={Journal of Computational Physics},
  volume={264},
  pages={1--22},
  year={2014},
  publisher={Elsevier}
}

@article{nguyen2001boundary,
  title={A boundary condition capturing method for incompressible flame discontinuities},
  author={Nguyen, Duc Q and Fedkiw, Ronald P and Kang, Myungjoo},
  journal={Journal of Computational Physics},
  volume={172},
  number={1},
  pages={71--98},
  year={2001},
  publisher={Elsevier}
}

@article{long2024edge,
  title={An Edge-based Interface Tracking (EBIT) Method for Multiphase Flows with Phase Change},
  author={Long, Tian and Pan, Jieyun and Zaleski, St{\'e}phane},
  journal={arXiv preprint arXiv:2402.13677},
  year={2024}
}
~~~
*/

