/**
# Momentum-conserving advection of velocity

This file implements momentum-conserving VOF advection of the velocity
components for the [two-phase](/src/two-phase.h) Navier--Stokes
solver.

On trees, we first define refinement and restriction functions which
guarantee conservation of each component of the total momentum. Note
that these functions do not guarantee conservation of momentum for
each phase. */

#if TREE
static void momentum_refine (Point point, scalar u) {
  refine_bilinear (point, u);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f[])*u[];
  double du = u[] - rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[]));
  foreach_child()
    u[] += du;
}

static void momentum_restriction (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f[])*u[];
  u[] = rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[]));
}
#endif // TREE

/**
We switch-off the default advection scheme of the [centered
solver](centered.h). */

event defaults (i = 0)
{
  stokes = true;

#if TREE

  /**
  On trees, the refinement and restriction functions above rely on the
  volume fraction field *f* being refined/restricted before the
  components of velocity. To ensure this, we move *f* to the front of
  the field list (*all*). */

  int i = 0;
  while (all[i].i != f.i) i++;
  while (i > 0 && all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f;
    
  /**
  We then set the refinement and restriction functions for the
  components of the velocity field. The boundary conditions on
  $\mathbf{u}$ now depend on those on $f$. */
  
  foreach_dimension() {
    u.x.refine = u.x.prolongation = momentum_refine;
    u.x.restriction = momentum_restriction;
    u.x.depends = list_add (u.x.depends, f);
  }
#endif
}

/**
We need to overload the stability event so that the CFL is taken into
account (because we set stokes to true). */

event stability (i++)
  dtmax = timestep (uf, dtmax);

/**
We will transport the two components of the momentum, $q_1=f \rho_1
\mathbf{u}$ and $q_2=(1 - f) \rho_2 \mathbf{u}$. We will need to
impose boundary conditions which match this definition. This is done
using the functions below. */

foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q1, bool * data)
{
  return clamp(f[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static double boundary_q2_x (Point neighbor, Point point, scalar q2, bool * data)
{
  return (1. - clamp(f[],0.,1.))*rho2*u.x[];
}

/**
Similarly, on trees we need prolongation functions which also follow
this definition. */

#if TREE
foreach_dimension()
static void prolongation_q1_x (Point point, scalar q1) {
  foreach_child()
    q1[] = clamp(f[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static void prolongation_q2_x (Point point, scalar q2) {
  foreach_child()
    q2[] = (1. - clamp(f[],0.,1.))*rho2*u.x[];
}
#endif

/**
We overload the *vof()* event to transport consistently the volume
fraction and the momentum of each phase. */

static scalar * interfaces1 = NULL;

event vof (i++) {

  /**
  We allocate two temporary vector fields to store the two components
  of the momentum and set the boundary conditions and prolongation
  functions. The boundary conditions on $q_1$ and $q_2$ depend on the
  boundary conditions on $f$. */
  
  vector q1[], q2[];
  for (scalar s in {q1,q2}) {
    s.depends = list_add (s.depends, f);
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  }
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1.x.boundary[i] = boundary_q1_x;
      q2.x.boundary[i] = boundary_q2_x;
    }
#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2.x.prolongation = prolongation_q2_x;
  }
#endif

  /**
  We split the total momentum $q$ into its two components $q1$ and
  $q2$ associated with $f$ and $1 - f$ respectively. */

  foreach()
    foreach_dimension() {
      double fc = clamp(f[],0,1);
      q1.x[] = fc*rho1*u.x[];
      q2.x[] = (1. - fc)*rho2*u.x[];
    }

  /**
  Momentum $q2$ is associated with $1 - f$, so we set the *inverse*
  attribute to *true*. We use the same slope-limiting as for the
  velocity field. */

  foreach_dimension() {
    q2.x.inverse = true;
    q1.x.gradient = q2.x.gradient = u.x.gradient;
  }

  /**
  We associate the transport of $q1$ and $q2$ with $f$ and transport
  all fields consistently using the VOF scheme. */

  scalar * tracers = f.tracers;
  f.tracers = list_concat (tracers, (scalar *){q1, q2});
  vof_advection ({f}, i);
  free (f.tracers);
  f.tracers = tracers;
  
  /**
  We recover the advected velocity field using the total momentum and
  the density */

  foreach()
    foreach_dimension()
      u.x[] = (q1.x[] + q2.x[])/rho(f[]);

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */
  
  interfaces1 = interfaces, interfaces = NULL;
}

/**
We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
}
