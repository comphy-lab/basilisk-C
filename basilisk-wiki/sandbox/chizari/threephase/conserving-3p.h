/**
# Momentum-conserving advection of velocity

This file implements momentum-conserving VOF advection of the velocity
components for the three-phase Navier--Stokes
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
    rhou += cm[]*rho(1.0-f2[]-f3[],f2[],f3[])*u[];
  double du = u[] - rhou/((1 << dimension)*cm[]*rho(1.0-f2[]-f3[],f2[],f3[]));
  foreach_child()
    u[] += du;
}

static void momentum_restriction (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(1.0-f2[]-f3[],f2[],f3[])*u[];
  u[] = rhou/((1 << dimension)*cm[]*rho(1.0-f2[]-f3[],f2[],f3[]));
}
#endif // TREE

/**
We switch-off the default advection scheme of the [centered
solver](centered.h). */

event defaults (i = 0) {
  stokes = true;

#if TREE

  /**
  On trees, the refinement and restriction functions above rely on the
  volume fraction field *f* being refined/restricted before the
  components of velocity. To ensure this, we move *f* to the front of
  the field list (*all*). */

  int i = 0;
  while (all[i].i != f3.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f3;

  i = 0;
  while (all[i].i != f2.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f2;

  i = 0;
  while (all[i].i != f1.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f1;

  /**
  We then set the refinement and restriction functions for the
  components of the velocity field. */
  
  foreach_dimension() {
    u.x.refine = u.x.prolongation = momentum_refine;
    u.x.restriction = momentum_restriction;
  }
#endif
}

/**
We will transport the two components of the momentum, $q_i=f_i \rho_i
\mathbf{u}$. We will need to
impose boundary conditions which match this definition. This is done
using the functions below. */

foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q1)
{
  return clamp(f1[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static double boundary_q2_x (Point neighbor, Point point, scalar q2)
{
  return clamp(f2[],0.,1.)*rho2*u.x[];
}

foreach_dimension()
static double boundary_q3_x (Point neighbor, Point point, scalar q3)
{
  return clamp(f3[],0.,1.)*rho3*u.x[];
}

/**
Similarly, on trees we need prolongation functions which also follow
this definition. */

#if TREE
foreach_dimension()
static void prolongation_q1_x (Point point, scalar q1) {
  foreach_child()
    q1[] = clamp(f1[],0.,1.)*rho1*u.x[];
}

foreach_dimension()
static void prolongation_q2_x (Point point, scalar q2) {
  foreach_child()
    q2[] = clamp(f2[],0.,1.)*rho2*u.x[];
}

foreach_dimension()
static void prolongation_q3_x (Point point, scalar q3) {
  foreach_child()
    q3[] = clamp(f3[],0.,1.)*rho3*u.x[];
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
  functions. */
  
  vector q1[], q2[], q3[];
  for (scalar s in {q1,q2,q3})
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1.x.boundary[i] = boundary_q1_x;
      q2.x.boundary[i] = boundary_q2_x;
      q3.x.boundary[i] = boundary_q3_x;
    }
#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2.x.prolongation = prolongation_q2_x;
    q3.x.prolongation = prolongation_q3_x;
  }
#endif

  /**
  We split the total momentum $q$ into its two components $q1, q2$ and
  $q3$ associated with $f1, f2$ and $f3$, respectively. */

  foreach()
    foreach_dimension() {
      q1.x[] = clamp(f1[],0,1)*rho1*u.x[];
      q2.x[] = clamp(f2[],0,1)*rho2*u.x[];
      q3.x[] = clamp(f3[],0,1)*rho3*u.x[];
    }
  boundary ((scalar *){q1,q2,q3});

  /**
  We use the same slope-limiting as for the
  velocity field. */

  foreach_dimension() {
    q1.x.gradient = q2.x.gradient = q3.x.gradient = u.x.gradient;
  }

  /**
  We associate the transport of $q1, q2, q3$ with $f1, f2, f3$ and transport
  all fields consistently using the VOF scheme. */

  f1.tracers = (scalar *){q1};
  f2.tracers = (scalar *){q2};
  f3.tracers = (scalar *){q3};
  vof_advection ({f1,f2,f3}, i);

  /**
  We recover the advected velocity field using the total momentum and
  the density */

  foreach()
    foreach_dimension()
      u.x[] = (q1.x[] + q2.x[] + q3.x[])/rho(1.0-f2[]-f3[],f2[],f3[]);
  boundary ((scalar *){u});

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