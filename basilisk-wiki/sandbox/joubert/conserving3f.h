/**
# Momentum-conserving advection of velocity for three phase flow.

This file is a modified version of [conserving](/src/navier-stokes/conserving.h) and implements momentum-conserving VOF advection of the velocity components for the [navier--Stokes](/src/navier-stokes/centered.h) using a [three-phase](three-phase.h) flow formulation. 

On trees, we first define refinement and restriction functions which
guarantee conservation of each component of the total momentum. Note
that these functions do not guarantee conservation of momentum for
each phase.

We use the *interfaces* list to iterates through all the volume
fractions.*/

extern scalar * interfaces;

#if TREE
static void momentum_refine (Point point, scalar u) {
  refine_bilinear (point, u);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhom(f1[],f2[],f3[])*u[];
  double du = u[] - rhou/((1 << dimension)*cm[]*rhom(f1[],f2[],f3[]));
  foreach_child()
    u[] += du;
}

static void momentum_restriction (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rhom(f1[],f2[],f3[])*u[];
  u[] = rhou/((1 << dimension)*cm[]*rhom(f1[],f2[],f3[]));
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
    components of velocity. To ensure this in a three-phase flow we move all
    the volume fractions in the *interfaces* list to the front of
    the field list (*all*). */

  int i = 0;
  for (scalar c in interfaces) {
    while (all[i].i != c.i) i++;
    while (i > 0 && all[i].i)
      all[i] = all[i-1], i--;
    all[i] = c;
  }

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
  We need to overload the stability event so that the CFL is taken into
  account (because we set stokes to true). */
  
  event stability (i++)
  dtmax = timestep (uf, dtmax);

/**
  We will transport the three components of the momentum, $q_2=f2 \rho_2
  \mathbf{u}$ and $q_3=f3 \rho_3 \mathbf{u}$ and $q_1=f1 \rho_1 \mathbf{u}$  We
  will need to impose boundary conditions which match this definition. This is
  done using the functions below. */


  foreach_dimension()
static double boundary_q2_x (Point neighbor, Point point, scalar q2, void * data)
{
  return clamp(f2[],0.,1.)*rho2*u.x[];
}

  foreach_dimension()
static double boundary_q3_x (Point neighbor, Point point, scalar q3, void * data)
{
  return clamp(f3[],0.,1.)*rho3*u.x[];
}

  foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q1, void * data)
{
  return  clamp(f1[],0.,1.)*rho1*u.x[];
}

/**
  Similarly, on trees we need prolongation functions which also follow
  this definition. */

#if TREE
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

foreach_dimension()
  static void prolongation_q1_x (Point point, scalar q1) {
    foreach_child()
      q1[] = clamp(f1[],0.,1.)*rho1*u.x[];
  }
#endif

/**
  We overload the *vof()* event to transport consistently the volume
  fraction and the momentum of each phase. */

static scalar * interfaces1 = NULL;

event vof (i++) {

  /**
    We allocate three temporary vector fields to store the three components
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
    We split the total momentum $q$ into its three components $q2$,$q3$ and $q1$
    associated with $f2$, $f3$ and $f1$ respectively. */

  foreach()
    foreach_dimension() {
      double fc1 = clamp(f1[],0,1);
      double fc2 = clamp(f2[],0,1);
      double fc3 = clamp(f3[],0,1);
      q2.x[] = fc2*rho2*u.x[];
      q3.x[] = fc3*rho3*u.x[];
      q1.x[] = fc1*rho1*u.x[];

    }
  boundary ((scalar *){q1,q2,q3});

  /**
    We use the same slope-limiting as for the
    velocity field. */

  foreach_dimension() {
    q1.x.gradient = q2.x.gradient = q3.x.gradient = u.x.gradient;
  }

  /**
    We associate the transport of $q1$ and $q2$ with $f$ and transport
    all fields consistently using the VOF scheme. */
  f1.tracers = (scalar *){q1};
  f2.tracers = (scalar *){q2};
  f3.tracers = (scalar *){q3};
  vof_advection ({f1, f2, f3}, i);

  /**
    We recover the advected velocity field using the total momentum and
    the density */

  foreach()
    foreach_dimension()
    u.x[] = (q1.x[] + q2.x[]+ q3.x[])/rhom(f1[],f2[],f3[]);
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
