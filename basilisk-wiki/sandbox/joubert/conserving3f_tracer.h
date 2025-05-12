/**
# Momentum-conserving advection of velocity for three phase flow
with "VOF tracer".

This file is a modified version of [conserving3f.h](http://www.basilisk.fr/sandbox/joubert/conserving3f.h) and implements
momentum-conserving VOF advection of the velocity
components for the [navier--Stokes](/src/navier-stokes/centered.h) 
compatible with the usage of "VOF tracers" using a [three-phase](three-phase.h) flow formulation.
For now, the vof event originally in [conserving](/src/navier-stokes/conserving.h) have to be placed directly in the code not here.
The interface would be advected two times otherwise.

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
  We do not overload the *vof()* event to transport consistently the volume
  fraction, the momentum and tracers of each phase here but directly
  in the code. This make the code heavier, but seems to work.*/
