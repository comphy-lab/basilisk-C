#define VELOCITY_JUMP 1

#include "common.h"
#include "utils.h"
#include "aslam.h"

scalar jump[];
face vector jumpf[];
vector n[], ur[];
face vector nf[];
scalar pg[], ps[];
bool _boiling = false;

vector u1g[], u2g[];
face vector uf1g[], uf2g[];

extern scalar f;

/**
## Projection Function

We define the function that performs the projection step with the
volume expansion term due to the phase change or due to density
changes. */

#include "poisson.h"

extern scalar * drhodtlist, * plist;
extern face vector * uflist;
mgstats mgpl, mgpe;

trace
mgstats project_jump (face vector * uflist, scalar * plist,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar p1 = plist[1], p2 = plist[0];
  face vector uf1 = uflist[1], uf2 = uflist[0];

  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  scalar drhodt1 = drhodtlist[1], drhodt2 = drhodtlist[0];
  foreach() {
    double div1 = 0., div2 = 0.;
    foreach_dimension() {
      div1 += uf1.x[1] - uf1.x[];
      div2 += uf2.x[1] - uf2.x[];
    }
    div1 /= dt*Delta;
    div2 /= dt*Delta;
    div1 += drhodt1[]/dt;
    div2 += drhodt2[]/dt;
    div[] = (f[] > 0.5) ? div1 : div2;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p1, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach()
    p2[] = p1[];

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face() {
    uf1.x[] -= dt*alpha.x[]*face_gradient_x (p1, 0);
    uf2.x[] -= dt*alpha.x[]*face_gradient_x (p2, 0);
  }

  return mgp;
}


trace
mgstats project_ghost (face vector uf, face vector ufg, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[], drhodt = drhodtlist[1];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
    div[] += drhodt[]/dt;
  }

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
mgstats project_liquid (face vector uf, face vector ufg, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[], drhodt = drhodtlist[1];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += ufg.x[1] - ufg.x[];
    div[] /= dt*Delta;
    div[] += drhodt[]/dt;
  }

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
    uf.x[] = ufg.x[] - dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

#define project(...) project_jump(__VA_ARGS__)
#include "navier-stokes/low-mach.h"
#undef project

event defaults (i = 0) {
  for (int d = 0; d < nboundary; d++) {
    pg.boundary[d] = p.boundary[d];
    pg.boundary_homogeneous[d] = p.boundary_homogeneous[d];
    ps.boundary[d] = p.boundary[d];
    ps.boundary_homogeneous[d] = p.boundary_homogeneous[d];
  }
}

scalar ls[];

event advection_term (i++,last) {
  assert (nv == 2);

  vector u1 = ulist[1], u2 = ulist[0];
  vector uf1 = uflist[1], uf2 = uflist[0];

  gradients ({ls}, {n});
  foreach() {
    double mag = 0.;
    foreach_dimension()
      mag += sq (n.x[]);
    mag = sqrt (mag);
    foreach_dimension()
      n.x[] /= (mag + 1e-10);
  }

  foreach_face()
    nf.x[] = face_value (n.x, 0);

  /**
  We compute the value of the centered ghost velocities. */

  foreach() {
    foreach_dimension() {
      if (_boiling) {
        //u2g.x[] = 0.;
        u1g.x[] = u2.x[] + jump[]*n.x[];
      }
      else {
        //u1g.x[] = 0.;
        u2g.x[] = u1.x[] - jump[]*n.x[];
      }
    }
  }

  /**
  We compute the value of the face ghost velocities. The interface
  normal on the face is computed using the level set. */

  foreach_face() {
    if (_boiling) {
      uf2g.x[] = 0.;
      uf1g.x[] = uf2.x[] + jumpf.x[]*nf.x[]*fm.x[];
    }
    else {
      uf1g.x[] = 0.;
      uf2g.x[] = uf1.x[] - jumpf.x[]*nf.x[]*fm.x[];
    }
  }

  foreach() {
    foreach_dimension() {
      u1.x[] = (f[] > 0.5) ? u1.x[] : u1g.x[];
      u2.x[] = (f[] < 0.5) ? u2.x[] : u2g.x[];
    }
  }

  foreach_face() {
    double ff = face_value (f, 0);
    uf1.x[] = (ff > 0.5) ? uf1.x[] : uf1g.x[];
    uf2.x[] = (ff < 0.5) ? uf2.x[] : uf2g.x[];
  }
}

event viscous_term (i++,last) {
  //setjump (f, ulist[1], ulist[0], jump, n);
  //setjumpf (f, uflist[1], uflist[0], jumpf, nf);
}

event acceleration (i++,last) {
  //setjump (f, ulist[1], ulist[0], jump, n);
  //setjumpf (f, uflist[1], uflist[0], jumpf, nf);
}

event projection (i++,last) {
  if (_boiling) {
    vector uf2 = uflist[0];
    project_ghost (uf2, uf2g, pg, alpha, dt, mgpl.nrelax);

    vector gpg[], u2 = ulist[0];
    centered_gradient (pg, gpg);
    foreach()
      foreach_dimension()
        u2g.x[] = u2.x[] + dt*gpg.x[];
  }
  else {
    vector uf1 = uflist[1];
    project_ghost (uf1, uf1g, pg, alpha, dt, mgpl.nrelax);

    vector gpg[], u1 = ulist[1];
    centered_gradient (pg, gpg);
    foreach()
      foreach_dimension()
        u1g.x[] = u1.x[] + dt*gpg.x[];
  }
}

event end_timestep (i++,last) {
  vector u1 = ulist[1], u2 = ulist[0];

  if (!_boiling) {
    face vector uf1 = uflist[1], uf2 = uflist[0];

    foreach() {
      foreach_dimension() {
        u1g.x[] = u1g.x[];
        u2g.x[] = u1.x[] - jump[]*n.x[];
      }
    }

    foreach_face() {
      uf1g.x[] = uf1g.x[];
      uf2g.x[] = uf1.x[] - jumpf.x[]*nf.x[]*fm.x[];
    }

    face vector ufs[];
    foreach_face() {
      double ff = face_value (f, 0);
      ufs.x[] = (ff > 0.5) ? uf1.x[] : uf1g.x[];
    }
    project_liquid (uf1, ufs, ps, alpha, dt, mgpe.nrelax);

    foreach_face()
      uf1g.x[] = uf1.x[];

    foreach() {
      foreach_dimension() {
        u1.x[] = (f[] > 0.5) ? u1.x[] : u1g.x[];
        u2.x[] = (f[] < 0.5) ? u2.x[] : u2g.x[];
      }
    }

    foreach_face() {
      double ff = face_value (f, 0);
      uf1.x[] = (ff > 0.5) ? uf1.x[] : uf1g.x[];
      uf2.x[] = (ff < 0.5) ? uf2.x[] : uf2g.x[];
    }
  }

  // Reconstruction
  foreach()
    foreach_dimension()
      ur.x[] = u1.x[]*f[] + u2.x[]*(1. - f[]);
}

