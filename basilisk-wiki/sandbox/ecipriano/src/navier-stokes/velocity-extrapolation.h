/**
# Velocity Extrapolation

We want to extrapolate the one-field velocity $\mathbf{u}$, in order
to obtain a liquid velocity $\mathbf{u}_l$ and a gas phase velocity
$\mathbf{u}_g$, which are both continuous and divergence-free at the
interface. The extrapolations are perfomed using Aslam's
PDE based approach. After the extrapolation step, the
[velocity-potential.h](velocity-potential.h) approach can
be used to force the divergence of the velocity to be zero.
This method was proposed by [Palmore et al., 2019](#palmore2019volume).
*/

#include "aslam.h"
#include "redistance.h"

/**
## Field Allocation

We create the fields corresponding to the liquid and gas
phase face velocities. */

face vector ufext1[], ufext2[];
scalar ps1[], ps2[];
mgstats mgpdiv1, mgpdiv2;
extern scalar f;
int nl1 = 1, nl2 = 1;

/**
## Helper Functions

We define a function that converts the vof fraction
to the level set field. */

void vof_to_ls (scalar f, scalar ls) {
  double deltamin = L0/(1 << grid->maxdepth);
  foreach()
    ls[] = -(2.*f[] - 1.)*deltamin*0.75;
#if TREE
  restriction({ls});
#endif
  redistance (ls, imax = 3);
}

/**
We define a function that performs a projection step to
correct the divergence of the velocity field. */

trace
mgstats project_div1 (face vector ufs, scalar ps,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar div[];
  foreach() {
    ps[] = 0.;
    div[] = 0.;
    foreach_dimension()
      div[] += ufext1.x[1] - ufext1.x[];
    div[] /= dt*Delta;
  }

#ifdef PS_IMPLICIT_SOURCE
  scalar marker[];
  mapregion (marker, f, inverse=true, interface=true, nl=2);
  foreach()
    marker[] *= 1./sq(dt);

  mgstats mgp = poisson (ps, div, alpha, lambda=marker,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
#else
  mgstats mgp = poisson (ps, div, alpha,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
#endif

  foreach_face()
    ufs.x[] = -dt*alpha.x[]*face_gradient_x (ps, 0);

  return mgp;
}

trace
mgstats project_div2 (face vector ufs, scalar ps,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar div[];
  foreach() {
    ps[] = 0.;
    div[] = 0.;
    foreach_dimension()
      div[] += ufext2.x[1] - ufext2.x[];
    div[] /= dt*Delta;
  }

#ifdef PS_IMPLICIT_SOURCE
  scalar marker[];
  mapregion (marker, f, inverse=false, interface=true, nl=2);
  foreach()
    marker[] *= 1./sq(dt);

  mgstats mgp = poisson (ps, div, alpha, lambda=marker,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
#else
  mgstats mgp = poisson (ps, div, alpha,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
#endif

  foreach_face()
    ufs.x[] = -dt*alpha.x[]*face_gradient_x (ps, 0);

  return mgp;
}

/**
## Extrapolations

At the end of the solution of the Navier-Stokes equations,
we apply the extrapolation procedure. */

vector uext1[], uext2[];

event end_timestep (i++)
{
  /**
  We need to start from the face one-field velocity $\mathbf{u}_f$
  and to return the two extended velocities on the faces,
  in order to be used for the advection steps. However,
  the extrapolation procedure must be applied on
  cell-centered fields. Therefore, we interpolate from face
  to cell. */

  //vector uext1[], uext2[];
  foreach() {
    foreach_dimension() {
      uext1.x[] = 0.5*(uf.x[1] + uf.x[])*f[];
      uext2.x[] = 0.5*(uf.x[1] + uf.x[])*(1. - f[]);
    }
  }

  /**
  We reconstruct a signed distance field from the volume
  fraction, since the extrapolations need interface normals
  defined over the whole domain of extrapolation. */

  scalar f1[], f2[], ls1[], ls2[];
  foreach() {
    f1[] = (f[] > 1.e-10) ? f[] : 0.;
    f2[] = (1. - f[] > 1.e-10) ? (1. - f[]) : 0.;
  }
  vof_to_ls (f1, ls1);

  foreach()
    ls2[] = -ls1[];

  /**
  We apply the constant extrapolations on the two velocities.
  */

  constant_extrapolation (uext1.x, ls1, 0.5, 10, c=f1, nl=nl1);
  constant_extrapolation (uext1.y, ls1, 0.5, 10, c=f1, nl=nl1);
  constant_extrapolation (uext2.x, ls2, 0.5, 10, c=f2, nl=nl2);
  constant_extrapolation (uext2.y, ls2, 0.5, 10, c=f2, nl=nl2);

  /**
  Finally, we reconstruct the face velocities from the
  extrapolated colocated velocities by linear interpolation.
  */

  foreach_face() {
    ufext1.x[] = 0.5*(uext1.x[] + uext1.x[-1]);
    ufext2.x[] = 0.5*(uext2.x[] + uext2.x[-1]);
  }

  /**
  The extrpolation procedure does not guarantee the final
  velocity to be divergence-free. Therefore, we cancel
  errors on the velocity divergence by solving an additional
  Poisson equation:
  $$
  \nabla \cdot \left( \alpha \nabla \phi \right)
  =
  \nabla \cdot u_f^E \\
  $$
  which is used to correct the extrapolated velocity
  according to:
  $$
  \mathbf{u}_f^E += \nabla \phi
  $$
  */

  face vector ufs1[], ufs2[];

  mgpdiv1 = project_div1 (ufs1, ps1, alpha, dt, mgpdiv1.nrelax);
  mgpdiv2 = project_div2 (ufs2, ps2, alpha, dt, mgpdiv2.nrelax);

  foreach_face() {
    ufext1.x[] -= ufs1.x[];
    ufext2.x[] -= ufs2.x[];
  }
}

/**
## References

~~~bib
@article{palmore2019volume,
  title={A volume of fluid framework for interface-resolved simulations of vaporizing liquid-gas flows},
  author={Palmore Jr, John and Desjardins, Olivier},
  journal={Journal of Computational Physics},
  volume={399},
  pages={108954},
  year={2019},
  publisher={Elsevier}
}
~~~
*/

