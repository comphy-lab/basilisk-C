/**
# Eulerian elasticity

In this file we define and advect the tensors that track the elongation of a
region of space $\Gamma$ occupied by an elastic solid. It is meant to be
combined with [capsule.h](capsule.h) to provide a description of biological
membranes.

It was shown by [Barthès-Biesel and Rallison](#berthes1981time) that the stress
tensor of an elastic membrane can be described in an Eulerian fashion with the
left Cauchy-Green deformation tensor $\bm{B}$ and the local projector onto the
membrane $\mathbf{P} = \mathbf{I} - \mathbf{n}\mathbf{n}$, with
$\mathbf{n}$ denoting the interface unit normal vector.

Recently, it was demonstrated by [Ii et al.](#ii2012full) that this
Eulerian expression of the stress tensor can be successfully implemented in order to
simulate elastic capsules without meshing their membrane and without requiring a
boundary-fitted fluid mesh or an immersed boundary method (IBM). They report
that it is more stable to advect the *modified* left Cauchy-Green deformation
tensor $\mathbf{G} = J^{-1}\mathbf{B}$, where
$J = \sqrt{(tr(\mathbf{B})^2 - tr(\mathbf{B}^2))/2}$ is the surface Jacobian
of the membrane, i.e. the local membrane area change [](#ii2012computational).

Both $J$ and $\mathbf{G}$ are only defined in the capsule region $\Gamma$, and
in this region they follow the advection equations:

$$ \partial_t J + \mathbf{u} \cdot \nabla J = (\nabla_s \cdot \mathbf{u}) J $$
$$ \partial_t \mathbf{G} + \mathbf{u} \cdot \nabla \mathbf{G} =
\mathbf{G} \cdot \nabla_s \mathbf{u} +  (\nabla_s \mathbf{u})^T \cdot \mathbf{G}
- (\nabla_s \cdot \mathbf{u}) \mathbf{G}
$$
$$\text{with} \quad \nabla_s = \mathbf{P} \cdot \nabla$$
denoting by $\mathbf{n}$ the interface normal, $\nabla_s$ the surface gradient and
$\nabla_s \cdot$ the surface divergence.
*/

/**
We have the option of defining a pre-inflated capsule by setting $J_0$ to a
value greater than 1. By default, there is no stress on the capsule in the
reference configuration.
*/
#ifndef J0
  #define J0 1.
#endif

// /**
// We can also advect the Eulerian quantities in the sole region of the membrane
// (instead of going through all grid cells and computing zero fluxes when we are
// out of the membrane area). In that case, we only need one layer of initialized
// cells close the the membrane
// */
// #ifndef LOCAL_BCG
//   #define LOCAL_BCG 0
// #endif
//
// #if (LOCAL_BCG)
//   #include "local_bcg.h"
// #endif

#ifndef EXTEND_MB_ATTRIBUTES
  #define EXTEND_MB_ATTRIBUTES 1
#endif
#ifndef SWITCH_GRADU_CONVENTION
  #define SWITCH_GRADU_CONVENTION 0
#endif
#ifndef SWITCH_SGRADU_TRANSPOSE
  #define SWITCH_SGRADU_TRANSPOSE 0
#endif


/**
We now declare the centered Eulerian quantities: the modified left Cauchy-Green
deformation tensor $\mathbf{G}$ and its source term in its advection equation
$\mathbf{\text{Source}_G}$, the velocity gradient and surface gradient,
a vector defining normals to the interface in the whole membrane region $\Gamma$,
the surface Jacobian $J$ and its source term in its advection equation
$\text{Source}_J$.
*/

scalar J[], sJ[];
symmetric tensor G[], sG[];

event defaults (i = 0) {
  for (scalar s in {G, sG, J, sJ}) {
      s.v.x.i = -1;
      foreach_dimension() {
        s[left] = neumann(0.);
        s[right] = neumann(0.);
      }
  }
  foreach() {
    J[] = J0;
    foreach_dimension() {
      G.x.x[] = 1.;
    }
    G.x.y[] = 0.;
    G.x.z[] = 0.;
    G.y.z[] = 0.;
  }
}

event init (i = 0) {
  boundary((scalar *) {J, sJ, G, sG});
}

event tracer_advection (i++) {
  foreach() {
    /**
    We loop through the membrane cells, and prepare construction of RHS of
    the advection equations of $J$ and $\bm{G}$:
    */
    sJ[] = 0.;
    foreach_dimension() {
      sG.x.x[] = 0.;
    }
    sG.x.y[] = 0.;
    sG.x.z[] = 0.;
    sG.y.z[] = 0.;

    #if EXTEND_MB_ATTRIBUTES
      if (IS_INTERFACE_CELL(point,f)) {
    #else
      if (GAMMA) {
    #endif
      /**
      We construct the velocity gradient tensor at the center of the cells.
      */
      pseudo_t grad_u, sgrad_u;
      foreach_dimension() {
        grad_u.x.x = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);
        grad_u.x.y = (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta);
        grad_u.x.z = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      }

      foreach_dimension() {
      /**
      We then construct the velocity surface gradient and surface divergence
      at the center of the cells.
      */
        sgrad_u.x.x = (1 - sq(extended_n.x[]))*grad_u.x.x
                        - extended_n.x[]*extended_n.y[]*grad_u.x.y
                        - extended_n.x[]*extended_n.z[]*grad_u.x.z;
        sgrad_u.x.y = (1 - sq(extended_n.y[]))*grad_u.x.y
                        - extended_n.x[]*extended_n.y[]*grad_u.x.x
                        - extended_n.y[]*extended_n.z[]*grad_u.x.z;
        sgrad_u.x.z = (1 - sq(extended_n.z[]))*grad_u.x.z
                        - extended_n.x[]*extended_n.z[]*grad_u.x.x
                        - extended_n.y[]*extended_n.z[]*grad_u.x.y;
        sJ[] += sgrad_u.x.x;
      }
      /** At this point $sJ$ contains the surface divergence of u, which
      we use to compute $sG$*/
      foreach_dimension() {
        sG.x.x[] = 2*(sgrad_u.x.x*G.x.x[] + sgrad_u.x.y*G.y.x[]
          + sgrad_u.x.z*G.z.x[]) - sJ[]*G.x.x[];
      }
      sG.x.y[] = sgrad_u.x.x*G.x.y[] + sgrad_u.x.y*G.y.y[]
                + sgrad_u.x.z*G.z.y[] + G.x.x[]*sgrad_u.y.x
                + G.x.y[]*sgrad_u.y.y + G.x.z[]*sgrad_u.y.z - sJ[]*G.x.y[];
      sG.x.z[] = sgrad_u.x.x*G.x.z[] + sgrad_u.x.y*G.y.z[]
                + sgrad_u.x.z*G.z.z[] + G.x.x[]*sgrad_u.z.x
                + G.x.y[]*sgrad_u.z.y + G.x.z[]*sgrad_u.z.z - sJ[]*G.x.z[];
      sG.y.z[] = sgrad_u.y.x*G.x.z[] + sgrad_u.y.y*G.y.z[]
                + sgrad_u.y.z*G.z.z[] + G.y.x[]*sgrad_u.z.x
                + G.y.y[]*sgrad_u.z.y + G.y.z[]*sgrad_u.z.z - sJ[]*G.y.z[];
      sJ[] *= J[];
    } // end if the cell is in the membrane region
  } // foreach

  /** We now use the Bell-Colella-Glaz scheme to compute the material
  derivative */
  boundary((scalar *){J, G});
  advection((scalar *) {J, G}, uf, dt);

  foreach() {
    #if EXTEND_MB_ATTRIBUTES
      if (IS_INTERFACE_CELL(point,f)) {
    #else
      if (GAMMA) {
    #endif
      /** We advance in time $J$ and \mathbf{G} using their respective
      source terms*/
      J[] += dt*sJ[];
      foreach_dimension() {
        G.x.x[] += dt*sG.x.x[];
      }
      G.x.y[] += dt*sG.x.y[];
      G.x.z[] += dt*sG.x.z[];
      G.y.z[] += dt*sG.y.z[];
    }
  }
  boundary((scalar *){J, G});

  /**Then we extend $J$ and \mathbf{G} in the normal direction from the interface
  in order to force these quantities to be constant on the interface. In this
  process, we re-use $sJ$ and $sG$ as temporary field values.*/
  #if EXTEND_MB_ATTRIBUTES
    normal_scalar_extension({J, G}, (scalar *) {sJ, sG}, nb_iter_extension);
  #endif

  foreach() {
    if (GAMMA) {
      /**Finally, we enforce $\bm{n}\cdot\bm{G} = \bm{G}\cdot\bm{n} = 0$ by replacing $\bm{G}$ by $\bm{P}\cdot\bm{G}\cdot\bm{P}$ */
      pseudo_t P, GP;
      foreach_dimension() {
        P.x.x = 1 - sq(extended_n.x[]);
        P.x.y = -extended_n.x[]*extended_n.y[];
        P.x.z = -extended_n.x[]*extended_n.z[];
      }
      foreach_dimension() {
        GP.x.x = G.x.x[]*P.x.x + G.x.y[]*P.y.x + G.x.z[]*P.z.x;
        GP.x.y = G.x.x[]*P.x.y + G.x.y[]*P.y.y + G.x.z[]*P.z.y;
        GP.x.z = G.x.x[]*P.x.z + G.x.y[]*P.y.z + G.x.z[]*P.z.z;
      }
      foreach_dimension() {
        G.x.x[] = P.x.x*GP.x.x + P.x.y*GP.y.x + P.x.z*GP.z.x;
      }
      G.x.y[] = P.x.x*GP.x.y + P.x.y*GP.y.y + P.x.z*GP.z.y;
      G.x.z[] = P.x.x*GP.x.z + P.x.y*GP.y.z + P.x.z*GP.z.z;
      G.y.z[] = P.y.x*GP.x.z + P.y.y*GP.y.z + P.y.z*GP.z.z;
    }
    else {
      J[] = J0;
      foreach_dimension() {
        G.x.x[] = 1.;
      }
      G.x.y[] = 0.;
      G.x.z[] = 0.;
      G.y.z[] = 0.;
    }
  }
  boundary((scalar *){J, G});
} // tracer_advection
