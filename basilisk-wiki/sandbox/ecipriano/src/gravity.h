/**
# Gravity

We add gravity using a method which is similar to the one used in
[reduced.h](/src/reduced.h), but without the hypotesis of constant
density.

Gravity force in the Navier--Stokes equations can be re-writtes as:
$$
  \rho\mathbf{g}
  = \rho\nabla\left(\mathbf{g}\cdot\mathbf{x}\right)
  = \nabla\left(\rho\mathbf{g}\cdot\mathbf{x}\right)
  - \mathbf{g}\cdot\mathbf{x}\nabla\rho
$$

Therefore, the RHS of the momentum equation can be re-written as:
$$
\partial_t\rho\mathbf{u}+\nabla\cdot(\rho\mathbf{u}\otimes\mathbf{u}) = 
\nabla\cdot(2\mu\mathbf{D})
- \nabla p_d {\color{blue}-\mathbf{g}\cdot\mathbf{x}\nabla\rho}
$$

where $p_d$ is the dynamic pressure, while the blue term is the contribution
added by this module. For a two-phase system with constant
properties, the density is a function of the phase indicator $H$ only.
Therefore:
$$
  \nabla\rho = \left[\rho\right]_\Gamma \nabla H
$$

and this method becomes the implementation of [reduced.h](/src/reduced.h).
If variable properties are considered, the gravitational force is not
reduced at the gas-liquid interface, but it will be non-null in the whole
domain depending on the gradients of the density field.
*/

coord G = {0.,0.,0.}, Z = {0.,0.,0.};

#include "curvature.h"

/**
Interfacial forces are a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

/**
The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  coord G1;
  foreach_dimension()
    G1.x = G.x;

  scalar phig[];
  position (f, phig, G1, Z, add = false);

#if TREE
  for (scalar f in {interfaces}) {
    f.prolongation = p.prolongation;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  scalar rhovar[];
  for (scalar f in {interfaces})
    foreach()
      rhovar[] = rho1v[]*f[] + rho2v[]*(1. - f[]);

#if TREE
  rhovar.prolongation = p.prolongation;
  rhovar.dirty = true;
#endif

  face vector av = a, sth[];
  foreach_face() {
    sth.x[] = 1.;
    if (f[] != f[-1] && fm.x[] > 0.) {
      double phif =
        (phig[] < nodata && phig[-1] < nodata) ?
        (phig[] + phig[-1])/2. :
        phig[] < nodata ? phig[] :
        phig[-1] < nodata ? phig[-1] :
        0.;
      sth.x[] = (phif == 0.) ? 1. : 0.;

      av.x[] -= alpha.x[]/(fm.x[] + SEPS)*phif*(rhovar[] - rhovar[-1])/Delta;
    }
  }

  /**
  Far from interfacial faces, we use the coordinates of the centroids
  of the cells intead of the interface centroids. */

  foreach_face() {
    coord o = {x,y,z};
    double phiof = 0.;
    foreach_dimension()
      phiof += (o.x - Z.x)*G1.x;
    phiof *= sth.x[];

    av.x[] -= alpha.x[]/(fm.x[] + SEPS)*phiof*(rhovar[] - rhovar[-1])/Delta;
  }

#if TREE
  for (scalar f in {interfaces}) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif
}
