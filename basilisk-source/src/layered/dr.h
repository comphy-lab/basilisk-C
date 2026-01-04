/**
# Boussinesq buoyancy

The momentum equation of the [multilayer solver](hydro.h) becomes
$$
\begin{aligned}
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
  {\color{blue} - \mathbf{{\nabla}} (h q)_k + \left[ q 
  \mathbf{{\nabla}} z \right]_k}
\end{aligned}
$$
where the terms in blue have been added and $q$ is the (hydrostatic)
pressure deviation due to (small) density variations.

The density variations are described by a scalar field $T$, which is
advected by the flow, and an associated "equation of state"
$\Delta\rho(T)$ defined by the user. This gives the additional
equations
$$
\begin{aligned}
\partial_t \left( h T \right)_k + \mathbf{{\nabla}} \cdot \left(
h \mathbf{u}  T \right)_k & = 0,\\
q(z) & = \int_0^z g \Delta \rho(T) dz
\end{aligned}
$$
*/

scalar T;

event defaults (i = 0)
{
  T = new scalar[nl];
  tracers = list_append (tracers, T);
}

event cleanup (t = end, last)
{
  delete ({T});
}

/**
The pressure gradient terms in blue are added to the acceleration of
the [multilayer solver](hydro.h). */

event acceleration (i++)
{

  /**
  The hydrostatic pressure deviation $q$ is stored on the interfaces
  between layers (consistently with the Keller box scheme
  discretisation, see [Popinet,
  2020](/Bibliography#popinet2020)). This gives the following vertical
  discrete integration scheme. */

  scalar q = new scalar[nl];
  
  foreach() {
    double ph = 0.;
    for (point.l = nl - 1; point.l >= 0; point.l--) {
      ph += G*drho(T[])*h[];
      q[] = ph;
    }
  }

  /**
  Once the pressure deviation is known, the terms in blue above are
  added to the face acceleration field `ha`, using the
  [pressure-gradient macro](hydro.h#horizontal-pressure-gradient). */
  
  foreach_face()
    hpg (pg, q, 0,
	 ha.x[] += pg);

  delete ({q});
}

/**
## See also

* [Boussinesq buoyancy for isopycnal layers](isopycnal.h)
*/
