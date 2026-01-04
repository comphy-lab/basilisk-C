/**
# Axisymmetric stream function

The axisymmetric or [Stokes stream
function](https://en.wikipedia.org/wiki/Stokes_stream_function) $\psi$
verifies
$$
u_r = -\frac{1}{r}\partial_z\psi
$$
$$
u_z = +\frac{1}{r}\partial_r\psi
$$
with $u_r$ and $u_z$ the radial and longitudinal velocity components
of an incompressible axisymmetric flow.

Given the definition of the vorticity $\omega$
$$
\omega = \partial_zu_r - \partial_ru_z
$$
$\psi$ verifies the variable-coefficient Poisson equation
$$
\partial_z\left(\frac{1}{r}\partial_z\psi\right) +
\partial_r\left(\frac{1}{r}\partial_r\psi\right) = - \omega
$$
This is not well-conditioned due to the divergence of the coefficients
at $r = 0$ and can be rewritten instead as
$$
\frac{\partial^2\psi}{\partial z^2} + 
\frac{\partial^2\psi}{\partial r^2} - \frac{1}{r}\partial_r\psi 
= - \omega r
$$
which is better behaved.

This equation can be inverted using the multigrid solver combined with
the relaxation and residual functions defined below. */

#include "poisson.h"

static void relax_psi (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  foreach_level_or_leaf (l)
    a[] = (- sq(Delta)*b[] - Delta*(a[0,1] - a[0,-1])/(2.*y) +
	   a[1] + a[-1] + a[0,1] + a[0,-1])/4.;
}

static double residual_psi (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = face_gradient_x (a, 0);
  foreach (reduction(max:maxres)) {
    res[] = b[] + (a[0,1] - a[0,-1])/(2.*y*Delta);
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + (a[0,1] - a[0,-1])/(2.*y*Delta);
    foreach_dimension()
      res[] += (face_gradient_x (a, 0) -
		face_gradient_x (a, 1))/Delta;  
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  return maxres;
}

/**
The function *axistream()* takes as input the $u$ velocity field (with
$x=z$ and $y=r$) and returns the corresponding axisymmetric stream
function. */

mgstats axistream (vector u, scalar psi)
{
  scalar omega[];
  foreach() {
    omega[] = y*(u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    psi[] = 0.;
  }
  psi[bottom] = dirichlet(0);
  return mg_solve ({psi}, {omega}, residual_psi, relax_psi, NULL, 0, NULL, 1);
}

/**
## See also

* [Axisymmetric metric](axi.h)
*/
