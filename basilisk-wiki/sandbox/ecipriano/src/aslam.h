/**
# Aslam Extrapolations

This module implements the PDE-based Aslam extrapolation
methods ([Aslam 2003](#aslam2004partial)). Using these
functions, a discontinuous field can be extended from the
liquid phase to the gas-phase and vice-versa using
constant and extrapolations.

These methods are expensive, since the PDE must be solved
at steady-state. However, most of the times they are
used just for a few steps in order to extend the field for
a few cells across the interface, without covering the whole
domain.
*/

#include "mapregion.h"
#include "redistance.h"

void vof_to_ls (scalar f, scalar levelset, int imax = 3) {
  double deltamin = L0/(1 << grid->maxdepth);
  foreach()
    levelset[] = -(2.*f[] - 1.)*deltamin*0.75;
#if TREE
  restriction({levelset});
#endif
  redistance (levelset, imax);
}

/**
## Constant Extrapolation

The constant extrapolation of a field *f*, existing only in
a portion of space defined by *ls*, can be extrapolated in
the whole domain solving the following PDE:

$$
\dfrac{\partial f}{\partial t}
+ H(\phi)\hat{\mathbf{n}}\cdot\nabla f = 0
$$

where $H(\phi)$ is the Heaviside function, used to define the
region where the field should be extrapolated and where it
should not be modified:

$$
H(\phi) =
\begin{cases}
  1 & \text{if } \phi > 0,\\
  0 & \text{if } \phi \leq 0.
\end{cases}
$$

while $\hat{\mathbf{n}}$ is the unit normal.
*/

trace
void constant_extrapolation (
  scalar f,                   // field to extrapolate
  scalar ls,                  // level set field
  double cfl = 0.5,           // CFL number
  int nmax,                   // number of maximum time steps
  (const) scalar s = {-1},    // source term, default zero
  (const) scalar c = {-1},    // vof field, optional
  int nl = 0,                 // from which layer of cells (optional, min=0, max=2)
  int nointerface = 0         // remove the interface from the heaviside (default false)
)
{
  scalar H[];
  vector n[], gf[];

  if (s.i < 0)
    s = zeroc;

  /**
  We compute the gradients of the level set for
  the calcation of the interface normals.
  */

  gradients ({ls}, {n});

  /**
  We compute the normals and the heaviside function H,
  which is non-null in the region where the field must
  be extrapolated. In case of vof field, the user can
  decide to extrapolate the field from a layer of cells
  which in not adjacent to the interface. This can be
  specified setting the variables *nl*, 0 by default,
  it can be set to 1 or 2. */

  if (c.i > 0)
    mapregion (H, c, nl=nl, nointerface=nointerface);
  else
    foreach()
      H[] = (ls[] <= 0.) ? 0. : 1.;

  foreach() {
    double maggf = 0.;
    foreach_dimension()
      maggf += sq (n.x[]);
    maggf = sqrt (maggf);
    foreach_dimension()
      n.x[] /= (maggf + 1.e-10);
  }

  /**
  We solve the PDE inside a loop over the maximum number
  of time steps that, in principle, should ensure to arrive
  at steady-state. In practice, less time-steps will be
  sufficients to extrapolate the fields over the interface. */

  int ts = 0;

  while (ts < nmax) {

    /**
    The gradients of the extrapolated functions are updated
    using an upwind scheme. */

    foreach()
      foreach_dimension()
        gf.x[] = (n.x[] <= 0.) ? (f[1] - f[])/Delta : (f[] - f[-1])/Delta;

    /**
    We solve a single step of the PDE. */

    foreach() {
      double dt = cfl*Delta;
      double nscalargf = 0.;
      foreach_dimension()
        nscalargf += n.x[]*gf.x[];
      f[] -= dt*H[]*(nscalargf - s[]);
    }
    ts++;
  }
}

/**
## Linear Extrapolation

The linear extrapolation of the field *f* along the normal
direction, can be achieved from the solution of the
following PDE:

$$
\dfrac{\partial f}{\partial t}
+ H(\phi)\left(\hat{\mathbf{n}}\cdot\nabla f - f_n \right) = 0
$$

Which is similar to the equation resolved with the constant
extrapolation procedure, except for the source term $f_n$,
which is the directional derivative of $f$ in the normal
direction, defined as:

$$
f_n = \hat{\mathbf{n}}\cdot \nabla f
$$

Since this function is defined only in the region where $f$
is defined, we apply the constant extrapolation in order to
obtain a field $f_n$ defined on the whole domain:

$$
\dfrac{\partial f_n}{\partial t}
+ H(\phi)\hat{\mathbf{n}}\cdot\nabla f_n = 0
$$

Therefore, the linear interpolation requires the solution of
two different extrapolation PDEs, and the same logic applies
if higher order extrapolations have to be solved.
*/

void linear_extrapolation (
  scalar f,                   // field to extrapolate
  scalar ls,                  // level set field
  double cfl,                 // CFL number
  int nmax,                   // number of maximum time steps
  (const) scalar s = {-1},    // source term, default zero
  (const) scalar c = {-1},    // vof field, optional
  int nl = 0,                 // from which layer of cells (optional, min=0, max=2)
  int nointerface = 0         // remove the interface from the heaviside (default false)
)
{
  scalar H[], fn[];
  vector n[], gf[];

  if (s.i < 0)
    s = zeroc;

  /**
  We compute the gradients of the level set for
  the calcation of the interface normals, and the
  gradient of *f* for the calculation of the directional
  derivative. */

  gradients ({ls, f}, {n, gf});

  /**
  We compute the normals and the heaviside function H,
  which is non-null in the region where the field must
  be extrapolated. In case of vof field, the user can
  decide to extrapolate the field from a layer of cells
  which in not adjacent to the interface. This can be
  specified setting the variables *nl*, 0 by default,
  it can be set to 1 or 2. */

  if (c.i > 0)
    mapregion (H, c, nl = (nl == 0.) ? 1. : min (nl+1, 2),
        nointerface=nointerface);
  else
    foreach()
      H[] = (ls[]+Delta <= 0.) ? 0. : 1.;

  foreach() {
    double maggf = 0.;
    foreach_dimension()
      maggf += sq (n.x[]);
    maggf = sqrt (maggf);
    foreach_dimension()
      n.x[] /= (maggf + 1.e-10);
  }

  /**
  We compute the directional derviative *fn*. */

  foreach()
    foreach_dimension()
      fn[] += n.x[]*gf.x[];

  /**
  We solve the constant extrapolation for extending
  the directional derivative. */

  int ts = 0;

  while (ts < nmax) {

    /**
    The gradients of the extrapolated functions are updated
    using an upwind scheme. */

    foreach()
      foreach_dimension()
        gf.x[] = (n.x[] <= 0.) ? (fn[1] - fn[])/Delta : (fn[] - fn[-1])/Delta;

    /**
    We solve a single step of the PDE. */

    foreach() {
      double dt = cfl*Delta;
      double nscalargf = 0.;
      foreach_dimension()
        nscalargf += n.x[]*gf.x[];
      fn[] -= dt*H[]*(nscalargf - s[]);
    }
    ts++;
  }

  /**
  We solve a constant extrapolation with source equal to
  the directional derivative. */

  if (c.i > 0)
    constant_extrapolation (f, ls, cfl, nmax, fn, c, nl, nointerface);
  else
    constant_extrapolation (f, ls, cfl, nmax, fn);
}

/**
## References

~~~bib
@article{aslam2004partial,
  title={A partial differential equation approach to multidimensional extrapolation},
  author={Aslam, Tariq D},
  journal={Journal of Computational Physics},
  volume={193},
  number={1},
  pages={349--355},
  year={2004},
  publisher={Elsevier}
}
~~~
*/

