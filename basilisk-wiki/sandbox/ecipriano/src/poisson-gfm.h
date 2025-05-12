/**
# Poisson--Helmholtz solver with Ghost Fluid Method

We want to solve a Poisson-–Helmholtz equation:

$$
  \nabla\cdot\left(\alpha\nabla a\right) + \lambda a = b
$$

for a two-phase systems including interface jump conditions using the GFM
approach:

$$
  \left[a\right]_\Gamma = a_\Gamma
$$

$$
  \left[\alpha\nabla a \cdot\mathbf{n}_\Gamma\right]_\Gamma = b_\Gamma
$$

According to [Liu et al. 2000](#liu2000boundary), the final form of the GFM
discretization is similar to the standard 2nd order discretization of the
Poisson equation. However, we have to modify the source terms according to
the interface jump condition.

Therefore, the final discretization of the pressure equation in 1D reads:

$$
\left[\beta_{i+1/2}\left(\dfrac{p_{i+1} - p_{i}}{\Delta}\right)
    - \beta_{i-1/2}\left(\dfrac{p_{i} - p_{i-1}}{\Delta}\right)
\right]\dfrac{1}{\Delta}
=
f_i + \sum_f \dfrac{\beta_f a_{\Gamma, f}}{\Delta^2}
    + \sum_f \dfrac{\beta_f b_{\Gamma, f}\theta_f}{\beta^+ \Delta}
$$

Where the sign of the interface jumps changes depending on the side of the
interface.
*/

#include "poisson.h"

/**
## Field Allocations

We allocate global variables useful to impose th interface boundary conditions:
*/

face vector intface[], intside[], reldist[], absdist[];

/**
## VOF or Level Set formulations

The global face vectors can be computed in different ways, mainly depending on
the approach used: VOF-based or Level Set-based. We declare function pointers
and corresponding functions that control the policy used to compute the global
vectors (VOF by default, but they can be overwritten).
*/

/**
Using VOF, we define interfacial faces according to the following condition:

$$
  \left(f_i - 0.5\right)\left(f_{i-1} - 0.5\right) < 0.
$$
*/

void intface_vof (scalar f) {
  foreach_face()
    intface.x[] = ((f[] - 0.5)*(f[-1] - 0.5) < 0.) ? 1. : 0.;
}

/**
Using VOF, we distinguish the side of the interface as:

$$
  f = \begin{cases}
    > 0.5 & \text{liquid} \\
    < 0.5 & \text{gas}
  \end{cases}
$$
*/

void intside_vof (scalar f) {
  foreach_face()
    intside.x[] = (f[-1] >= 0.5) ? 1. : -1.;
}

/**
Using VOF, we define the relative position of the interface on the current
face as:

$$
  \lambda = \dfrac{f_{i-1} - 0.5}{f_{i-1} - f_{i}}
$$
*/

void reldist_vof (scalar f) {
  foreach_face()
    reldist.x[] = (intface.x[] == 1.) ? (f[-1] - 0.5)/(f[-1] - f[]) : 0.;
}

/**
Using VOF, we define the absolute position of the interface on the current
face as:

$$
  \mathbf{x}_\Gamma = \mathbf{x}[-1] + \lambda \Delta
$$
*/

void absdist_vof (scalar f) {
  vector xp[];
  foreach() {
    coord o = {x,y,z};
    foreach_dimension()
      xp.x[] = o.x;
  }

  foreach_face()
    absdist.x[] = xp.x[-1] + reldist.x[]*Delta;
}

/**
Similar functions can be defined for a Level Set approach. */

void intface_ls (scalar d) {
  foreach_face()
    intface.x[] = (d[]*d[-1] < 0.) ? 1. : 0.;
}

void intside_ls (scalar d) {
  foreach_face()
    intside.x[] = (d[-1] > 0.) ? 1. : -1.;
}

void reldist_ls (scalar d) {
  foreach_face()
    reldist.x[] = (intface.x[] == 1.) ? fabs(d[-1])/(fabs(d[]) + fabs(d[-1])) : 0.;
}

void absdist_ls (scalar d) {
  absdist_vof (d);
}

/**
We link the function pointers to their default functions. */

void (* intfacefun)(scalar) = intface_vof;
void (* intsidefun)(scalar) = intside_vof;
void (* reldistfun)(scalar) = reldist_vof;
void (* absdistfun)(scalar) = absdist_vof;

/**
## *poisson_ghost()*: Extension of the *poisson()* function including the GFM

The same arguments used by *poisson()* must be provided, together with the
following additional inputs:

* *ajump*: jump in the variable being solved: $\left[a\right]_\Gamma$
* *bjump*: jump in the derivative of the variable being solved: $\left[\beta\nabla a\cdot \mathbf{n}_\Gamma\right]_\Gamma / \beta_k$
*/

mgstats poisson_ghost (scalar a, scalar b,
     (const) scalar f = {-1},
     (const) face vector ajump = {{-1}},
     (const) face vector bjump = {{-1}},
     (const) face vector alpha = {{-1}},
     (const) scalar lambda = {-1},
     double tolerance = 0.,
     int nrelax = 4,
     int minlevel = 0,
     scalar * res = NULL,
     double (* flux) (Point, scalar, vector, double *) = NULL)
{
  /**
  We unpack the optional inputs for convenience, setting the jump conditions to
  zero by default. */

  if (f.i < 0)
    f = zeroc;
  if (ajump.x.i < 0)
    ajump = zerof;
  if (bjump.x.i < 0)
    bjump = zerof;
  if (alpha.x.i < 0)
    alpha = unityf;
  if (lambda.i < 0) {
    const scalar zeroc[] = 0.; // fixme
    lambda = zeroc;
  }

  /**
  We compute the face vectors with the interfacial faces, the side of the
  interface, the relative and absolute distances of the interface. */

  intfacefun (f);
  intsidefun (f);
  reldistfun (f);
  absdistfun (f);

  /**
  The know term `b` is modified according to the jump conditions
  discretization. */

  face vector aj[];
  foreach_face() {
    if (intface.x[]) {
      if (intside.x[] > 0.)
        aj.x[] = alpha.x[]*ajump.x[]/sq(Delta);
      else
        aj.x[] = -alpha.x[]*ajump.x[]/sq(Delta);
    }
    else
      aj.x[] = 0.;
  }

  face vector bj[];
  foreach_face() {
    if (intface.x[]) {
      if (intside.x[] > 0.)
        bj.x[] = alpha.x[]*bjump.x[]*reldist.x[]/Delta;
      else
        bj.x[] = alpha.x[]*bjump.x[]*reldist.x[]/Delta;
    }
    else
      bj.x[] = 0.;
  }

  foreach() {
    double js = 0.;
    foreach_dimension()
      js -= aj.x[1] - aj.x[];
    foreach_dimension()
      js += bj.x[1] + bj.x[];
    b[] += js;
  }

  /**
  We call the multigrid poisson solver. */

  mgstats mgp = poisson (a, b, alpha, lambda, tolerance,
      nrelax, minlevel, res, flux);

  return mgp;
}

/**
## Notes and Improvements

1. The derivative jump `bjump` must be provided by the user already divided by
$\beta^+$ or $\beta^-$ for variable coefficient $\beta$ simulations, depending
on the side of the interface.

2. The mechanism to include the derivative jump should be reviewed. However, it
is not necessary for applications to surface tension.


## References

~~~bib
@article{vukvcevic2017implementation,
  title={Implementation of the ghost fluid method for free surface flows in polyhedral finite volume framework},
  author={Vuk{\v{c}}evi{\'c}, Vuko and Jasak, Hrvoje and Gatin, Inno},
  journal={Computers \& fluids},
  volume={153},
  pages={1--19},
  year={2017},
  publisher={Elsevier}
}

@article{liu2000boundary,
  title={A boundary condition capturing method for Poisson's equation on irregular domains},
  author={Liu, Xu-Dong and Fedkiw, Ronald P and Kang, Myungjoo},
  journal={Journal of computational Physics},
  volume={160},
  number={1},
  pages={151--178},
  year={2000},
  publisher={Elsevier}
}
~~~
*/
