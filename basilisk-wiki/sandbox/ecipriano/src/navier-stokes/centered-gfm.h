/**
# Incompressible Navier--Stokes solver with Ghost Fluid Method for surface tension (centered formulation)

This extension of the [centered.h](/src/navier-stokes/centered.h)
Navier--Stokes equations solver implements the GFM according to the VOF-based
model proposed by [Vukvcevic et al. 2017](#vukvcevic2017implementation).

We solve the incompressible Navier--Stokes equations:

$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The difference with respect to the classic centered formulation is that the
surface tension is not added as an acceleration term in the momentum equation,
but it is directly included in the discretization of the Poisson equation
at the projection step.
*/

/**
## Ghost Fluid Method Properties

We declare a structure containing properties which are useful for the GFM, such
as the surface tension coefficient `sigma`, the volume fraction `f`, and
possibly gravity `G`, which can be accounted for in the pressure jump of the
GFM formulation. The gravity `G` must be provided multiplied by the density jump,
due to the order of inclusions of the headers:
$$
  \mathbf{G} = \mathbf{g}[\rho] = \mathbf{g}(\rho_1 - \rho_2)
$$
A constant curvature can be imposed simply providing a value for the attribute
`kappa` of the following structure.
*/

struct GhostFluidMethod {
  scalar f;
  double sigma;
  double kappa;
  coord G;
};

struct GhostFluidMethod gfm = {
  .sigma = 0.,
  .kappa = 0.,
  .G = {0.,0.,0.},
};

/**
## Curvature Calculation

We compute the curvature using the height--functions, just like in
[curvature.h](/src/curvature.h), but we need to overwrite the function
`interfacial()`, because we want to compute the curvature on a wider stencil,
according to the different definition of interfacial cell used here ad adapted
from [Vukvcevic et al. 2017](#vukvcevic2017implementation).
*/

#include "fractions.h"
#include "curvature.h"
#include "poisson-gfm.h"

scalar kappag[];

static inline int interfacial_ghost (Point point, scalar c) {
  int res = 0;
  foreach_dimension() {
    if (intface.x[1] || intface.x[])
      res++;
  }
  return res;
}

trace
cstats curvature_ghost (scalar c, scalar kappa,
      double sigma = 1.[0], bool add = false)
{
  int sh = 0, sf = 0, sa = 0, sc = 0;

  /**
  On trees we set the prolongation and restriction functions for
  the curvature. */
  
#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

#if dimension > 1
  
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);

  foreach(reduction(+:sh) reduction(+:sf)) {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */

    if (!interfacial_ghost (point, c))
      k[] = nodata;

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_curvature (point, c, h)) != nodata)
      sh++;
    else if ((k[] = height_curvature_fit (point, c, h)) != nodata)
      sf++;
  }
  
  foreach (reduction(+:sa) reduction(+:sc)) {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial_ghost (point, c)) {

      /**
      ...or the average of the curvatures in the $3^{d}$ neighborhood
      of interfacial cells. */

      double sk = 0., a = 0.;
      foreach_neighbor(1)
        if (k[] < nodata)
          sk += k[], a++;
        if (a > 0.)
          kf = sk/a, sa++;
        else

  /**
  Empty neighborhood: we try centroids as a last resort. */

  kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = nodata;

    /**
    We add or set *kappa*. */

    if (kf == nodata)
      kappa[] = nodata;
    else if (add)
      kappa[] += sigma*kf;
    else
      kappa[] = sigma*kf;
  }

#else // dimension == 1
  foreach() {
    if (!interfacial_ghost (point, c))
      kappa[] = nodata;
    else {
      double r = x + sign(c[-1] - c[1])*(clamp(c[],0.,1.) - 0.5)*Delta;
      double p = r > 0. ? - 2.*sigma/r : 0.;
      if (add)
        kappa[] += p;
      else
        kappa[] = p;
    }
  }
#endif // dimension == 1

  return (cstats){sh, sf, sa, sc};
}

/**
## Ghost Fluid Method Discretization

According to [Liu et al. 2000](#liu2000boundary), the final form of the GFM
discretization is similar to the standard 2nd order discretization of the
Poisson equation. However, we have to modify the source terms according to
the pressure jump at the interface:

$$
\left[p\right]_\Gamma = \sigma\kappa
$$

We consider that the jump of pressure gradient at the interface is null:

$$
\left[\partial_n p\right] = 0
$$

Therefore, the final discretization of the pressure equation in 1D reads:

$$
\left[\beta_{i+1/2}\left(\dfrac{p_{i+1} - p_{i}}{\Delta}\right)
    - \beta_{i-1/2}\left(\dfrac{p_{i} - p_{i-1}}{\Delta}\right)
\right]\dfrac{1}{\Delta}
=
f_i + \dfrac{\beta_{i-1/2}\left[p\right]_\Gamma}{\Delta^2}
    - \dfrac{\beta_{i+1/2}\left[p\right]_\Gamma}{\Delta^2}
$$

Where the pressure jump $\left[p\right]_\Gamma$ changes sign according to the
side of the interface.
*/


attribute {
  double sigma;
}

vector ubackup[];
face vector pjump[];
face vector betac[];

trace
mgstats project_ghost (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  /**
  We recover the volume fraction field. */

  scalar f = gfm.f;

  /**
  We store the current centered velocity, this is required to avoid to make
  modifications to the centered solver. */

  extern vector u;
  foreach()
    foreach_dimension()
      ubackup.x[] = u.x[];

  /**
  We compute interface properties. */

  intfacefun (f);
  intsidefun (f);
  reldistfun (f);
  absdistfun (f);

  /**
  We compute the curvature field for the GFM: `kappag`. */

  curvature_ghost (f, kappag, sigma = gfm.sigma);

  if (gfm.kappa)
    foreach()
      kappag[] = gfm.kappa*gfm.sigma;

  /**
  And we interpolate the curvature on the faces based on the relative position
  of the interface. Other interpolation methods can also be tested. */

  face vector kappaf[];
  foreach_face() {
    if (intface.x[])
      kappaf.x[] = reldist.x[]*kappag[] + (1. - reldist.x[])*kappag[-1];
    else
      kappaf.x[] = 0.;
  }

  /**
  We compute the pressure jump on the faces. */

  foreach_face()
    pjump.x[] = kappaf.x[];

  /**
  We add the gravity contribution in the pressure jump. The term `gfm.G.x` must
  be multiplied by the density jump at the interface outside this module, in
  order to avoid problems with the variables `rho1` and `rho2` which are
  defined in [two-phase.h](/src/two-phase.h) and not always included. */

  foreach_face() {
    double Gdist = 0.;
    foreach_dimension()
      Gdist += gfm.G.x * absdist.x[];
    pjump.x[] -= fm.x[]*Gdist;
  }

  /**
  The diffusion coefficient $\beta=1/\rho$ is computed like the field `alpha`
  using an arithmetic average based on the volume fraction. This is different
  with respect to the approach described by [Vukvcevic et al. 2017](#vukvcevic2017implementation)
  but it seems to work better here. */

  foreach_face()
    betac.x[] = alpha.x[];

  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson_ghost (p, div, f = f, ajump = pjump, alpha = betac,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  /**
  We recover the face velocity values from the pressure gradient and the
  pressure jump. */

  foreach_face() {
    if (intface.x[]) {
      if (intside.x[] > 0.)
        uf.x[] -= dt*betac.x[]*(p[] - p[-1] + pjump.x[])/Delta;
      else
        uf.x[] -= dt*betac.x[]*(p[] - p[-1] - pjump.x[])/Delta;
    }
    else
      uf.x[] -= dt*betac.x[]*(p[] - p[-1])/Delta;
  }

  return mgp;
}

/**
In order to avoid code duplication, we overwrite the projection function in the
centered solver, using the new projection defined here which exploits the GFM.
*/

#define project(...) project_ghost(__VA_ARGS__)
#include "navier-stokes/centered.h"
#undef project

/**
The stability condition is copied from [tension.h](/src/tension.h). */

event stability (i++)
{

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  double sigma = 0.;
  //for (scalar c in interfaces)
    //sigma += c.sigma;
  sigma += gfm.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}

/**
At the end of the time step we reconstruct the centered velocity field using
thhe pressure gradient and the interface pressure jump. */

event end_timestep (i++, last) {

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face() {
    if (intface.x[]) {
      if (intside.x[] > 0.)
        gf.x[] = fm.x[]*a.x[] - betac.x[]*(p[] - p[-1] + pjump.x[])/Delta;
      else
        gf.x[] = fm.x[]*a.x[] - betac.x[]*(p[] - p[-1] - pjump.x[])/Delta;
    }
    else
      gf.x[] = fm.x[]*a.x[] - betac.x[]*(p[] - p[-1])/Delta;
  }

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);

  foreach()
    foreach_dimension()
      u.x[] = ubackup.x[] + dt*g.x[];
}

/**
## Notes and Improvements

1. The adaptive grid works only if the interfacial cells used by the GFM are
never unrefined.

2. Reduce code duplication

3. Increasing the density ratio, a static circular droplet shows unphysical
displacements after some simulation time. This problem disappears if the
curvature is constant, and it seems to happen also with
[tension.h](/src/tension.h) or with [integral.h](/src/integral.h), therefore it
is probably a known issue.


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


