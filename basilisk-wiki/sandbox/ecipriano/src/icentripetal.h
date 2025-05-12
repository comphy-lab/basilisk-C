/**
# Interfacial Suspending Force

We re-express the suspending force in [centripetal.h](centripetal.h)
as an [interfacial force](/src/iforce.h). The suspending
force is expressed as a function of the gradient of a potential:
$$
\mathbf{f}_m = \rho \alpha \nabla \xi
$$
As explained by [Popinet 2018](#popinet2018), a source term expressed
as a gradient of a potential can be discretized using a well-balanced
scheme, similarly to the surface tension force.  Setting $\mathbf{u}=0$
in the Navier-Stokes equations:
$$
-\nabla p + \rho\nabla\xi =
-\nabla p + \nabla(\rho\xi) - \xi\nabla\rho =
-\nabla p' - \xi\nabla\rho
$$
with $\rho$ which jumps from a constant value $\rho_l$
in liquid phase, to the value $\rho_g$ in gas phase
according to:
$$
\rho(H) = (\rho_l - \rho_g)H +\rho_l = [\rho]H + \rho_l
$$
whose gradient reads:
$$
\nabla\rho = [\rho]\nabla H
$$
Therefore, the body force $\mathbf{f}_m$ is replaced by
an interfacial force, using an approach similar to the
one used in [tension.h](/src/tension.h) and [reduced.h](/src/reduced.h):
$$
-\nabla p' - \xi\nabla\rho =
-\nabla p' - \phi\nabla H
$$
where the term $\phi = \xi[\rho]$ is computed here and added
to the surface forces by [iforce.h](iforce.h). The potential
is computed as:
$$
\xi = \dfrac{\xi_0}{|\mathbf{x}_\Gamma - \mathbf{x}_c|}
$$
where $\xi_0$ is an arbitrary constant that regulates the
intensity of the force; $\mathbf{x}_\Gamma$ is the center of
the interface element; $\mathbf{x}_c$ is the point where the
center of the droplet should remain fixed.
*/

#define CENTRIPETAL

/**
We need the interfacial force module as well as some
functions to compute the position of the interface. */

#include "iforce.h"
#include "curvature.h"

/**
## User Data

We need to define the coordinate of the point where the center of
the droplet should remain *p*, and a parameter that controls the
intensity of the suspending force *eps*.
*/

struct SuspendingForceModel {
  coord p;
  double eps;
  double sigma;
};

struct SuspendingForceModel sfm = {
  .p = {0., 0., 0.},
  .eps = 1.25e-4,
  .sigma = 0.,
};

/**
# Position of an interface

This is similar to reduced but instead of computing the position
of the interface we calculate the potential:
$$
\xi = \dfrac{\xi_0}{|\mathbf{x}_\Gamma - \mathbf{x}_c|}
$$

This is defined only in interfacial cells. In all the other cells it
takes the value *nodata*.

We first need a function to compute the magnitute of the distance
between the interface and the suspending center $|\mathbf{x}_\Gamma
- \mathbf{x}_c|$ of an interface. For accuracy, we first try to
use height functions. */

foreach_dimension()
static double pos_centripetal_x (Point point, vector h)
{
  if (fabs(height(h.x[])) > 1.)
    return nodata;
  coord o = {x, y, z};
  o.x += height(h.x[])*Delta;
  double pos = 0.;
  foreach_dimension()
    pos += sq (o.x - sfm.p.x);
  return sqrt (pos);
}

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the position. This is done by the function below which returns
the HF position given a volume fraction field *f* and a height function
field *h*.*/

static double height_position_centripetal (Point point, scalar f, vector h)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *f*) and corresponding HF position function *pos*
  (defined above). */

  typedef struct {
    double n;
    double (* pos) (Point, vector);
  } NormPos;
  struct { NormPos x, y, z; } n;
  foreach_dimension()
    n.x.n = f[1] - f[-1], n.x.pos = pos_centripetal_x;

  /**
  We sort these pairs in decreasing order of $|n|$. */

  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormPos, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormPos, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormPos, n.y, n.z);
#endif

  /**
  We try each position function in turn. */

  double pos = nodata;
  foreach_dimension()
    if (pos == nodata)
      pos = n.x.pos (point, h);

  return pos;
}

void position_centripetal (scalar f, scalar pos, bool add = false)
{

  /**
  On trees we set the prolongation and restriction functions for
  the position. */
  
#if TREE
  pos.refine = pos.prolongation = curvature_prolongation;
  pos.restriction = curvature_restriction;
#endif

  vector fh = f.height, h = automatic (fh);
  if (!fh.x.i)
    heights (f, h);
  foreach() {
    if (interfacial (point, f)) {
      double hp = height_position_centripetal (point, f, h);
      if (hp == nodata) {
        coord n = mycs (point, f), o = {x,y,z}, c;
        double alpha = plane_alpha (f[], n);
        plane_area_center (n, alpha, &c);

        hp = 0.;
        foreach_dimension()
          hp += sq (o.x + Delta*c.x - sfm.p.x);
        hp = sqrt (hp);
      }

      if (add)
        pos[] += 1./hp*sfm.eps*(rho2v[] - rho1v[]);
      else
        pos[] = 1./hp*sfm.eps*(rho2v[] - rho1v[]);
    }
    else
      pos[] = nodata;
  }
}

/**
## Stability condition

We add the possibility to compute the stability conditions based on
the time step required for the time-explicit discretization of the
surface tension force. This allows to have a reliable time-step when
using this module without surface tension. The coefficient $\sigma$
is set as *sfm.sigma*. */

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

  double sigma = sfm.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}

/**
We overload the acceleration() event to add the contribution of
gravity to the interfacial potential $\phi$.

If $\phi$ is already allocated, we add the suspending force potential,
otherwise we allocate a new field and set it to the contribution of
the suspending force potential. */

event acceleration (i++)
{
  scalar phi = f.phi;
  if (phi.i)
    position_centripetal (f, phi, add = true);
  else {
    phi = new scalar;
    position_centripetal (f, phi, add = false);
    f.phi = phi;
  }
}

/**
## Notes

I noticed that this model does not work properly with AMR, except if
the adapt_wavelet_leave_interface function is used in order to avoid
coarsening of interface cells. Maybe different refinement and coarsening
functions must be used for pos?

## References

~~~bib
@hal{popinet2018, hal-01528255}
@article{saufi2019dropletsmoke++,
  title={DropletSMOKE++: a comprehensive multiphase CFD framework for the evaporation of multidimensional fuel droplets},
  author={Saufi, Abd Essamade and Frassoldati, Alessio and Faravelli, T and Cuoci, A},
  journal={International Journal of Heat and Mass Transfer},
  volume={131},
  pages={836--853},
  year={2019},
  publisher={Elsevier}
}
~~~
*/
