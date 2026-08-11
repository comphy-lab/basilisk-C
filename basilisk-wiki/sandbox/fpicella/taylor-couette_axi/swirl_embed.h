/**
# Azimuthal velocity for axisymmetric flows -- EMBED-safe version

This is an embed-safe version of Basilisk's standard
navier-stokes/swirl.h.

Coordinates:
    x = axial coordinate z
    y = radial coordinate r
    w = azimuthal velocity u_theta

The modification with respect to the standard implementation is
the evaluation of d(mu)/dy in tracer_diffusion(): divisions by
fm.y are performed only when the corresponding face metrics are
strictly positive.
*/

#include "tracer.h"
#include "diffusion.h"

scalar w[], * tracers = {w};

/*
Azimuthal velocity is zero on the symmetry axis.
*/
w[bottom] = dirichlet (0.);


/**
Allocate the acceleration field if necessary.
*/

event defaults (i = 0)
{
  if (is_constant (a.x)) {
    a = new face vector;

    foreach_face()
      a.x[] = 0.;
  }
}


/**
Centrifugal acceleration:

        a_r = w^2/r

The velocity is averaged onto radial faces.
*/

event acceleration (i++)
{
  face vector av = a;

  foreach_face (y)
    av.y[] +=
      y > 0. ?
      sq(w[] + w[0,-1])/(4.*y) :
      0.;
}


/**
Diffusion of the azimuthal velocity.

For axisymmetric swirl,

  theta = rho
  D     = mu

with rho and mu already containing the axisymmetric metric.

The only difference from standard swirl.h is the protected
evaluation of d(mu)/dy.
*/

event tracer_diffusion (i++)
{
  scalar beta[], theta[];

  foreach() {

    theta[] = rho[];

    double muc =
      (mu.x[] +
       mu.x[1] +
       mu.y[] +
       mu.y[0,1])/4.;

    /*
    ------------------------------------------------------
    EMBED-safe evaluation of d(mu)/dy
    ------------------------------------------------------

    Standard swirl.h effectively evaluates

        mu.y/fm.y

    on both neighboring radial faces.

    With embedded boundaries fm.y may vanish.

    If either face is blocked, do not evaluate the
    derivative across that solid boundary.
    */

    double dymu = 0.;

    if (fm.y[] > 0. &&
        fm.y[0,1] > 0.) {

      dymu =
        (mu.y[0,1]/fm.y[0,1] -
         mu.y[]    /fm.y[]) /
        Delta;
    }

    /*
    Axisymmetric source term.

    Protect the symmetry axis explicitly.
    */

    if (y > 0.)
      beta[] =
        -(rho[]*u.y[] + muc/y)/y
        - dymu;
    else
      beta[] = 0.;
  }

  diffusion (w, dt, mu,
             theta = theta,
             beta = beta);
}