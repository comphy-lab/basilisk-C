double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2) {
    mu = new face vector;
    muAx = new face vector;
    muAy = new face vector;
  }

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

#include "fractions.h"


/** 
We need to define the new time-evolving operators 
counting for the modified pressure and viscosity terms:
$$
\textcolor{purple}{ \overline{\textbf{\textsf{A}}} }
\boldsymbol{\nabla}
=
\begin{pmatrix}
\mathrm{e}^{7 \tau /5} \, \partial_x \\
\mathrm{e}^{11 \tau /5} \, \partial_y \\
\end{pmatrix}
$$
for the pressure, and:

$$ 
\boxed{ 
  \textcolor{red}{\mathcal{D}_{\overline{u}_\xi}}
  =
  A_\xi \partial_\xi (2 \mu \partial_\xi \overline{u}_\xi) 
  + A_\eta \partial_\eta (\mu \partial_\eta \overline{u}_\xi) 
  + A_\eta \partial_\eta (\mu \partial_\xi \overline{u}_\eta) 
} 
$$

$$ 
\boxed{ 
  \textcolor{red}{\mathcal{D}_{\overline{u}_\eta}}
  = 
  A_\eta \partial_\xi(\mu \partial_\xi \overline{u}_\eta) 
  + A_\eta \partial_\xi (\mu \partial_\eta \overline{u}_\xi) 
  + A_\eta \partial_\eta(2 \mu \partial_\eta \overline{u}_\eta) 
} 
$$

with:
$$ 
  \begin{aligned} 
    A_\xi(\tau)
    &= 
    \mathrm{e}^{3\tau/5},
    \\ 
    A_\eta(\tau)
    &=
    \mathrm{e}^{7\tau/5}. 
  \end{aligned} 
$$
for the viscosity.

To do so, we take profit of the existing *face vectors* 
`alphav` and `muv`, that we redefine following the two 
directions $\xi$ and $\eta$:

 Required coefficient fields (all face-centred):

  | On x-faces                         | On y-faces                      | 
  | ---------------------------------- | ------------------------------- |
  | `alphav.x = exp(1.4*t)*alphav.x`   | `alphav.y = exp(2.2*t)*alphav.y`|
  | `muv.x`   = `Ax*muv.x`             | `muv.y`   = `Ay*muv.y`          |
  | `muAxv.x` = `Ax*muv.x`             | `muAxv.y` = `Ax*muv.y`          |
  | `muAyv.x` = `Ay*muv.x`             | `muAyv.y` = `Ay*muv.y`          |
*/


event properties (i++)
{
  foreach_face(x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = exp(1.4*t)*fm.x[]/rho(ff);
    // alphav.x[] = fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      face vector muAxv = muAx;
      face vector muAyv = muAy;
      muv.x[] = fm.x[]*mu(ff)*exp(.6*t);
      muAxv.x[] = fm.x[]*mu(ff)*exp(.6*t);
      muAyv.x[] = fm.x[]*mu(ff)*exp(1.4*t);
      // muv.x[] = fm.x[]*mu(ff);
    }
  }
  foreach_face(y) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.y[] = exp(2.2*t)*fm.y[]/rho(ff);
    // alphav.y[] = fm.y[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      face vector muAxv = muAx;
      face vector muAyv = muAy;
      muv.y[] = fm.y[]*mu(ff)*exp(1.4*t);
      muAxv.y[] = fm.y[]*mu(ff)*exp(.6*t);
      muAyv.y[] = fm.y[]*mu(ff)*exp(1.4*t);
      // muv.y[] = fm.y[]*mu(ff);
    }
  }  
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
