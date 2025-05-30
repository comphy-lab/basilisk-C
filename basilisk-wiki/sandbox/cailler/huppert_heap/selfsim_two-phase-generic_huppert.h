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
  
  if (mu1 || mu2)
    mu = new face vector;

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
\left[
\textcolor{red}{ \overline{\textbf{\textsf{M}}} }
.\boldsymbol{\nabla}
\right] \cdot
=
{}^{\mathrm{t}}\left(
  \begin{pmatrix}
    \mathrm{e}^{3 \tau /5} \, \partial_x &
    \mathrm{e}^{7 \tau /5} \, \partial_y 
  \end{pmatrix}
  {}^{\mathrm{t}} ( \,\, )
\right)
$$
for the viscosity. 
In practice, for the viscosity term, we have the following 
equivalence:
$$
\left[
\textcolor{red}{ \overline{\textbf{\textsf{M}}} }
.\boldsymbol{\nabla}
\right] \cdot \left(
  2 \, \overline{\textbf{\textsf{S}}}
\right) 
= 
\left(
  \textcolor{red}{ \overline{\mathbf{M}}' }
  \cdot 
  \boldsymbol{\nabla}^2
\right) 
\overline{\mathbf{u}}, 
\quad \text{where } \,
\textcolor{red}{ \overline{\mathbf{M}}' 
=
\begin{pmatrix}
  \mathrm{e}^{3 \tau /5} \\
  \mathrm{e}^{7 \tau /5}
\end{pmatrix}
}
$$

To do so, we take profit of the existing *face vectors* 
`alphav` and `muv`, that we redefine following the two 
directions $\xi$ and $\eta$:
*/


event properties (i++)
{
  foreach_face(x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = exp(1.4*t)*fm.x[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff)*exp(.6*t);
    }
  }
  foreach_face(y) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.y[] = exp(2.2*t)*fm.y[]/rho(ff);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.y[] = fm.y[]*mu(ff)*exp(1.4*t);
    }
  }  
  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
