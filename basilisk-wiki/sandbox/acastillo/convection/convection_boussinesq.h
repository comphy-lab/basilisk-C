/**
\renewcommand*{\vec}[1]{\boldsymbol{#1}}

# Boussinesq equations

We use the dimensionless Boussinesq equations written in the following form:

$$
\begin{aligned}
\nabla \cdot \vec{u} =& 0
\\
\frac{\partial \vec{u}}{\partial t} + \nabla \cdot \left(\vec{u} \otimes \vec{u}
\right) =& -\nabla p + \text{Pr}\text{Ra}^{0.5} \nabla^2 \vec{u} + \text{Pr}\theta\vec{e}_y
\\
\frac{\partial\theta}{\partial t} + \nabla \cdot \left( \vec{u} \theta \right)
=& \text{Ra}^{-0.5} \nabla^2 \theta \end{aligned}
$$

Under the Boussinesq approximation our system is fully defined in terms of two
dimensionless parameters: the *Rayleigh* and *Prandtl* numbers, in addition to
the geometry and boundary conditions.

Instead of writing new code, we combine the (incompressible) Navier-Stokes
equations with an advection-diffusion problem for the temperature field.
*/

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar T[];
scalar * tracers = {T};
mgstats mgT;

face vector D[];

double Ra, Pr;


/** 
  We use the *Rayleigh* and *Prandtl* numbers to initalize the viscosity field
*/

event init (i=0) {

  a = new face vector;
  mu = new face vector;

  face vector muc = mu;
  foreach_face(){
    muc.x[] = fm.x[]*Pr/sqrt(Ra);
    D.x[] = fm.x[]*1./sqrt(Ra);
  }
}

/**
  and include the diffusion of the temperature field
*/

event tracer_diffusion (i++) {
  mgT = diffusion (T, dt, D);
}

/**
Finally, we add the bouyancy to the acceleration term
*/
event acceleration (i++) {
  face vector av = a;
  #if dimension == 2
  foreach_face(y)
    av.y[] += Pr*(T[] + T[0,-1])/2.;
  #else
  foreach_face(z)
    av.z[] += Pr*(T[] + T[0,0,-1])/2.;
  #endif
}

/** and clean-up the variables we allocated manually */
event cleanup (i = end) {
  delete ((scalar*){a,mu});
}