/**
This code demonstrates a bug where the no-slip condition is violated at the embedded boundary in pressure-driven viscous microchannel flow. It is similar to [Antoon van Hooft's code for scalar mixing in a tube flow](http://basilisk.fr/sandbox/Antoonvh/pois-front.c) where the key difference is here I numerically compute the flow field instead of prescribing the exact solution.

![This plot shows the spatial distribution of the solute at the final simulated time step.](marginal_embed_2/c.png)

![This plot shows the spatial distribution of the x-component of the velocity at the final simulated time step. I expect a fully developed parabolic flow with zero slip at the embedded boundaries.](marginal_embed_2/ux.png)

![This plot shows the spatial distribution of the pressure field at the final simulated time step.](marginal_embed_2/p.png)

![This plot shows the spatial distribution of the y-component of the velocity at the final simulated time step.](marginal_embed_2/uy.png)
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "diffusion.h" 
#include "run.h"

#define tfinal 1
#define time_step 0.1
#define length 5
#define aspect_ratio 10

scalar c[], * tracers = {c};

int main() {
  /**
  The goal is to simulate viscous Stokes flow, although the convective term of the Navier-Stokes equations can be turned on (leaving the default stokes = false) for debugging purposes.
  */
  # if 1
  stokes = true;
  //TOLERANCE = HUGE;
  NITERMIN = 2;
  DT = time_step;
  # endif
  
  /**
  The suggested workaround for this [already-documented bug](https://groups.google.com/g/basilisk-fr/c/BTIpES3crlY/m/d5IMVjN7BQAJ) to adjust the domain size or origin did not do the trick.
  */

  L0 = length; origin (0, 0);
  run();
}

event init (t = 0) {
  /** You can define microchannel geometry here. As suggested in a [previous post](https://groups.google.com/g/basilisk-fr/c/BTIpES3crlY/m/d5IMVjN7BQAJ), moving the boundary (try adjusting the multiplying factor in yoffset) did not help.
  */
  double width = L0/aspect_ratio; double yoffset = length*0.5;
  solid (cs, fs, -(width/2 - y + yoffset)*(-width/2 - y + yoffset));
  refine (cs[] && level < 10);
  fractions_cleanup(cs, fs);

  // boundary conditions
  c[left] = dirichlet(0); c[right] = dirichlet(0); c[embed] = neumann(0);
  
  /**
  Both BC's #1 and #2 led to some instabilities near the channel inlet on the left, and have very different flow fields despite similar pressure fields.
  */
  // BC #1: assign a pressure gradient
  //p[left] = dirichlet(1)*cs[]; p[right] = dirichlet(0);
  
  // BC #2: uniform inflow and outflow
  //u.n[left] = dirichlet(1)*cs[]; u.n[right] = dirichlet(1)*cs[];
  
  /** The no-slip condition is implemented here but seems to be ignored by the solver. */
  // no-slip at embedded boundaries
  u.n[embed] = dirichlet(0); u.t[embed] = dirichlet(0);

  /** Alternatively to either BC #1 or #2, you can drive flow through this acceleration term (making sure to select BC #1). In this case, the no-slip condition is observed, however if the advective term of the solute is turned on, the concentration c goes to zero everywhere. */
  # if 1
  // using acceleration term to assign dpdx
  const face vector g[] = {-1.,0.};
  a = g;
  mu = fm;
  # endif

  /** Initialize the solute as a thin plug. */
  // intialize solute cloud
  double startpoint = length/4;
  foreach() 
    c[] = (startpoint < x)*(x < startpoint*1.1)*cs[];
}

/** The advection and diffusion of the scalar are enforced here, although while debugging the flow field I comment out the advection call. */
// simple advection diffusion of passive scalar c
event advection_and_diffusion (i++) {
  face vector kappa[];
  foreach_face()
    kappa.x[] = fs.x[]*0.01;
  diffusion (c, dt, kappa);
  //advection({c}, u, DT);
}

// Generate images at the final timestep.
event stop (t = tfinal) {
  output_ppm (c, file = "c.png", n=400);
  output_ppm (u.x, file = "ux.png", n=400);
  output_ppm (u.y, file = "uy.png", n=400);
  output_ppm (p, file = "p.png", n=400);
}
