/**
This code demonstrates viscous Stokes flow in a microchannel with an advective-diffusive tracer. The flow field in numerically computed for using the navier-stokes centered solver, rather than imposing an analytical expression.

![This plot shows the x-direction velocity in the microchannel.](microchannel_viscous_flow_with_solute/ufx.png)

![This movie shows the solute advection-diffusion over time in the microchannel.](microchannel_viscous_flow_with_solute/c.mp4)
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "diffusion.h" 
#include "run.h"
/**
You must be careful in your choice of time-step which requires some trial and error. 
*/
#define tfinal 60
#define time_step 0.1
/**
For a simple microchannel geometry, define the length and aspect ratio. The higher the aspect ratio, the thinner the channel. 
*/
#define length 10
#define aspect_ratio 30

scalar c[], * tracers = {c};
double width; double yoffset;

int main() {
/**
Set stokes = true to turn off the convective term in the Navier-Stokes solver. Forcing a minimum number of iterations and a specific time-step can be necessary to prevent the solver from being too coarse.
*/
  stokes = true;
  NITERMIN = 2;
  DT = time_step;

  L0 = length; origin (0, 0);
  width = L0/aspect_ratio; yoffset = length*0.5;

  run();
}

event init (t = 0) {
  /**
  Define the microchannel by placing solid boundaries. It can be arbitrarily shaped, and here is a sine wave.
  */
  solid (cs, fs, -(width/2 - y + yoffset + 0.2*width*sin(x/width))*(-width/2 - y + yoffset + 0.2*width*sin(x/width)));
  refine (cs[] && level < 10);
  fractions_cleanup(cs, fs);

  /** 
  Set dirichlet boundary conditions for the solute at the edges of the domain and no-flux conditions at the embedded boundary.
  */
  c[left] = dirichlet(0); c[right] = dirichlet(0); c[embed] = neumann(0);
  p[left] = dirichlet(0)*cs[]; p[right] = dirichlet(0);
  u.n[embed] = dirichlet(0); u.t[embed] = dirichlet(0);
  /**
  Drive the flow with an acceleration term.
  */
  const face vector g[] = {10,0.};
  a = g;
  mu = fm;
  /**
  Intialize the solute as a thin puck to observe dispersion.
  */
  double startpoint = length/4;
  foreach() 
    c[] = (startpoint < x)*(x < startpoint*1.1)*cs[];
}
/**
Assign the advection and diffusion of the solute here. Note that the advection could be automatically done if the tracer module were included at the top, although you would still need to override the diffusion which is zero by default, and furthermore you may wish to impose additional advective effects on the solute.
*/
event advection_and_diffusion (i++) {
  face vector kappa[];
  foreach_face()
    kappa.x[] = fs.x[]*0.01;
  diffusion (c, dt, kappa);
  advection({c}, uf, DT);
}
/**
Create the outputs here.
*/
event solute_movie (i += 2) {
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (c, file = "c.mp4", n=800, mask = m,
    box = {{0, yoffset - width},{length, yoffset + width}});
}

event stop (t = tfinal) {
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (uf.x, file = "ufx.png", n=800, mask = m,
    box = {{0, yoffset - width},{length, yoffset + width}});
}
