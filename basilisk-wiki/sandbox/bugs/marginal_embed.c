/**
# No flux across the embedded boundary for a marginal case.

Following [this discussion](https://groups.google.com/forum/#!topic/basilisk-fr/SLnKb8U0BXs), the embud_flux() function returns zero for marignal cases. 

## Howto Reproduce: 

Modify the [karman.c](/src/examples/karman.c) to use no slip everywhere; 

The channel's walls remain free slip: 

![Animation of the vorticity field.](marginal_embed/vort.mp4)(loop)

![Animation of the tracer field.](marginal_embed/f.mp4)(loop)
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

scalar f[];
scalar * tracers = {f};
face vector muv[];

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  run(); 
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
}

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (0.5 - y, 0.5 + y);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

event movies (i += 4; t <= 15.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, 9, 4);
}