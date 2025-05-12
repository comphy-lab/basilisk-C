/**
![An all-liquid planetary core cools down at the edge. Image via
[Express
news](www.express.co.uk/).](https://cdn.images.express.co.uk/img/dynamic/151/590x/Earth-core-leaking-geology-study-Earth-molten-core-mantle-exchange-tungsten-1160314.jpg)

## Liquid core convection in *a* planet.

A sphere with radius $R$ contains a boussinesq fluid with viscosity
$\nu$ and reference density $\rho_0$. The fluid is subject to its own
gravity potential $\Phi$ and therefore has a quadratic variation with
the radial coordinate ($r$), $\Phi = A r^2$, with $A$ a constant
($2\pi G\rho_0/$`dimension`). For generaltiy we use a $b$ ($\langle b\rangle =
t^{-2}$) field to describe the effect of gravity:

$$\mathbf{F_g} = b \mathbf{r},$$

which eliminates $\rho_0$ as a system parameter.

The fluid is cooled at the sphere's edge and thereby becomes slighly
more dense. We use the (nagative) $b$-flux ($B$) to
describe the cooling. The medium has a $b$-diffussivity $\kappa$.

With system parameters $R, \nu, B,\kappa$ and the
dimensionality we can construct the following dimensionless groups:

$$\Pr = \frac{\nu}{\kappa},$$  

$$Re = \frac{\|B\|^{1/3}R^{8/3}}{\nu},$$  

$$\Pi = \texttt{dimension}.$$    
  
We choose $Pr = 1$ and $Re = 30$. 

![The 2D solution for the $b$ field and some particle flow
 tracers](core/mov2D.mp4)

![The 3D solution for the $b$ field (on a hemisphere) and some
 particle flow tracers](core.mov3D.mp4)

Note that the earth has a solid core in a liquid mantle.
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
//#include "view2.h"     //does not draw_vof facets for z > 0 in 3D
#include "view.h"
#define BVIEW 1
#include "particles.h"

#define RADIUS (sqrt(sq(x) + sq(y) + sq(z)))

scalar b[], * tracers = {b};
face vector muv[], av[];
double muf = 1./30., R = 1., B = -1.;

b[embed]   = neumann (B/muf);
u.t[embed] = dirichlet (0.);
u.n[embed] = dirichlet (0.);
#if dimension == 3
u.r[embed] = dirichlet (0.);
#endif

// Relaxation time (non physical)
double tau = 10.;

int main() {
  L0 = R*2.1;
  X0 = Y0 = Z0 = -L0/2.;
  mu = muv;
  a = av;
  run();
}

event init (t = 0) {
  scalar phi[];
  foreach_vertex()
    phi[] = R - RADIUS;
  fractions (phi, cs, fs);
  init_particles_2D_square_grid (6, 0, 0, R);
}

event properties (i++) {
  foreach_face()
    muv.x[] = fs.x[]*muf;
  boundary ((scalar*){muv});  
}

event acceleration (i++) {
  foreach_face() {
    coord f = {x, y, z};
    av.x[] = face_value(b,0)*f.x;
  }
}

event tracer_diffusion (i++)
  diffusion (b, dt, muv);
/**
We modify the buoyancy to omit large negative values. This "heating"
does not influence the dynamics and should balance the cooling in a
quasi steady state.
 */
event no_run_away (i++) {
  double tendency = dimension*B/R;
  double mean = statsf(b).sum/statsf(b).volume;
  foreach() 
    b[] -= (cs[] > 0)*dt*mean/tau;
  boundary ({b});  
  printf ("%g %g cf.: %g\n",t, mean, tendency*tau);
}

#if TREE
event adapt (i++) 
  adapt_wavelet ({b, u}, (double[]){0.4, 0.2, 0.2, 0.2}, 6, 4);
#endif

event movie (t += .1) {
  view (fov = 19);
#if dimension == 2
  scatter (loc);
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("b", linear = true);
#else
  view (theta = sin(t/10.), phi = cos(t/15));
  scatter (loc, s = 30);
  scalar phi[];
  foreach_vertex()
    phi[] = R - RADIUS;
  boundary ({phi});
  scalar c[];
  face vector f[];
  fractions (phi, c, f, 0.05);
  draw_vof ("c", "f", color = "b", linear = true);
#endif
  char mnane[99];
  sprintf (mnane, "mov%dD.mp4", dimension);
  save (mnane);
}

event stop (t = 100);
