/**
# A mode 3 instability 

Unstable vortex following Flierl 

~~~gnuplot Inital profile
set xr [0:2]
set yr [-0.1: 1.1]
set size square
set ylabel 'v_{theta} [-]'
set xlabel 'r [-]'
set grid
plot 'prof' u 1:2 w l lw 2 t 'Piecewise',\
         '' u 1:3 w l lw 3 t 'Initialized'
~~~

![Vorticity and material line](mode3/mov.mp4)

~~~gnuplot Length of material line
set xr [0:25]
set yr [-0.1: 15.1]
set size square
set ylabel 'Time [-]'
set xlabel 'L(t) - L(0) [-]'
set grid
plot 'out' u 1:($2- 2.45*3.1415) w l lw 2 t 'Length'
~~~

We can observe an errouneous mode 4 instability
 */
#include "navier-stokes/centered.h"
#include "diffusion.h" // it is a Helmholtz filter
#include "profile6.h"
#include "higher-order.h"
#define interpolate_linear interpolate_quartic
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"
Particles parts;
long unsigned int np = 1e5;
double P = 1e-5, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameters

int maxlevel = 8;

double length (Particles p) {
  double lenp = 0;
  coord p1 = {pi, exp(1)}, pp = {0};
  foreach_particle_in(p) {
    if (p1.x == pi && p1.y == exp(1))
      p1 = (coord){x, y};
    else {
      double d = 0;
      foreach_dimension()
	d += sq(p().x - pp.x);
      lenp += sqrt(d);
    }
    pp = (coord){x, y};
  }
  double d = 0;
  foreach_dimension()
    d += sq(p1.x - pp.x);
  lenp += sqrt(d);
  return lenp;
}

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 8.;
  origin (-pi*4./3., -2*exp(1) + 1.2); // Not centered
  N = 1 << maxlevel;
  run();
}

#define RAD (sqrt((sq(x) + sq(y))))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M))))

event init (t = 0) {
  parts = new_tracer_particles (np);
  double Theta = 0;
  foreach_particle_in(parts) {
    p().x = (b + 1.)/2.*cos(Theta);
    p().y = (b + 1.)/2.*sin(Theta);
    Theta += 2*pi/(double)(np + 1);
  }
  TOLERANCE = 1e-4;
  double betav = 1./(1 - sq(b)), alphav = -betav*sq(b);
  refine (RAD < 1.2*b && level < maxlevel + 2);
  foreach() {
    double rp = RAD*RADP(P,m), vr = 0;
    if (rp <= 1.)
      vr = rp;
    else if (rp > 1 && rp <= b) 
      vr = alphav/rp + betav*rp;
    u.x[] = -y/rp*vr;
    u.y[] =  x/rp*vr;
  }
  /**
The Helmholtzfilter is applied to smoothen the discontinuous flow
profile.
  */
  scalar vtb[], vta[], * ab = {vtb, vta};
  foreach() 
    vtb[] = x/RAD*u.y[] - y/RAD*u.x[];
  const face vector muc[] = {sq(0.4/(2*pi)), sq(0.4/(2*pi))};
  foreach_dimension()
    diffusion (u.x, 1., muc);
  foreach() 
    vta[] = x/RAD*u.y[] - y/RAD*u.x[];
  profile (ab, sqrt(sq(x) + sq(y)), "prof");
  unrefine (level >= maxlevel);
}

event line_len (t += 0.1) 
  printf ("%g %g\n", t , length(parts));

event mov (t += .5) {
  scalar omg[];
  vorticity (u, omg);
  output_ppm (omg, file = "o.mp4", n = 300, linear = true,
	      min = -1.5, max = 1.5, map = cool_warm);
  squares ("omg", linear = true, map = cool_warm);
  scatter (parts);
  save ("mov.mp4");
}

event stop (t = 25);
/**
# Reference

G. R. Flierl, "On the instability of geostrophic vortices", *J. Fluid
Mech. 197, 349*
*/
