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

![Vorticity and meterial line](mode3f4/mov.mp4)

~~~gnuplot Length of material line
set xr [0:25]
set yr [-0.1: 15.1]
set size square
set xlabel 'Time [-]'
set ylabel 'L(t) - L(0) [-]'
set grid
plot 'log' u 2:($9- 2.45*3.1415) w l lw 2 t 'Length'
~~~

The `TOLERANCE` controls the divergence, which can be diagnosed
without any discrete approximation.

~~~gnuplot Maximum divergence
set xr [0:25]
set yr [1e-9: 1e-5]
set logscale y
set size square
set xlabel 'Time [T]'
set ylabel 'Divergence [T^{-1}]'
set grid
plot 'log' u 2:10 w l lw 4 t 'Max Divergence', '' u 2:11 w l lw 2 t 'mgp2.resa'
~~~

![Last snapshot](mode3f4/img.png)
 */
#define RKORDER 3
#include "nsf4t.h"
#include "profile6.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

scalar * tracers = NULL;
Particles parts;

long unsigned int np = 1e5;
double P = 1e-5, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameter

int maxlevel = 8;

double length (Particles p) {
  double lenp = 0;
  coord p1 = {0,0}, pp = {0};
  foreach_particle_in(p) {
    if (p1.x == 0 && p1.y == 0)
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
  TOLERANCE = 1e-4;
  parts = new_tracer_particles (np);
  double Theta = 0;
  foreach_particle_in(parts) {
    p().x = (b + 1.)/2.*cos(Theta);
    p().y = (b + 1.)/2.*sin(Theta);
    Theta += 2*pi/(double)(np + 1);
  }
  refine (sq(x) + sq(y) < sq(b*1.2) && level < maxlevel + 2);
  printf ("# Initializing using %ld cells...\n", grid->tn);
  double betav = 1./(1 - sq(b)), alphav = -betav*sq(b);
  vector uc[], ucf[];
  foreach() {
    double rp = RAD*RADP(P,m), vr = 0;
    if (rp <= 1.)
      vr = rp;
    else if (rp > 1 && rp <= b) 
      vr = alphav/rp + betav*rp;
    uc.x[] = -y/rp*vr;
    uc.y[] =  x/rp*vr;
    ucf.x[] = ucf.y[] = 0;
  }
  boundary ((scalar*){uc, ucf});
  /**
     The Helmholtzfilter is applied to smoothen the piecewise flow
     profile.
  */
  const face vector alphaf[] = {-sq(0.4/(2.*pi)), -sq(0.4/(2.*pi))};
  foreach_dimension()
    poisson (ucf.x, uc.x, alphaf, unity);  
  foreach_face()
    u.x[] = ((-ucf.x[-2] + 7*(ucf.x[-1] + ucf.x[]) - ucf.x[1])/12.);
  scalar vtb[], vta[], * ab = {vtb, vta};
  foreach() {
    vtb[] = x/RAD*uc.y[] - y/RAD*uc.x[];
    vta[] = x/RAD*ucf.y[] - y/RAD*ucf.x[];
  }
  profile (ab, sqrt(sq(x) + sq(y)), "prof");
  boundary ((scalar*){u});
  project (u, p2);
  unrefine (level >= maxlevel);
}

event logger (i++) {
  boundary_flux ({u});
  scalar div[];
  foreach() {
    div[] = 0;
    foreach_dimension()
      div[] += (u.x[1] - u.x[0]);
    div[] = fabs(div[]/Delta);
  }
  stats ds = statsf (div);
  fprintf (stderr, "%d %g %d %d %d %d %ld %d %g %g %g\n", i, t, mgp.i, 
           mgp.nrelax, mgp2.i, mgp2.nrelax, grid->tn,
	   grid->maxdepth, length (parts), ds.max, mgp2.resa);
}

event mov (t += 0.5) {
  scalar omg[];
  vorticityf (u, omg);
  squares ("omg", linear = true, map = cool_warm);
  scatter (parts);
  save ("mov.mp4");
}

event stop (t = 25) {
  clear();
  scalar omg[];
  vorticityf (u, omg);
  view (fov = 9.39975, samples = 4);
  squares ("omg", map = cool_warm,
           min = -2, max = 2, linear = true);
  save ("img.png");
}

/**
# Reference

G. R. Flierl, "On the instability of geostrophic vortices", *J. Fluid
Mech. 197, 349*
*/
