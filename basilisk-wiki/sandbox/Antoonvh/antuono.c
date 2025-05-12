/**
# A Steady fully three-dimensional solution to the NS eq.

Matteo Antuono (2020) found a steady 3D solution to the Navier-Stokes
equations for incompressible flow. Here we visualize it.

![Tracers and low-pressure isosurfaces](antuono/mov.mp4)

![Curve stretching](antuono/movc.mp4)

~~~gnuplot
set xlabel 'time'
set ylabel 'Curve length'
set grid
plot 'out' w l lw 2
~~~
 */
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h" 

double k = 1;
Particles parts, curve;

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 2*pi/k;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.05;
  N = 32;
  run();
}

event init (t = 0) {
  TOLERANCE = 1e-4;
  double pref = 4*sqrt(2)/(3*sqrt(3));
  foreach() {
    coord cc = {x, y ,z};
    foreach_dimension() {
      u.x[] = (sin(k*cc.x - 5.*pi/6)*cos(k*cc.y - pi/6)*sin(k*cc.z) -
	       cos(k*cc.z - 5.*pi/6)*sin(k*cc.x - pi/6)*sin(k*cc.y));
      u.x[] *= pref;
    }
  }
  parts = init_tp_square (20);
  int np = 1e5;
  curve = new_tracer_particles (np);
  double th = 0;
  foreach_particle_in(curve) {
    p().x = sin(th);
    p().z = cos(th);
    p().y = 0;
    th += 2*pi/np;
  }
}

event curve_length (t < 5; i++) {
  printf ("%g %g\n", t, plength(curve));
}
event mov (i++) {
  stats ps = statsf (p);
  scalar p2[];
  foreach()
    p2[] = p[] - ps.sum/ps.volume;
  boundary ({p2});
  view (theta = sin(t/10), phi = sin(t/10) + 0.5);
  scatter (parts);
  isosurface ("p2", -0.35);
  box();
  save ("mov.mp4");
  if (t < 5) {
    clear();
    draw_curve (curve);
    box();
    save ("movc.mp4");
  }
}

event stop (t = 10);

/**
## Reference

Antuono, Matteo. "Tri-periodic fully three-dimensional analytic
solutions for the Navier–Stokes equations." *Journal of Fluid
Mechanics* 890
(2020). [PDF](https://www.researchgate.net/publication/339989054_Tri-periodic_fully_three-dimensional_analytic_solutions_for_the_Navier-Stokes_equations)
*/
