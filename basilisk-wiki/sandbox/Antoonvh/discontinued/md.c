/**
# A 2D Leonnard-Jones solid

![Closed packed atoms](md/mov.mp4)
 */
#include "grid/multigrid.h"
#include "../inertial-particles.h"
#include "../particle_reference.h"
#include "view.h"
#include "../scatter2.h"

scalar parts[];
Particles part;

coord G = {0, 0};
double sqdistance (particle a1, particle a2) {
  double d = 0;
  foreach_dimension()
    d += sq(a1.x - a2.x);
  return d > 0 ? d : 0;
}

double force(double qd) {
  return (48*(1./pow(qd, 6.5) - 0.5/pow(qd, 3.5)));
}

coord p_acc (struct PA_inp inp) {
  particle p1 = inp.p;
  coord a = G;
  Point point = locate (p1.x, p1.y, p1.z);
  foreach_neighbor(1) {
    foreach_particle_point(parts) {
      double d = sqdistance (p1, p());
      if (d) {
	coord dir;
	foreach_dimension() {
	  dir.x = d < sq(L0)/4 ? (p1.x - p().x) : 1.;
	}
	normalize (&dir);
	foreach_dimension()
	  a.x += dir.x*force(sqdistance (p1, p()));
      }
    }
  }
  return a;
}

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 70;
  X0 = Y0 = -L0/2.;
  DT = 0.01;
  N = 4;
  run();
}

event init (t = 0) {
  part = init_ip_square (n = 15, l = 16.5);
  foreach_particle() {
    foreach_dimension()
      p().x += 0.01*noise();
  }
}

event timestep(i++) {
  particle_boundary (part);
  assign_particles (part, parts);
  dt = dtnext(DT);
}

event mov (i+= 20) {
  scatter(part);
  box();
  //cells();
  save ("mov.mp4");
}

event stop (t = 100) {
  ;
}
