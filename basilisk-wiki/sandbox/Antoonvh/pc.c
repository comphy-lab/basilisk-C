/**
# Flow tracers and inertial particles

![Heavy magenta particles are spun out of vortices](pc/parts.mp4)

*/
#include "navier-stokes/stream.h"
#include "tracer-particles.h"
#include "stokes-particles.h"
#include "view.h"
#include "scatter2.h"

Particles flow, heavy;

// We do not use the values of these fields:
face vector mu[];
scalar rho[];

int main() {
  X0 = Y0 = -L0/2;
  N = 64;
  DT = 0.25;
  run();
}

event init (t = 0) {
  new_tracer_particles (0);
  new_inertial_particles (0);
  foreach() 
    omega[] = noise();
  refine (level < 7);
}

event add_p (t = 1) {
  flow  = init_tp_cells();
  heavy = init_ip_cells();
  foreach_particle_in(heavy) 
    p().u2.z = 1.;  //Set relaxation timescale 
}

event bviewer (i += 4) {
  if (flow) {
    squares ("omega");
    scatter (heavy, s = 5, pc = {1, 0, 1});
    translate (z = 0.02) {
      scatter (flow , s = 2);
    }
    save ("parts.mp4");
  }
}

event stop (t = 250);

