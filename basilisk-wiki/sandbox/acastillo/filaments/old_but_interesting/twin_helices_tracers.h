#include "tracer-particles.h"
#include "scatter2.h"

int n_particles=256;
Particles flow_tracers;

event init (t = 0){
  new_tracer_particles (0);
}

event add_particles (i = 5) {
  if (pid() == 0) {
    flow_tracers  = new_tracer_particles (n_particles);
    int j = 0;
    while (j < n_particles) {
      pl[flow_tracers][j].x  = L0/4.0*(noise()-0.5);
	    pl[flow_tracers][j].y  = L0/4.0*(noise()-0.5);
	    pl[flow_tracers][j].z  = L0/4.0*(noise()-0.5);
	    j++;
    }
  } else {
    flow_tracers  = new_tracer_particles (0);
  }
  particle_boundary (flow_tracers);
  set_particle_attributes (flow_tracers);

}

event movie (t += tsample1) {
  view(camera="iso", fov=4*L0);
  squares ("u.z", linear = false, alpha=-L0/2);
  scatter (flow_tracers,  s = 20);
  box();
  save ("lambda2_particles.mp4");
}
