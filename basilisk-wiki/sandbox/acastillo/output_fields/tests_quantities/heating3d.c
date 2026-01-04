/**
# A rising plume 

Simulates a thermal plume from a Gaussian heat source using buoyancy-driven flow to illustrate the reference height field.

Based on [heating.c](/sandbox/Antoonvh/heating.c) by Antoon van Hooft.

![Buoyancy field](heating/field_b.mp4)
![Reference height](heating/field_yref.mp4)

*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "acastillo/output_fields/available_potential.h"
#include "view.h"

// Fields
scalar b[];               // Buoyancy field
scalar yref[];            // Reference height      
scalar * tracers = {b};
face vector av[];         

// Parameters
double siz = 0.3;         // Heat source size
double muv = 2e-4;        // Diffusivity coefficient
double be = 2e-3, ue = 0.05;  // Adaptation tolerances
int maxlevel = 7, minlevel = 4;

int main() {
  L0 = 6.0;
  X0 = Y0 = Z0 = -L0/2;        // Center domain at origin
  a = av;                      
  DT = 0.1;
  N = 1 << minlevel;
  const face vector muc[] = {muv, muv, muv};
  mu = muc;
  run();
}

#define RAD2(x,y,z) (sq(x) + sq(y) + sq(z))

event init (t = 0) {
  b.gradient = minmod2;
  
  // Nested refinement around heat source
  refine (RAD2(x,y,z) < sq(20*siz) && level < maxlevel - 2);
  refine (RAD2(x,y,z) < sq(10*siz) && level < maxlevel - 1);
  refine (RAD2(x,y,z) < sq(2*siz) && level < maxlevel);

  foreach(){
    b[] = exp(-RAD2(x,y,z)/sq(siz));  // Gaussian buoyancy source
  }
}

event acceleration(i++) {
  double tau = 1.0;
  
  foreach(){
    // Continuous heat source
    b[] += (1.0 - b[]) * dt/tau * exp(-RAD2(x,y,z)/sq(siz));
    b[] = clamp(b[], 0.0, 1.0);
  }
  
  // Buoyancy force
  foreach_face(y) 
    av.y[] = (b[] + b[0,-1])/2.0;
}
// Thermal diffusion
event tracer_diffusion (i++) {
  diffusion (b, dt, mu);
}

// Update reference height for available potential energy
event update_yref (i++) {
  reference_height(yref, b, 0.0, 1.0, store=false, reverse=false);
}

// Adaptive mesh refinement
event adapt (i++) {
  adapt_wavelet({b, u, yref}, (double[]){be, be, ue, ue, ue}, maxlevel, minlevel);
}

// Movie output
event mov (t += 0.02; t <= 5.0) {
  squares("b", linear = false, spread = -1);
  save ("field_b.mp4");

  squares("yref", linear = false, spread = -1);
  save ("field_yref.mp4");
}

#undef RAD2