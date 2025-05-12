/**
# Greenhouse convection

We present a simple flow setup to study the flow topology of
convection in a two-dimensional greenhouse gemoetry.
 
 ![the resulting movie](greenhouse/bb.mp4)
 */

#include "embed.h"
#include "navier-stokes/centered.h"
/**
   Buoyancy is used as the themodynamic variable
 */
#include "tracer.h"
#include "diffusion.h"
#include "view.h"

scalar b[], *tracers = {b};

/**
Values for the viscosity and diffusivity are set and the related
fields are delared together with a field to store the buoyancy
acceleration.
 */

double nu = 1./500., kappa = 1./500.;

face vector nuf[], kappaf[], av[];

/**
   A heat-flux boundary condition for the bottom boundary, and a fixed
   "temperature" at the roof.
*/
double B = 0.01;
b[bottom] = neumann (B/kappa);
b[embed] = dirichlet (0);

/**
Furthermore, a no slip condition is imposed for the flow
 */
u.t[bottom] = dirichlet (0);
u.t[embed] = dirichlet (0); // actually u.x
u.n[embed] = dirichlet (0); // actually u.y

int main() {
  periodic (left);
  L0 = 10;  // domain size
  DT = 0.1; // Initial timestep 
  mu = nuf; // set viscosity
  a = av;   // set Buoyancy acceleration
  run(); 
}

event init (t = 0) {
  solid (cs, fs, 0.5*cos(6*pi/L0*x) + -y + 3); // the roof geometry 
  foreach() 
    b[] = TOLERANCE*noise()*cs[]; // Perturbation
}

/**
The details of the embedded boundary implementation require the face
fields to be zero outside the domain.
 */
event properties (i++) {
  foreach_face() {
    nuf.x[] = nu*fs.x[];
    kappaf.x[] = kappa*fs.x[];
  }
}

/**
The buoyancy field is directly related to the vertical body force
 */
event acceleration (i++) {
  coord g = {0, 1};
  foreach_face()
    av.x[] = g.x*(b[] + b[-1])/2.;
}
/**
Advection of the buoyancy field is automatically traken care of, but
the diffusion step needs to written here
 */
event trace_diffusion (i++) {
  diffusion (b, dt, kappaf);
}

/**
The mesh is adapted following the wavelet-based estimation for the
discretization error.
 */
event adapt (i++) {
  adapt_wavelet ({cs, b, u}, (double[]){0.001, 0.01, 0.01, 0.01, 0.01}, 7);
}

/**
## Output

The output consists of two movies:

* A simple one
 */

event movie (t += 0.5; t <= 400) {
  output_ppm (b, file = "b.mp4", n = 400,
	      linear = true, min = 0, max = 1);
}

/**
 * and a fancy one (using `bview`)
 */

event bmovie (t += 0.5) {
  view (tx = -0.5, ty = -0.2, width = 500, height = 300, fov = 15);
  draw_vof ("cs", "fs", filled = -1, fc = {1, 1, 1});
  draw_vof ("cs", "fs", lw = 2, lc = {0, 0, 0});
  squares ("b", linear = true, min = 0, max = 1);
  save ("bb.mp4");
}

