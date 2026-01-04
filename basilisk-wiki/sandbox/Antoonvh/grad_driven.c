/**
# Modified tracer-gradient flow

We consider Stokes flow driven by a force field $\vec{a} =
\vec f(\nabla \phi)$. I.e. a (non-linear) vector-valued function of $\nabla
\phi$, where $\phi$ is a tracer field. More specifically,

$$\vec{a} = -|\nabla \phi|\ \nabla \phi.$$

Results:

![Tracer field: After a chaotic initial period we see the emergence of (quasi) stable "blobs"](grad_driven/phi.mp4)

![Vorticity field: These blobs appear to be associated with vorticity](grad_driven/vort.mp4)
*/

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "utils.h"

scalar phi[] , * tracers = {phi};
face vector av[];

const face vector muc[] = {0.001, 0.001};

int j = 0;

int main() {
  periodic (left);
  L0 = 10;
  X0 = Y0 = -L0/2.;
  a = av;
  mu = muc;
  N = 128;
  phi.gradient = minmod2;
  for (j = 0; j < 2; j++) 
    run();
}
/** The CFL condition is relevant for the tracer advection.*/
event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}
/** We study two cases (identified with `j`) */
event init (t = 0) {
  TOLERANCE = 1e-5;
  stokes = true;
  foreach()
    phi[] = (j ? y < -L0/4. : exp(-sq(x) - sq(y))) + noise()/100.;
}
/**
We implement $\vec{a} = |\nabla \phi |\ \nabla \phi = |\nabla \phi|^2 \frac{\nabla \phi}{|\nabla \phi|}$.
 */
event acceleration(i++) {
  vector gr[];
  gradients ({phi},{gr});
  foreach() {
    double r  = 0;
    foreach_dimension()
      r += sq(gr.x[]);
    if (r > 0) {
      foreach_dimension()
	gr.x[] = r*gr.x[]/(sqrt(r));
    } else // 0 = 0
      foreach_dimension()
	gr.x[] = 0;
  }
  foreach_face() 
    av.x[] = -face_value (gr.x, 0);
}

// Diffussion with (Pr = 1)
event tracer_diffusion(i++) {
  diffusion (phi, dt, mu);
}

// Output movies
event mov (t += 0.5) {
  output_ppm (phi, file = "phi.mp4", n = 300, min = -0.75, max = 0.75);
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "vort.mp4", n = 300, min = -5, max = 5, map = cool_warm);
}

event stop (t = 100);
