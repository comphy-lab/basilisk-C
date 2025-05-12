/**
# Internal Waves 

In a stratified fluid internal waves can exist (also reffered to as
gravity waves). An interesting feature of these waves is the so-called
dispersion relation between the angle of wave propagation ($\theta$),
stratification strength ($N^2$) and the freqency of the wave
($\omega$), according to,

$$ \omega = N^2 \cos(\theta).$$

## Numerical set-up

The Navier-Stokes equantions under Boussinesq approximation are solved
on a $256 \times 256$ miltigrid. In the centre of the domain an
oscillating force exites the internal waves with a freqency
corresponding to $\theta = 45^o$.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"

scalar b[];
scalar * tracers = {b};
face vector av[];
double sqN = 1., omega;

b[top]    = neumann (sqN);
b[bottom] = neumann (-sqN);

int main() {
  omega = sqrt(1./2.);
  L0 = 30;
  X0 = Y0 = -L0/2;
  a = av;
  TOLERANCE = 1e-4;
  DT = 0.2/omega;
  N = 256;
  run();
}
/**
   The initial stratification is set.
*/
event init (t = 0) {
  foreach()
    b[] = sqN*y;
  boundary ({b});
}

/**
## Acceleration

We apply gravity and the localized oscillarory force.
 */
event acceleration (i++) {
  coord del = {0, 1};
  foreach_face() 
    av.x[] = del.x*((b[] + b[-1])/2. +
		    0.1*(sin(omega*t)*((sq(x) + sq(y)) < 1)));
}
/**
## Output 

We output a .mp4 file showing the evolution of the magnitude of the
gradient of the buoyancy field (|$\nabla b$|).
*/
event output (t += 0.5; t <= 75) {
  scalar grb[];
  foreach() {
    grb[] = 0;
    foreach_dimension()
      grb[] += sq((b[1] - b[-1])/(2*Delta));
    grb[] = sqrt(grb[]);
  }
  output_ppm (grb, file = "grb.mp4", min = 0.8, max = 1.2);
}
/**
## Results

The dispersion relation appears to be statisfied.

![Visualization of the internal waves](internalwacesMG/grb.mp4)

The next step is to perform this simulation using adaptive grids, [See
here](internalwacesAMR.c).
*/
