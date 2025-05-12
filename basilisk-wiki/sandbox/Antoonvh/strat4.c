/**
# Interal waves and the dispersion relation

Internal waves can exist in a stratified fluid. An interesting feature
of these waves is the so-called dispersion relation between the angle
of wave propagation ($\theta$), stratification strength ($N^2$) and
the freqency of the wave ($\omega$), according to,

$$ \omega = N^2 \cos(\theta).$$

The Navier-Stokes equations under the Boussinesq approximation are
 solved. In the centre of the domain an oscillating force exites the
 internal waves with a freqency corresponding to $\theta = 45^o$.

![A pass for the 4th order solver](strat4/grb.mp4)
*/
#include "nsf4t.h"
#include "view.h"
scalar b[], zeros[], * tracers = {b};
b[bottom] = dirichlet_vert_bottom; 
b[top] = dirichlet_vert_top (L0);
vector av;

double omega;

int main() {
  omega = sqrt(2)/2.;
  periodic (left);
  av.x = zeros;
  av.y = b;
  L0 = 30;
  a  = av;
  TOLERANCE = 1e-5;
  DT = 0.2;
  N = 128;
  run();
}

event init (t = 0) {
  foreach_vert()
    b[] = y;
  boundary ({b});
}

event accel (i++) {
  foreach()
    if (sq(x - L0/2) + sq(y - L0/2) < 1)
      u.y[] += 0.1*dt*sin(omega*t); // Not accurate, but flexible...  
}

event output (t += 0.5) {
  scalar grb[];
  foreach() {
    grb[] = 0;
    foreach_dimension()
      grb[] += sq((b[1] - b[-1])/(2*Delta));
    grb[] = sqrt(grb[]);
  }
  output_ppm (grb, file = "grb.mp4", min = 0.8, max = 1.2, n = 300);
}

event stop (t = 75);
