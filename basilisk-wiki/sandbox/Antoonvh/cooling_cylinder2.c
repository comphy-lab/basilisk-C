/**
# Cooling cylinder

We study the convergence of the heat in a cooling cylinder with increasing resolution.

![Evolution of the Buoyancy/temperature at various resolutions](cooling_cylinder2/b.mp4)

~~~gnuplot
set size square
set grid
set xlabel 't'
set ylabel 'heat in cylinder'
plot 'out' palette
~~~
 */
#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar b[], * tracers = {b};

vertex scalar dist[];

double muv = 0.01;
double kappa = .1, kappaf = 5.;

face vector muc[], kap[], av[];

u.t[embed] = dirichlet (0);
u.n[embed] = dirichlet (0);

int main() {
  L0 = 10;
  X0 = -L0/2.;
  Y0 = -L0/4.;
  mu = muc;
  a = av;
  DT = 2e-1;
  for (N = 32; N <= 1025; N *= 2)
    run();
}


event init (t = 0) {
  TOLERANCE = 1e-4;
  foreach_vertex()
    dist[] = sqrt(sq(x) + sq(y));
  fractions (dist, cs, fs, 1);
  double DELT = L0/N;
  int w = 5;
  /**
The diffusivity jump is smeared out over a transition region that is `w` cells wide.
   */
  foreach_face() {
    double d = (dist[] + dist[0,1])/2. ;
    double v = (kappa - kappa/kappaf)/(2*w*DELT)*(1 - d) + (kappa + kappa/kappaf)/2;
    kap.x[] = max(min(v, kappa), kappa/kappaf);
  }
  foreach_face()
    muc.x[] = fm.x[]*muv;
  double b0  = 1;
  foreach() {
    double d = (dist[] + dist[0,1] + dist[1,1] + dist[1])/4.;
    b[] = b0*erf((1 - d)/0.25)/2;
  }
}

event diffution (i++) {
  diffusion (b, dt, kap);
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event movies (t += 0.1) {
  output_ppm (b, file = "b.mp4", n = 300, max = 0.6, min = -0.6);
}

event log (i++) {
  double tb = 0;
  foreach(reduction(+:tb)) {
    tb += sq(Delta)*(1 - cm[])*b[];
  }
  if (pid() == 0)
    printf ("%g %g %d\n", t, tb, N);
}

event stop (t = 10);
