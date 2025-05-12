/**
# Stokes flow past a Corona Cell

![Local cell curcature, tracer particles, cells, and flow
 speed](cor/cor.mp4)

The [STL file](https://www.thingiverse.com/thing:4166787) by [Greg
Bejtlich](https://www.thingiverse.com/gregbejtlich/about) via
[Thinkiverse](www.thingiverse.com) was simplified with
[meshlab](http://www.meshlab.net/).
 */
#include "grid/octree.h"
#include "embed.h"
#include "fractions.h"
#include "distance.h"
#include "curvature.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "tracer-particles.h"
#include "scatter2.h"

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
u.r[embed] = dirichlet (0.);

int maxlevel  = 7;
face vector muc[];
Particles P1;

int main() {
  mu = muc;
  periodic (back);
  const face vector av[] = {0, 0, 0.1};
  a = av;
  run();
}

event init (t = 0) {
  stokes = true;
  coord * p = input_stl (fopen ("cor.stl", "r"));
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  size (2.*maxl);
  X0 = (max.x + min.x)/2. - L0/2;
  Z0 = (max.z + min.z)/2. - L0/2;
  Y0 = (max.y + min.y)/2. - L0/2;
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1e-4*L0}, maxlevel).nf);
  boundary ({d});
  vertex scalar phi[];  foreach_vertex()
    phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
	      d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  
  P1 = init_tp_square (15, l = L0/2, zp = Z0 + L0/10.);
  
  DT = 0.2;
}

event defaults (i++) {
  foreach_face()
    muc.x[] = 50*fs.x[];
  boundary ((scalar*){muc});
}

event adapt (i++)
  adapt_wavelet ({cs, u}, (double[]){0.01, 0.2, 0.2, 0.2}, maxlevel);

event movie (t += 0.5) {
  scalar kappa[];
  curvature (cs, kappa);
  view (theta = 0.5, phi = 0.5);
  draw_vof("cs", "fs", color = "kappa");
  box();
  translate (x = -L0/2)
    cells (n = {1,0,0});
  translate (y = -L0/2)
    squares ("u.z", n = {0, 1, 0}); 
  scatter (P1);
  save ("cor.mp4");
}

event stop (t = 150);
