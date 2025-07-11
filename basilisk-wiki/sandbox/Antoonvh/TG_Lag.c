
/**
# Taylor-Green vortex test

![Works pretty good. Well done Lagrangian advection!](http://antoonvanhooft.nl/media/TG_Lag.mp4)

![Convergence in OK!](http://antoonvanhooft.nl/media/conv.png)
)
 */

#define RKORDER 3

#include "view.h"
#define VIEW 1
#include "master-omgpsi_LAG.h"
#include "higher-order.h"


int maxlevel = 6;
int main() {
  foreach_dimension() {
    periodic(left);
    periodic_x = true;
  }
  L0 = 1.;
  X0 = Y0 = -0.5;
  
  for (maxlevel = 4; maxlevel <= 7; maxlevel++) {
    N = 1 << (maxlevel);
    DT = 2*L0/N;
    run();
  }
}

event init (t = 0) {
  slave_init (L0, (coord){X0, Y0}, maxlevel, false);
  p_omg = new_particles(0);
  foreach() {
    foreach_child() {
      particle pn;
      pn.x = x;
      pn.y = y;
      pn.s = 2.*sin(2.*pi*x)*sin(2.*pi*y);
      pn.ti = t;
      add_particle(pn, p_omg);
    }
  }
}

#if 1
event movie (t += 0.1) {
  if (N == 64) {
    view (width = 800, height = 800, bg = {0,0.5, 0.5});
    scatter_color (p_omg, s = 2, map = blue_white_red);
    save ("TG.mp4");
    save ("TG.png");
    cells();
    squares ("omega", min = -1, max = 1 );
    save ("TG-grid.mp4");
    slave_level();
  }
}
#endif

event adapt (i++) {
  compute_omega(p_omg, omega);
  adapt_number();
}

event stop (t = 2) {
  double err = 0, me = 0;
  int np = 0;
  foreach_particle_in(p_omg, reduction(+:err) reduction(+:np) reduction(max:me)) {
    if (p().ti == 0) {
      double el = fabs(p().s - 2*sin(2*pi*x)*sin(2*pi*y)); 
      err += el;
      if (el > me)
	me = el;
      np++;
    }
  }
  printf ("%d %g %g\n", N, err/np, me);
  return 1;
}
