#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "no-coalescence.h"

#include "view.h"
#include "profiling.h"

f[top] = 0.;
f[bottom] = 0.;
f[left] = 0.;
f[right] = 0.;

u.t[left]  = dirichlet(0);
u.t[right] = dirichlet(0);

int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  mu1 = mu2 = 0.01;
  rho1 = 0.1, rho2 = 1.;
  f.sigma = 0.1;
  N = 128;
  G.y = -1;
  run();
}

event init (t = 0)
{
  int nb = 7;
  double a = 2.*L0/N, R = 0.9999*(L0 - (nb + 1)*a)/(2.*nb);
  for (double xc = -L0/2. + R + a; xc <= L0/2. - (R + a); xc += 2.*R + a)
    for (double yc = -L0/2. + R + a; yc <= 0.; yc += 2.*R + a) {
      scalar f1[];
      fraction (f1, - (sq(x - xc) + sq(y - yc) - sq(R)));
      foreach()
	f[] += f1[];
    }
  boundary ({f});
}

event logfile (i++)
  fprintf (stderr, "%d %g %g\n", i, t, statsf(f).sum);

event movie (t += 0.1)
{  
  clear();
  double cmap[NCMAP][3];
  randomap (cmap);
  for (scalar c in interfaces) {
    color b = colormap_color (cmap, c.i, 0, 100);
    draw_vof (c.name, filled = 1, lc = {b.r/255., b.g/255., b.b/255.});
  }
  squares ("u.y", linear = true);
  box();
  save ("movie.mp4");
}

event end (t = 35)
  dump ("dump");
