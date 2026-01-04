/**
# Films do not break

![Interface and horizontal velocity component](df_rain/wider.mp4)(width="50%")

![Zoom on the mesh, interface and horizontal velocity component](df_rain/movie.mp4)(width="50%")
*/

#include "navier-stokes/centered.h"
#if CLSVOF
# include "two-phase-clsvof.h"
#else
# include "two-phase.h"
#endif
#include "reduced.h"
#if CLSVOF
# include "integral.h"
#else
# include "tension.h"
#endif
#include "navier-stokes/perfs.h"
#include "view.h"

int maxlevel = 11;

double tend = 13;
double We = 5.;
double Re = 1000;

double femax = 3e-2, uemax = 5e-2;

int main()
{
  periodic(right);
  
  L0 = 10.;
  X0 = - L0/2.;
  Y0 = - L0/2. - 1.;
#if CLSVOF
  const scalar sigma[] = 1./We;
  d.sigmaf = sigma;
#else
  f.sigma = 1./We;
  G.y = -1.;
#endif
  
  rho1 = 1, rho2 = 0.1;
  mu1 = mu2 = 1./Re;

  TOLERANCE = 1e-3 [*];

  N = 64;
  
  display_control (femax, 1e-2, 0.1);
  display_control (uemax, 1e-2, 0.1);
 
  run();
}

#if CLSVOF
event acceleration (i++)
{
  foreach_face (y)
    a.y[] -= 1.;
}
#endif

event init (t = 0) {
  refine ( sq(y) < 2. && level < 8 );
#if CLSVOF
  foreach()
    d[] = - sq(y + 0.01*sin(2*pi*x/L0*2)) + sq(0.5);
#else
  fraction (f, -sq(y + 0.01*sin(2*pi*x/L0*2) ) + sq(0.5) );  
#endif
}

event movies (i += 10; t <= tend)
{
  view (fov = 5,
	tx = -0.126, ty = -0.012, tz = -1.688,
	width = 1000, height = 1000);
  squares (color = "u.x", spread = -1);
  box ();
  draw_vof (c = "f", lc = {1,0,0});
  cells ();
  save ("movie.mp4");
  
  view (fov = 30, near = 0.01, far = 1000,
	tx = 0.019, ty = 0.095, tz = -1.632,
	width = 1000, height = 1000);
  squares (color = "u.x", spread = -1);
  box ();
  draw_vof (c = "f", lc = {1,0,0});
  save ("wider.mp4");
}

#if TREE
event adapt (i++) {
  vector ut[];
  foreach()
    foreach_dimension()
      ut.x[] = u.x[]*(1. + tanh((y - Y0 - L0/2.)*10.))/2.;
  adapt_wavelet ({f, ut}, (double[]){femax, uemax, uemax}, maxlevel, minlevel = 7);
}
#endif
