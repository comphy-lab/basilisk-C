/**
# 3D Sessile drop

This is the 3D equivalent of the [2D](sessile_embed.c) test case.

The volume of a [spherical
cap](https://en.wikipedia.org/wiki/Spherical_cap) of radius $R$ and
(contact) angle $\theta$ is
$$
V = \frac{\pi}{3}R^3(2+\cos\theta)(1-\cos\theta)^2
$$
or equivalently
$$
\frac{R}{R_0} = \left(\frac{1}{4}(2+\cos\theta)(1-\cos\theta)^2\right)^{-1/3}
$$
with $R_0$ the equivalent radius of the droplet
$$
R_0 = \left(\frac{3V}{4\pi}\right)^{1/3}
$$
To test this relation, a drop is initialised as a half-sphere
(i.e. the initial contact angle is 90$^\circ$) and the contact angle
is varied between 30$^\circ$ and 150$^\circ$. The drop oscillates and
eventually relaxes to its equilibrium position. The curvature along
the interface is close to constant.

![Relaxation toward a $60^\circ$ contact angle.](sessile-embed3D/movie.mp4)

 */

#include "grid/multigrid3D.h"
//#include "grid/octree.h"
#include "../embed-navier.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact-embed.h"
#include "view.h"

#define T 1.

double theta0=60, volume_vof_init;

#define MAXLEVEL 5

int main()
{
  size (1.);

  /**
  We shift the bottom boundary. */

  origin (0, -0.26, 0);
  //origin (-0.5, -0.26, -0.5);

  N = 1 << MAXLEVEL;
  
  // We use a constant viscosity. 
  mu1 = 0.1; mu2 = 0.1;
  rho1= 0.1; rho2 = 0.1;
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */
  for (theta0 = 30; theta0 <= 150; theta0 += 30) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}

event init (t = 0)
{
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
  phi[] = y;
  boundary ({phi});
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(0.25)));
}


event logfile (i++; t <= T)
{
  /**
  If the standard deviation of curvature falls below $10^{-2}$, we stop the computation
  (convergence has been reached). */  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-2)
    return true;
}


#if 1
event cleanfghost (i++, t<=T, last) { // clean fluid ghost cells
  foreach()
    if (cs[]<=0) f[]=0.;
}
#endif

#if 1
event movie (i += 5; t <= T)
{
 if (theta0 == 60) {
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 32, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.299, tz = -2.748, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
  mirror (n = {1,0,0}) {
  cells (n = {0,1,0});
  draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
    }
  mirror (n = {0,0,1}) {
  cells (n = {0,1,0});
  draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
    mirror (n = {1,0,0}) {
    cells (n = {0,1,0});
    draw_vof (c = "cs", s = "fs", filled = -1, fc = {0.8,0.8,0.8});
    draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
    draw_vof (c= "f", edges = true);
    }
  }
  save ("movie.mp4");
  }
}
#endif


/** We check the volume conservation */
#if 0
event volume (i++, t<=T){
  if (t==0) volume_vof_init = statsf (f).sum;
  char name[80];
  sprintf (name, "volume-angle%g-mesh%d.dat", theta0, N);
  static FILE * fp = fopen (name,"w");
  stats s = statsf (f);
  double erreur = ((volume_vof_init - s.sum)/volume_vof_init)*100;
  fprintf (fp, "%g %.5g\n", t, erreur); 
}
#endif



event end (t = end)
{
  /** At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/

  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);
}


#if TREE
event adapt (i++) {
#if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif


/**
We compare $R/R_0$ to the analytical expression.

The accuracy is not as good (yet) as that of the [sessile
test](/src/test/sessile3D.c).

~~~gnuplot


reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 30,1 to 150,1 nohead dt 2
kappa(theta) = 2.*((2. + cos(theta))*(1. - cos(theta))**2/4.)**(1./3.)
R0(V) = (3.*V/(4.*pi))**(1./3.)
set xtics 30,30,150
plot 2./(kappa(x*pi/180.)) lw 3 lt -1 dt 2 lc rgb "black" t 'analytical', \
     'log' u 2:(2.*$3/R0($5)) w p pt 4 ps 3 lt -1 lw 3 t 'numerical'

~~~

## See also

* [3D Sessile drop on a domain boundary](/src/test/sessile3D.c)

*/



