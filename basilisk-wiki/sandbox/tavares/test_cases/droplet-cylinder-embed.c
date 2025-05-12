/**
# droplet on an embedded cylinder


~~~gnuplot Equilibrium shapes for $30^\circ \leq \theta \leq 150^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1:1]
set yrange [0:]
f0(x) = sqrt(0.5**2 - x**2) + 0.575
f(x)  = sqrt(0.5**2 - x**2) - 0.575

plot 'out' w l lt -1 lw 3 lc rgb "blue" t 'Numerical solution',\
  f0(x) with filledcurves above y1 = 0.575 fc "black" t 'cylinder',\
 -1*f(x) with filledcurves below y1 = 0.575 fc "black" t ''
 
set term pop
~~~
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tavares/contact-embed.h"
#include "view.h"

#define R0 0.5
#define xc 0.
#define yc 0.575 
#define T 15

double theta0, volume_vof_init;
int LEVEL = 6;

int main() {
  size (2.);

  /**
  We set the origin */

  origin (-1, 0);

  init_grid (1 << LEVEL);
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
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

 //u.n[embed] = dirichlet( 0 ); 
 //u.t[embed] = dirichlet( 0 );

event init (t = 0)
{
  /**
  We define the cylinder and the initial (half)-circular
  interface. */
  solid (cs, fs, (sq(x - xc) + sq(y - yc) - sq(R0)));

  fraction (f, - (sq(x - xc) + sq(y - (yc+sqrt(2)/2)) - sq(R0)));
}

event logfile (i++; t <= T)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;
}

#if 0
event volume (i++, t<=T){
  if (t==0) volume_vof_init = statsf (f).sum;

  char name[80];
  sprintf (name, "volume-mesh%d-angle%g.dat", N, theta0);
  static FILE * fp = fopen (name,"w");
  stats s = statsf (f);
  double erreur = ((volume_vof_init - s.sum)/volume_vof_init)*100;
  fprintf (fp, "%g %.5g %.5g\n", t, erreur, dt); 
}
#endif


event end (t = end)
{
  /**
  At the end, we output the equilibrium shape. */

  output_facets (f, stdout);

 /**
  We compute the curvature only in full cells. */
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.5g %.3g\n", N, theta0, R, R/sqrt(V/pi), s.stddev);
}

event movie(i+=10,last){
  if (theta0 == 120) {
    view(fov=20, tx = 0, ty = -0.5);
  //view(fov=20, quat = {0,0,-0.707,0.707}, tx = 0, ty = -0.5);
    draw_vof ("f", lw=2);
    draw_vof("cs", "fs",filled=-1);
    squares("f", linear = true, min = 0, max = 1);
    save("movie.mp4");
  }
}
/**
![Relaxation toward a $120^\circ$ contact angle.](droplet-cylinder-embed/movie.mp4)
*/

// fixme: Comparison to theory is missing (add soon) 