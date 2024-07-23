/**
# Stress test for wetting and drying

Bumps uniformly covered with water create ponds and streams (i.e. a
lot of wetting/drying going on).

~~~gnuplot Free surface coloured with the norm of the velocity.
set term @PNG enhanced size 1024,512 font ",8"
set output 'eta.png'

! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out

unset key
set hidden3d
unset xtics
unset ytics
unset border
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

dry=1e-3
set view 28,55

set multiplot layout 2,2 scale 1,1.3
splot './eta-1' u 1:2:4 every 3:3 w l, \
      './eta-1' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-2' u 1:2:4 every 3:3 w l, \
      './eta-2' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-3' u 1:2:4 every 3:3 w l, \
      './eta-3' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-8' u 1:2:4 every 3:3 w l, \
      './eta-8' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
unset multiplot
~~~

~~~gnuplot Evolution of the level of refinement.
set output 'level.png'
dry = 0.
set multiplot layout 2,2 scale 1,1.3
splot './level-1' u 1:2:4 every 3:3 w l, \
      './level-1' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-2' u 1:2:4 every 3:3 w l, \
      './level-2' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-3' u 1:2:4 every 3:3 w l, \
      './level-3' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-8' u 1:2:4 every 3:3 w l, \
      './level-8' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
unset multiplot
~~~

~~~gnuplot Velocity field at the end of the simulation.
reset
set term @PNG enhanced size 1024,512 font ",8"
set output 'vectors.png'

set xrange [0:1000]
set yrange [0:1000]
set size ratio -1
unset key
dry=1e-3
plot './eta-8' u 1:2:($3>dry?$5*30.:0):($3>dry?$6*30.:0) every 2:2 w vec lc 0

! rm -f eta-? level-?
~~~
*/

// use KDT database rather than analytical function
#define KDT 1

#if ML
# include "layered/hydro.h"
# include "layered/nh.h"
#elif IMPLICIT
# include "saint-venant-implicit.h"
#else
# include "saint-venant.h"
#endif

#if KDT
# include "terrain.h"
#endif

#define LEVEL 7

int main()
{
  init_grid (1 << LEVEL);
  size (1000.);
  G = 9.81; 
#if IMPLICIT || ML
  DT = 10.;
# if ML  
  CFL_H = HUGE;
# endif
#endif
  run();
}

#define zb(x,y) ((cos(pi*x/L0)*cos(pi*y/L0) + \
		  cos(3.*pi*x/L0)*cos(3.*pi*y/L0)) - 2.*x/1000.)

#if !KDT
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = zb(x,y);
}
#endif

event init (i = 0)
{
#if !KDT
  zb.refine = refine_zb; // updates terrain
  foreach() {
    zb[] = zb(x,y);
    h[] = 0.1;
  }
#else
  FILE * fp = popen ("xyz2kdt ponds", "w");
  for (double x = 0.; x <= 1000; x += 1.)
    for (double y = 0.; y <= 1000.; y += 1.)
      fprintf (fp, "%g %g %g\n", x, y, zb(x,y));
  pclose (fp);
  terrain (zb, "ponds", NULL);
  foreach()
    h[] = 0.1;
#endif
}

event friction (i++) {
#if IMPLICIT && !ML
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(q)/sq(h[]);
    foreach_dimension()
      q.x[] /= a;
  }
#else
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
#endif
}

event logfile (i += 10) {
  stats s = statsf (h);
#if IMPLICIT && !ML
  scalar u[];
  foreach()
    u[] = h[] > dry ? q.x[]/h[] : 0.;
  norm n = normf (u);
#else
  norm n = normf (u.x);
#endif
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %.4g %g %g %g %g %g\n", 
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt);
#if !IMPLICIT
  assert (s.min > 0.);
#endif
}

event outputfile (t <= 1200.; t += 1200./8) {
#if ML || !IMPLICIT 
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h, zb, u}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({h, zb, l}, stdout, N);

  nf++;
#endif
}

// int event (t += 100)
//  output_matrix (h, stdout, N, true);

event adapt (i++) {
  // we do this so that wavelets use the default bilinear
  // interpolation this is less noisy than the linear + gradient
  // limiters used in Saint-Venant not sure whether this is better
  // though.
  scalar h1[];
  foreach()
    h1[] = h[];
  adapt_wavelet ({h1}, (double[]){1e-2}, LEVEL, 4); 
}
