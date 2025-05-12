/**
# Initializing a vortex ring in 3D

On this page we check if it is possible to initialize coherent vortex
ring without resorting to complex models. Or even doing proper
solenoidal initializaiton. For this purpose, we generate a compact
ring vortex by projecting the azimuthally rotated Kaufmann vortex onto
the space of divergence free vector fields.

![Quite a steady structure](compactring/ring.mp4)

We compare the translation velocity against a velocity scale:

$$V = \frac{\Gamma}{R_{major}},$$

~~~gnuplot The translation velocity
  set xlabel 'time'
  set ylabel 'Position'
  set key top right
  plot 'out', x*1.4 t '1.4*Vt'
~~~

and a profile,

~~~gnuplot
reset
set parametric
set tr [0:40]
c = 0
set cbr [-10:70]
set yr [0:40]
set xlabel 'u.z'
set ylabel 'r'
set grid
plot 'log' u 3:2:1 w l palette lw 2  t 'time', c,t lw 2 t ''
~~~

For more characterization see [here](vrk.c). But the axisymmetric setup does not seem to caputure the setup here.
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "lambda2.h"

#define RANKINE (d > R1 ? Gamma1/(2*pi*d) : d*Gamma1/(2*pi*sq(R1)))
#define KAUFMANN (d*Gamma1/((sq(R1) + sq(d))*2*pi)) 

/**
   Define some toroidal coordinates. 
*/
double major_rt = 10;
#define R      (sqrt(sq(x) + sq(y)))
#define SINPSI (y/R)
#define COSPSI (x/R)
#define radi   (R <= 0.01 ? major_rt				\
		: (sqrt(sq(x - major_rt*COSPSI) +		\
			sq(y - major_rt*SINPSI) + sq(z))))

/**
The domain is centered around {0,0,0}.
 */
int maxlevel = 9;
double R1, Gamma;
int main() {
  R1 = major_rt/6.;       // Compactness
  Gamma = major_rt*2.*pi; // V = 1
  foreach_dimension()
    u.x.refine = refine_linear;
  periodic (back);
  L0 = 120;
  X0 = Y0 = -L0/2;
  Z0 = -L0/4.;
  run();
}

/**
## Initialize the flow
*/
coord cross (coord a, coord b) {
  coord product;
  foreach_dimension() 
    product.x = a.y*b.z - a.z*b.y;
  return product;
}

coord induced_v (double Gamma1, double R1, Point point, coord pr) { 
  coord u, c = {x, y ,z};
  double d = 0;
  foreach_dimension()
    d += sq(c.x - pr.x);
  d = sqrt(d);
  double vtheta = KAUFMANN;
  coord rc_to_pc = {x - pr.x, y - pr.y, z - pr.z}; 
  coord tang = {-pr.y/d, pr.x/d, 0};
  coord n = cross (tang, rc_to_pc);
  normalize (&n);
  foreach_dimension()
    u.x = vtheta*n.x;
  return u;
}

event init (t = 0) {
  refine (radi < 3*R1 && level < maxlevel - 1);
  refine (radi < 1.5*R1 && level < maxlevel);
  foreach() {
    coord pr1 = {major_rt*COSPSI, major_rt*SINPSI, 0}; // Closest vortex center
    coord u1 = induced_v (Gamma, R1, point, pr1);
    foreach_dimension()
      pr1.x *= -1;
    coord u2 = induced_v (Gamma, R1, point, pr1);
    foreach_dimension()
      u.x[] += u1.x + u2.x;
  }
  boundary ((scalar*){u});
}

#if ADD
event new_vortex (t = ADD)
  event ("init"); //Leaf frogging
#endif

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.05, 0.05, 0.05}, maxlevel);

event position (t += 1) {
  double omgz = 0, omgt = 0, D;
  scalar l2[];
  lambda2 (u, l2);
  foreach(reduction(+:omgz) reduction(+:omgt)) {
    if (l2[] < 0.1) {
      omgz += dv()*fabs(l2[])*z;
      omgt += dv()*fabs(l2[]);
    }
  }
  printf ("%g %g\n", t, omgz /= omgt);
  fflush (stdout);
  if ((int)(t + 0.5) % 10 == 0) {
    for (double xp = X0 + 0.1; xp < X0 + L0 - 0.1; xp += D) {
      Point p = locate (xp, 0, omgz);
      if (p.level > 0) {
	D = Delta;
	fprintf (stderr, "%g %g %g\n", t, xp, interpolate (u.z, xp, 0, omgz));
      } else
	D = HUGE;
    }
  fputs ("\n", stderr);
  }
}

event bviewer (t += 0.25) {
  view (theta = 0.6, phi = 0.35);
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});
  translate (z = Z0) {
    box();
    isosurface ("l2", -0.1);
    translate (y = -L0/2)
      cells (n = {0,1,0});
  }
#if ADD
  save ("rings.mp4");
#else
  save ("ring.mp4");
#endif
}

event end (t = 60);
