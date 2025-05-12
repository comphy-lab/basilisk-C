/**
![With a little help from man kind, helical vortices exist in
 nature. Image via [FYFD](https://fuckyeahfluiddynamics.tumblr.com/)
 on
 [tumblr](https://www.tumblr.com/tagged/wind-turbine)](https://66.media.tumblr.com/5ea122d9202fdde557494ad87e4f2a2a/tumblr_n8b8lhNz6T1qckzoqo1_400.jpg)

# A helical vortex instability

![The movie shows a $\lambda_2$ isosurface, 100 flow tracers and a
sclice of the grid.](helical/helix.mp4)
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"

#include "view.h"
#include "tracer-particles.h"
#include "scatter2.h"
#include "lambda2.h"
/**
The vortex parameters include:
 */
int not = 3;                   // Number of windings in the domain
double Ls = 5;                 // That domain size
double major_rt = 1, R1 = 0.2; // Helical vortex radii
double Gamma1 = 10;            // Vortex strength

int maxlevel = 7 ;

Particles trac;

int main() {
  L0 = Ls;
  X0 = Y0 = -L0/2;
  periodic (back);
  run();
}
/**
## Helical vortex initialization.

Our helical vortex tube lies along the parameterized line 

$$\mathbf{x}(t) = \{R_m \mathrm{cos}(2\pi t),\ R_m\mathrm{sin}(2\pi t),\ \frac{tL_0}{\mathtt{not}}\},$$

~~~gnuplot
set parametric
L0 = 5
not = 3
set ur [0:not]
set view equal xyz
set view 345, 330
set view azimuth 30
splot cos(2*3.1415*u),sin(2*3.1415*u),u*L0/not
~~~

We consider the velocity induced by the line-vortex elements in the
plane with the cell-centre $\{x, y, z\}$:

$$t \in \{ t_c,\ t_c - 0.5,\ t_c + 0.5,\ t_c - 1,\ t_c - 1.5, ...\}.$$

$$\mathbf{x}_1(t_c) = \frac{xR_m}{r},$$
$$r = \sqrt{x^2 + y^2},$$
$$t_c = \frac{\mathrm{cos}^{-1}(\mathbf{x}_1/r)}{2 \pi}.$$

The final step is for the solver to project such divergent flow field.
*/
#define RAD (sqrt(sq(x) + sq(y)))
#define KAUFMANN (d*Gamma1/((sq(R1) + sq(d))*2*pi))

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
  coord tang = {-pr.y/d, pr.x/d, 0}; //Only approx.
  coord n = cross (tang, rc_to_pc);
  normalize (&n);
  foreach_dimension()
    u.x = vtheta*n.x;
  return u;
}

event init (t = 0) {
  TOLERANCE = 1e-4;
  int tr = 3; // number of ghost turns...
  refine  (RAD < major_rt + 3*R1   && RAD > major_rt - 3*R1   && level < maxlevel - 1);
  refine  (RAD < major_rt + 1.5*R1 && RAD > major_rt - 1.5*R1 && level < maxlevel);
  double Dr = L0/not;
  foreach() {
    double a = major_rt*x/RAD;
    double tc = acos(a)/(2*pi);
    if (y < 0)
      tc = -tc;
    for (double ti = tc - tr; ti <= tr + not ; ti += 0.5) {
      coord pr = {major_rt*cos(ti*2*pi), major_rt*sin(ti*2*pi), ti*Dr};
      coord un = induced_v (Gamma1, R1, point, pr);
      foreach_dimension()
	u.x[] += un.x;
    }
  }
  boundary ((scalar*){u});
  
  trac = init_tp_circle (100, 0, 0, major_rt, L0/2); 
}

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.2, 0.2, 0.2}, maxlevel); 
  
event movie (t += 0.01) { 
  scalar l2[]; 
  lambda2 (u, l2); 
  foreach() 
    l2[] = l2[] > -0.01 ? nodata : l2[];
  stats f = statsf (l2);
  foreach()
    l2[] = l2[] == nodata ? 0 : l2[];
  boundary ({l2});
  view (theta = 0.4, phi = 0.6);
  translate (z = -L0/2) {
    cells();
    isosurface ("l2", -f.stddev);
    scatter (trac);
  }
  save ("helix.mp4");
}

event stop (t = 5) 
  dump();

