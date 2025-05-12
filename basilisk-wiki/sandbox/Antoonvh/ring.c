/**
# Two opposing Gaussian vortex rings  

Two Gaussian vortex rings in a triply-periodic domain.

![The flow structure dissipates](ring/l2.mp4)

See [this example](two_rings.c) for a more pleasing result.
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "filaments.h"
#include "lambda2.h"
#include "view.h"

#define R      (sqrt(sq(x) + sq(y)))
#define radi(zp)   (R <= 0.01 ? major_rt				\
		    : (sqrt(sq(x - major_rt*(x/R)) +			\
			    sq(y - major_rt*(y/R)) + sq(z - zp))))

double ue = 5e-3;
int maxlevel = 9;

double major_rt = 5; // Ring' major radius (minor = 1)
double d1 = 6./2.;    // Initial distance
double muv = 5e-4;   // Fluid's viscosity
int main() {
  foreach_dimension()
    periodic(left);
  L0 = 125;
  origin (-L0/2 + pi, -L0/2 + exp(1), -L0/2 - sqrt(2));
  const face vector muc[] = {muv, muv, muv};
  mu = muc;
  run();
}
/**
The rings need to be defined as parametric curves
 */
coord ring1 (double t) {
  coord C;
  C.x = major_rt*cos(t);
  C.y = major_rt*sin(t);
  C.z = d1;
  return C;
}

coord ring2 (double t) {
  coord C;
  C.x =  major_rt*cos(t);
  C.y = -major_rt*sin(t);
  C.z = -d1;
  return C;
}

event init (t = 0) {
  refine ((radi(-d1) < 8 || radi(d1) < 8) && level < maxlevel - 2);
  refine ((radi(-d1) < 4 || radi(d1) < 4) && level < maxlevel - 1);
  refine ((radi(-d1) < 2 || radi(d1) < 2) && level < maxlevel);
  /**
## Initializing a periodic flow with a prescribed vorticity distribution

The vector potential $\mathbf{A}$ is related to the flow $\mathbf{u}$,

$$\mathbf{u} = \mathbf{\nabla} \times \vec{A}.$$

The vector potential is also related to the vorticity vector $\vec{\omega}$,

$$\mathbf{\nabla}^2 \vec{A} = \vec{\omega},$$

where $\mathbf{\nabla}^2$ is a component wise Laplacian operator.
   */
  vector A[], omg[], omg2[];
  int segs = 500;
  get_vor_vector (omg,  ring1, 0, 2*pi, segs, Lamb_Oseen);
  get_vor_vector (omg2, ring2, 0, 2*pi, segs, Lamb_Oseen);
  foreach() {
    foreach_dimension() {
      A.x[] = 0.;
      omg.x[] += omg2.x[];
    }
  }
  foreach_dimension() {
    stats o = statsf (omg.x);
    foreach() {
      omg.x[] -= o.sum/o.volume;
    }
  }
  boundary ((scalar*){A, omg});
  foreach_dimension()
    poisson (A.x, omg.x);
  foreach() { 
    foreach_dimension() {
      u.x[] = (A.z[0,1,0] - A.z[0,-1,0] -
	       A.y[0,0,1] + A.y[0,0,-1])/(2*Delta);
    }
  }
  boundary ((scalar*){u});
}

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){ue, ue, ue}, maxlevel);

event mov (i += 4) {
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 12, phi = 0.5, theta = 0.4);
  isosurface ("l2", -0.0001);
  cells (n = {0, 1, 0});
  char str[99];
  sprintf (str, " t = %.3g", t);
  draw_string (str, size = 20);
  save ("l2.mp4");
}

event stop (t = 1000);
