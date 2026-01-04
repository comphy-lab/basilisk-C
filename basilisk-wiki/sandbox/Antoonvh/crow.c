/**
![Condensation trails of airplanes are often associated with vortical motion and can be subject to the Crow instability between two coutner-rotating vortex tubes](https://www.researchgate.net/profile/Tufan-Kumar-Guha/publication/311977373/figure/fig23/AS:960082900688911@1605913009246/Contrails-showing-the-development-of-Crow-instabilities-in-the-wake-of-an-aircraft.jpg)

# Crow's instability

...Using a 4th order solver and bwatch.

![Volumetric rendering of the $\lambda_2$ vortex-detection field. The visualization shows the evolution of two counter-rotating vortices in a periodic box. These are initialized with a long-wave disturbance (compared to vortex-core size) such that they are subject to the Crow instability mechanism (Crow, 1970)[^1]](https://www.antoonvanhooft.nl/media/crow.mp4)
*/
#define RKORDER 3
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;
#include "filaments.h"
#include "lambda2.h"
#include "bwatch.h"

int n_seg = 1000;
double ue = 2e-4, vis = 5e-4;

int main(int argc, char ** argv) {
  if (argc > 1)
    ue = atof (argv[1]);
  if (argc > 2)
    vis = atof (argv[2]);
  foreach_dimension()
    periodic(left);
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0/2;
  N = 1 << 5;
  a_LO = 0.2;
  const scalar visc[] = vis;
  nu = visc;
  run();
}
/**
The vortex tubes are defined as parametric cuves
 */
double d = 0.8;
coord line1 (double theta) {
  coord C;
  C.x = X0 - L0 + L0*theta;
  C.y = d/2 + a_LO*sin(theta*2*pi)/2. ;
  C.z = 0;
  return C;
}

coord line2 (double theta) {
  coord C;
  C.x = X0 + 2*L0 - L0*theta;
  C.y = -d/2 + a_LO*sin(theta*2*pi)/2. ;
  C.z = 0;
  return C;
}


event init (t = 0) {
  TOLERANCE = 1e-5;
  vector omega[], psi[];
  int n_seg = 1000;
  foreach_dimension()
    psi.x.prolongation = refine_4th;
  foreach() {
    foreach_dimension()
      psi.x[] = 0;
  }
  boundary ((scalar*){psi});
  double oe = 0.02/sqrt(a_LO);
  do {
    get_vor_vector (omega, line1, 0 , 3, n_seg, Lamb_Oseen);
    get_vor_vector (omega, line2, 0 , 3, n_seg, Lamb_Oseen, true);
    printf ("#cells: %ld \n", grid->tn);
  } while (adapt_wavelet ((scalar*){omega}, (double[]){oe, oe, oe}, 8).nf > (grid->tn/100));
  foreach_dimension() {
    stats so = statsf (omega.x);
    foreach()
      omega.x[] -= so.sum/so.volume;
  }
  foreach_dimension()
    poisson (psi.x, omega.x);
  vector uc[];
  foreach_dimension()
    uc.x.prolongation = refine_4th;
  foreach(){
    foreach_dimension()
      uc.x[] = -((8*(psi.z[0,1] - psi.z[0,-1]) + psi.z[0,-2] - psi.z[0,2]) -
		 (8*(psi.y[0,0,1] - psi.y[0,0,-1]) + psi.y[0,0,-2] - psi.y[0,0,2]))/(12*Delta);
    uc.z[] += L0/40; //Co-moving reference frame.
  }
  boundary ((scalar *){uc});
  vector_to_face (uc);
  project (u, p);
}

#define FUNC(x) (exp(-x) + x - 1)

event mov (t += 0.1) {
  static FILE * fp = popen ("ppm2mp4 crow.mp4", "w");
  vector uc[];
  face_to_vector (uc);
  scalar l2[];
  lambda2 (uc, l2);
  foreach() 
    l2[] = l2[] < 0 ? FUNC(-l2[]): 0;
  boundary ({l2});
  watch (fov = 7, O = {1e-6, 5, 20},
	 poi = {1e-4, 1e-4, 1e-4},
	 nx = 1024, ny = (40*24));
  volume (l2, sc = .3, min = -15, max = 15,
	  cols = true, shading = 1, mval = 1e-3);
  store (fp);
  //store (fopen("test.ppm", "w"));
  plain();
}

event adapt (i++)
  adapt_flow (ue, 99, 1);

event logger (i++) {
  fprintf (stderr, "%d %g %d %d %d %d %ld %d\n", i, t, mgp.i, 
           mgp.nrelax, mgp2.i, mgp2.nrelax, grid->tn, grid->maxdepth);
}

event stop (t = 35);
/**
## Reference

[^1]: Crow, S. C. (1970). "Stability theory for a pair of trailing vortices". AIAA Journal. 8 (12): 2172–2179. 


*/
