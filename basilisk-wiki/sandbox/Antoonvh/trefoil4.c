/**
# Andres Castillo's Trefoil example...

...Using a 4th order solver and bwatch. See [Andres'
page](../acastillo/filaments/trefoil.c) for my inspiration.

<div class="figure">
<video controls="" preload="metadata" width="900">
<source src="https://surfdrive.surf.nl/files/index.php/s/vAfPbpVqwy3dfrn/download
" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Volumetric rendering of the negative $\lambda_2$ field (via surfdrive)
</p>
</div>

<div class="figure">
<video controls="" preload="metadata" width="900">
<source src="https://surfdrive.surf.nl/files/index.php/s/gVqk1Aw0yIPaCpg/download
" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Another volumetric rendering of the negative $\lambda_2$ field (via surfdrive)
</p>
</div>
*/
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;
#include "filaments.h"
#include "lambda2.h"
#include "bwatch.h"

int n_seg = 1000;
double ue = 1e-3, vis = 5e-4;

int main(int argc, char ** argv) {
  if (argc > 1)
    ue = atof (argv[1]);
  if (argc > 2)
    vis = atof (argv[2]);
  foreach_dimension()
    periodic(left);
  L0 = 16;
  X0 = Y0 = Z0 = -L0/2;
  N = 1 << 5;
  a_LO = 0.2;
  const scalar visc[] = vis;
  nu = visc;
  run();
}
/**
The Knot is defined as a parametric curve
 */
coord knot (double theta) {
  coord C;
  C.x = sin(theta) + 2.*sin(2.*theta);
  C.y = cos(theta) - 2.*cos(2.*theta);
  C.z = -sin(3.*theta) - 4;
  return C;
}

event init (t = 0) {
  refine (sq(x) + sq(y) + sq(z) < sq(4) && level < 6);
  vector omega[], psi[];
  int n_seg = 1000;
  foreach_dimension()
    psi.x.prolongation = refine_4th;
  foreach() {
    foreach_dimension()
      psi.x[] = 0;
  }
  boundary ((scalar*){psi});
  double oe = 0.01/sqrt(a_LO);
  do {
    get_vor_vector (omega, knot, 0 , 2*pi, n_seg, Lamb_Oseen);
    printf ("#cells: %ld \n", grid->tn);
  } while (adapt_wavelet ((scalar*){omega}, (double[]){oe, oe, oe}, 9).nf > (grid->tn/100));
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
  }
  boundary ((scalar *){uc});
  vector_to_face (uc);
  project (u, p);
}

#define FUNC(x) (exp(-x) + x - 1)

event mov (t += 0.1) {
  static FILE * fp = popen ("ppm2mp4 tref.mp4", "w");
  vector uc[];
  face_to_vector (uc);
  scalar l2[];
  lambda2 (uc, l2);
  foreach() 
    l2[] = l2[] < 0 ? FUNC(-l2[]): 0;
  boundary ({l2});
  watch (fov = 12, O = {18*sin(sin(t/50)), 5*sin(t/50), 20*cos(t/50)},
	 poi = {0, 0, 2},
	 nx = 1024, ny = (40*24));
  volume (l2, sc = .4, min = -15, max = 15,
	  cols = true, shading = 1, mval = 1e-6);
  store (fp);
  plain();
}

event adapt (i++)
  adapt_flow (ue, 99, 1);

event logger (i++) {
  fprintf (stderr, "%d %g %d %d %d %d %ld %d\n", i, t, mgp.i, 
           mgp.nrelax, mgp2.i, mgp2.nrelax, grid->tn, grid->maxdepth);
}

event stop (t = 35);

