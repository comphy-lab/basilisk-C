/**
# Elliptic instability

A dipolar vortex whose streamlines are elliptical can be subject to the so-called 3-dimensional "elliptic instability". It is demonstrated here with the Lamb-Chaplygin dipole model.

![Volumetric rendering of the $\lambda_2$ field. Side views are rendered from dumps in post processing](https://www.antoonvanhooft.nl/media/elliptical.mp4)(width=800px)

*/
#define RKORDER (3)
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;
#include "lambda2.h"
#include "bwatch.h"
#include "view.h"
double ue = 1e-4, vis = 1.5e-3;

int main(int argc, char ** argv) {
  if (argc > 1)
    ue = atof (argv[1]);
  if (argc > 2)
    vis = atof (argv[2]);
  foreach_dimension()
    periodic(left);
  L0 = 4;
  X0 = Y0 = Z0 = -L0/2;
  N = 1 << 5;
  const scalar visc[] = vis;
  nu = visc;
  run();
}

double R = 1, k = 3.8317;
double U = 1;

#define RAD (sqrt(sq(y) + sq(z)))
#define SINT (y/RAD) 

event init (t = 0) {
  vector omega[], psi[];
  foreach_dimension()
    psi.x.prolongation = refine_4th;
  foreach() {
    foreach_dimension()
      psi.x[] = 0;
  }
  boundary ((scalar*){psi});
  
   vector uc[];
    foreach_dimension()
      uc.x.prolongation = refine_4th;
   // The Lamb-Chaplygin dipole model is not even 4th order smooth...
  do {
    foreach() {
      foreach_dimension()
	omega.x[] = (RAD < R)*noise()/100.;
      omega.x[] =  ((RAD < R)*((-2*U*j1(k*RAD)*SINT/(k*j0(k)))))*sq(k);
    }
    foreach_dimension() {
      stats so = statsf (omega.x);
      foreach()
	omega.x[] -= so.sum/so.volume;
    }
    foreach_dimension() // Omega.x is the only non-trivial problem
      poisson (psi.x, omega.x);
    foreach(){
      foreach_dimension()
	uc.x[] = -((8*(psi.z[0,1] - psi.z[0,-1]) + psi.z[0,-2] - psi.z[0,2]) -
		   (8*(psi.y[0,0,1] - psi.y[0,0,-1]) + psi.y[0,0,-2] - psi.y[0,0,2]))/(12*Delta);
    }
    foreach()
      uc.z[] += 0.6*U;
     vector_to_face (uc); 
     multigrid_restriction((scalar*){u}); 
     printf ("#cells: %ld, depth: %d\n", grid->tn, grid->maxdepth); 
  } while (adapt_flow(ue, 99, 1).nf > (grid->tn/100));
  printf ("#cells: %ld, depth: %d\n", grid->tn, grid->maxdepth);
  project (u, p);
}

#define FUNC(x) (exp(-x) + x - 1)

event dumper (t = {7, 10, 15}) {
  vector uc[];
  face_to_vector (uc);
  scalar l2[];
  lambda2 (uc, l2);
  foreach() 
    l2[] = l2[] < 0 ? FUNC(-l2[]): 0;
  char fn[99];
  sprintf (fn, "dump_%d", (int)t);
  dump(fn);
}

event mov (t = 7; t += 0.02) {
  static FILE * fp = popen ("ppm2mp4 elliptic.mp4", "w");
  vector uc[];
  face_to_vector (uc);
  scalar l2[];
  lambda2 (uc, l2);
  foreach() 
    l2[] = l2[] < 0 ? FUNC(-l2[]): 0;
  boundary ({l2});
  watch (fov = L0/0.8, O = {0.0001, 0, -20},
	 poi = {0.001, 0, 0},
  	 nx = 1024, ny = (30*24));
  volume (l2, sc = 4, min = -30, max = 30,
	  cols = true, shading = 1, mval = 1e-3);
  store (fp);
  plain();
#if 0  
  static FILE * fpc = popen ("ppm2mp4 elliptic_colors.mp4", "w");
  
  vector omg[];
  vorticityf3(u, omg);
  // Direction only?
  foreach() {
    double OMG = 0;
    foreach_dimension() {
      OMG += sq(omg.x[]);
      omg.x[] = fabs(omg.x[]);
    }
    foreach_dimension()
      ;//omg.x[] = OMG ? fabs(omg.x[])/sqrt(OMG) : 0.5;
  }
  volume (l2, sc = 4, max = sq(k)/1.5,
	  cols = true, shading = 1, colorv = omg, mval = 1e-3);
  store (fpc);
  plain();
#endif
  
  isosurface ("l2", 2);
  cells();
  save ("iso_surf.mp4");
  save ("snapshot.png");
}

event adapt (i++)
  adapt_flow (ue, 99, 1);

event logger (i++) {
  fprintf (stderr, "%d %g %d %d %d %d %ld %d\n", i, t, mgp.i, 
           mgp.nrelax, mgp2.i, mgp2.nrelax, grid->tn, grid->maxdepth);
}

event stop (t = 15);

