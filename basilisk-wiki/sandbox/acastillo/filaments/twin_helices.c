#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"
#include "view.h"
#include "lambda2.h"

#define MINLEVEL 3
#define RAD (sqrt(sq(x) + sq(y)))

int n_seg;
double as;
double Hs1, Rs1, n_turns1;
double Hs2, Rs2, n_turns2;

u.n[top] = neumann (0); u.t[top] = dirichlet (0);
u.n[bottom] = neumann (0); u.t[bottom] = dirichlet (0);
u.n[left] = neumann (0); u.t[left] = dirichlet (0);
u.n[right] = neumann (0); u.t[right] = dirichlet (0);



int main() {
  L0 = 8; ;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.001;
  N = 1<<MINLEVEL;

  double reynolds= 16100;
  const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds};
  mu = muc;

  periodic(back);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);
}

#ifdef _INITIAL
#include "twin_helices_init.h"
event stop (i = 2)
  if (pid()==0) printf( "... Done generating initial conditions! \n" );
#else

event init (t = 0){
  if (pid()==0) printf( "Restoring from previous run \n" );
#ifdef _BETA
  char name[80];
  sprintf (name, "./case%d/level%d/0_init/dump", _BETA, MAXLEVEL);
  //sprintf (name, "./case%d/level%d/3_relax/dump", _BETA, MAXLEVEL);
  printf(name);
  restore(name);

  scalar l2[];
  lambda2 (u, l2);

  view (camera="iso");
  box();
  cells();
  squares ("u.x", linear = false, n = {1,0,0} );
  isosurface ("l2", 0);
  save ("twin_helices1.png");

#else
  restore("dump");
#endif

  if (pid()==0) printf( "... Done restoring! \n" );

}

/**
# Outputs
*/
double tsample1=0.005, tsample2=0.025, tsample3=0.125, tsample4=1.0;

#include "3d/ellipticity.h"
#include "twin_helices_export.h"
#include "twin_helices_moments.h"

event stop (t = 0.5)
  dump (list = (scalar *){u});
#endif

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif
