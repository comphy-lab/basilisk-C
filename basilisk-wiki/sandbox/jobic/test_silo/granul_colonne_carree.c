
#include "grid/quadtree.h"
#include "granular.h"
#include "utils.h"
#include "view.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

#define LEVEL 8

int    numFOut=0;
double H0,visco=1e-3,normeG,Lx=1.,milieu,W=0;

int main(int argc, char * argv[]) {

  L0 = 1.;
  H0=0.75;
  milieu=Lx/2.;
  W=Lx/3.;
  init_grid(128);
  
  const face vector mdpdx[] = {0.,-1.};
  a = mdpdx;
  mu = muv;
  rho = rhov;
  alpha = alphav;

  u.n[right] = dirichlet(0);
  u.t[right] = dirichlet(0);
  f[right] = neumann(0);
  
  u.n[left] = dirichlet(0);
  u.t[left] = dirichlet(0);

  u.n[top] = neumann(0);
  /*
  As suggested in granular.h*/
  p[top] = dirichlet(-RHOF*Lx);
  
  u.n[bottom] = neumann(0);
  u.t[bottom] = neumann(0);
  p[bottom] = dirichlet(0);
  
  DT=0.001;
  run(); 
}

event init (t = 0) {
  scalar phi[];
  
  foreach_vertex()
    phi[] = x < Lx ? H0 - y : 0.;
  fractions (phi, f);
  
  foreach() {
    u.x[] = 0.;
    p[] = y<0.05 ? 0 : max(H0 - y,0.);
  }

  /**/
  char name[80];
  sprintf(name, "./init_sol_at_%g",t);
  output_vtu((scalar *) {f,p,phi}, (vector *) {u}, name);
  /**/
  
}

/**
convergence outputs
*/
void mg_print(mgstats mg)
{
	if (mg.i > 0 && mg.resa > 0.)
		fprintf(stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
		exp(log(mg.resb / mg.resa) / mg.i));
}

event velocity (t+=0.1) {

  double V = 1;
  V = 0;
  foreach(reduction(+:V))
  	  V += f[] * Delta * Delta;
  if (t >= 0.) fprintf(stdout, "%lf %lf \n", t, V);
  fflush(stdout);

  char name[80];
  if (numFOut<10)
    sprintf(name, "./sol_at_00%1d",numFOut);
  else if (numFOut<100)
    sprintf(name, "./sol_at_0%2d",numFOut);
  numFOut++;
  output_vtu((scalar *) {f,p}, (vector *) {u}, name);

}


/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){3e-3,3e-3}, LEVEL, 5);

event stop (t = 1.2){

  /**/
  char name[80];
  sprintf(name, "./sol_FIN_%g",t);
  output_vtu((scalar *) {f,p}, (vector *) {u}, name);
  /**/

  return 1; 
}

#if 1
event movie (t += 0.05) {
    static FILE * fp1 = popen ("ppm2mpeg > level.mp4.ppm", "w");
    scalar l[];
    foreach()
    l[] = level;
    boundary ({l});
    output_ppm (l, fp1, min = 0, max = LEVEL,
                n = 512, box = {{0,0},{Lx,Lx}});
    
    foreach()
    l[] = f[]*(1 + sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    
    static FILE * fp2 = popen ("ppm2mpeg > velo.mp4.ppm", "w");
    output_ppm (l, fp2, min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{Lx,Lx}});
}
#endif

