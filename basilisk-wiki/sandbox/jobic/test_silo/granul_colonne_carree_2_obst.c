
#include "grid/quadtree.h"
#include "granular.h"
#include "utils.h"
#include "view.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

#define LEVEL 8

int    numFOut=0;
double H0=0.,visco=1e-3,normeG=0.,Lx=1.,milieu=0.,W=0.;
scalar champs_solide[];
coord center = {0.5,0.35,0.};

/**
We define here an obstacle, borrowed from aberny's sandbox (cgsBool.c)*/
double cubeF (double x, double y, double z, coord center, double size)
{

  /**
  Our cube is defined as the intersection of 4 orthogonal planes. We
  define the first two planes, $P1_{Plus}$ and $P1_{Minus}$.
  
  We then define P1 has: 
  $$ P1 = P1_{Plus}\cap -P1_{Minus}$$*/

  double P1_Plus = x - size/2. - center.x;
  double P1_Minus = x + size/2. - center.x;
  double P1 = max (P1_Plus, -P1_Minus);

  /**
  We apply the same process to otain P2. */

  double P2_Plus = y - size/2. - center.y;
  double P2_Minus = y + size/2. - center.y;
  double P2 = max (P2_Plus, -P2_Minus);

  /**
  The cube is finally given by:

  $$P1 \cap P2 $$*/

  double c=  max (P1, P2);
  return c;
}

double sphere (double x, double y, double z, coord center, double radius) 
{
  return (sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
          - sq (radius));
}


double geometry (double x, double y, double z)  {
  return cubeF (x, y, z, center, 0.11); 
}

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

  u.n[bottom] = (fabs(x - milieu) < W/2.) ? neumann(0) : dirichlet(0);
  u.t[bottom] = (fabs(x - milieu) < W/2.) ? neumann(0) : dirichlet(0);
  p[bottom] = (fabs(x - milieu) < W/2.) ? dirichlet(0) : neumann(0);
  pf[bottom] = (fabs(x - milieu) < W/2.) ? dirichlet(0) : neumann(0);
  
  DT=0.0001;
  run(); 
}

event init (t = 0) {
  scalar phi[];
  
  foreach_vertex()
    phi[] = x < Lx ? H0 - y : 0.;
  fractions (phi, f);

  foreach() {
    if (geometry(x,y,z)<0)
      champs_solide[]=0.;
    else
      champs_solide[]=1.;
  }

  foreach() {
    u.x[] = 0.;
    p[] = y<0.05 || champs_solide[]<1  ? 0 : max(H0 - y,0.);
    //p[] = 0.;
  }

  /**/
  char name[80];
  sprintf(name, "./init_sol_at_%g",t);
  output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);
  /**/
  
}

#if 0
event mydebug (i++) {
  char name[80];
  sprintf(name, "./sol_at_%d_iter",i);
  output_vtu((scalar *) {f,p}, (vector *) {u}, name);
}
#endif


/**
convergence outputs
*/
void mg_print(mgstats mg)
{
	if (mg.i > 0 && mg.resa > 0.)
		fprintf(stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
		exp(log(mg.resb / mg.resa) / mg.i));
}

event zeroesInSolids(i++) {
  foreach() {
    foreach_dimension(){
      u.x[]*=champs_solide[];
    }
    //p[]*=champs_solide[];
  }
}

event velocity (t+=0.1) {

  double V = 1;
  V = 0;
  foreach(reduction(+:V))
  	  V += f[] * Delta * Delta * champs_solide[];
  if (t >= 0.) fprintf(stdout, "%lf %lf \n", t, V);
  fflush(stdout);

}

event out_vtu (t+=0.1) {

  scalar press[];
  foreach()
    press[] = champs_solide[] * p[];

  char name[80];
  if (numFOut<10)
    sprintf(name, "./sol_at_00%1d",numFOut);
  else if (numFOut<100)
    sprintf(name, "./sol_at_0%2d",numFOut);
  numFOut++;
  output_vtu((scalar *) {f,p,press}, (vector *) {u}, name);

}
/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++)
  adapt_wavelet ((scalar*){f,u}, (double[]){.75,3e-3,3e-3}, LEVEL, 5);

event stop (t = 3){

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

