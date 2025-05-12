
#include "grid/quadtree.h"
#include "granular.h"
#include "navier-stokes/double-projection.h"
#include "utils.h"
#include "view.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

#define LEVEL 7

scalar champs_solide[];
coord center = {0.5,0.25,0.};
int    numFOut=0;
double H0,angle,amplitude=0.001,visco=1e-3,normeG,rect=0;

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

double geometry (double x, double y, double z)  {
  return cubeF (x, y, z, center, 0.11); 
}


int main(int argc, char * argv[]) {
  if (argc > 1)
    rect=atof(argv[1]);
  if (argc > 2)
    amplitude=atof(argv[1]);
  L0 = 1.;
  H0=0.5;
  angle=30.*M_PI/180;
  init_grid(128);
  
  const face vector mdpdx[] = {amplitude*sin(angle),-amplitude*cos(angle)};
  normeG=sqrt(sq(amplitude*sin(angle))+sq(amplitude*cos(angle)));
  a = mdpdx;

  /**
   Boundary conditions are periodic
  */
    periodic (right);

  /**
   slip at the bottom
  */
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
  p[top] = dirichlet(0);
  pf[top] = dirichlet(0);

  /**
   no slip at the bottom
  */
  u.n[bottom] = dirichlet(0);
  uf.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);

  DT=0.001;
  run(); 
}

event init (t = 0) {
  scalar phi[];
  
  foreach_vertex()
    phi[] = H0 - y;
  fractions (phi, f);
  
  if (rect) {
    foreach() {
      if (geometry(x,y,z)<0)
    	champs_solide[]=0.;
      else
    	champs_solide[]=1.;
    }
  } else {
    foreach() 
      champs_solide[]=1.;
  }
  
  /**
   to help convergence, we init a linear velocity porfile
  */
  double Vmax=normeG*H0*H0/8/(mug*100);
  foreach() {
    u.x[] = y>0.5 || champs_solide[] < 1 ? 0 : 2*Vmax*y;
    p[] = champs_solide[] < 1 ? 0 : max(H0 - y,0.);
  }

  char name[80];
  sprintf(name, "./init_sol_at_%g",t);
  output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);
}

scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}

event conv (t += 0.1; i <= 100000) {
  double du = change (u.x, un);
  fprintf(stdout,"t= %g %g %g %g\n",t,interpolate (u.x, 0.5, 0.25),interpolate (u.y, 0.5, 0.25),du);
  fflush(stdout);
  if (i > 1 && du < 1.0e-4) {
    scalar nu[];
    foreach()
      nu[] = norm(u);
     
    char name[80];
    view (bg = {1,1,1},
  	  width = 600, height = 600);
    squares ("nu", linear = true);
    //cells();
    sprintf (name, "sol.png");
    save (name);    

    sprintf(name, "./FIN_sol_at_%g",t);
    output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);

    return 1; /* stop */
  }
}

event zeroesInSolids(i++) {
  foreach() {
    foreach_dimension(){
      u.x[]*=champs_solide[];
    }
    //p[]*=champs_solide[];
  }
}


event defaults (t+=1) {
  char name[80];
  scalar press[];
  
  foreach()
    press[] = champs_solide[] * p[] * f[];


  if (numFOut<10)
    sprintf(name, "./sol_at_00%1d",numFOut);
  else if (numFOut<100)
    sprintf(name, "./sol_at_0%2d",numFOut);
  else if (numFOut<1000)
    sprintf(name, "./sol_at_%2d",numFOut);
  numFOut++;
  output_vtu((scalar *) {f,p,press,champs_solide}, (vector *) {u}, name);
}


/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++)
  adapt_wavelet ((scalar*){f,u}, (double[]){0.75,3e-3,3e-3}, LEVEL, 7);

event stop (t = 30){

  /**/
  char name[80];
  sprintf(name, "./FIN_sol_at_%g",t);
  output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);
  /**/

  return 1; 
}
