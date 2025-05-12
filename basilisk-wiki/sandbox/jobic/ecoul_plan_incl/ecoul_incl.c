
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "vof.h"
#include "utils.h"
#include "view.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

/**
passive fluid small density to preserve 0 pressure
and small viscocity
*/
#define RHOF 1e-4
#define mug  1e-5
#define LEVEL 8

scalar f[];
scalar * interfaces = { f };
face vector alphav[];
face vector muv[];
scalar rhov[];
scalar champs_solide[];
coord center = {-0.5,-0.25,0.};
double H0,angle,amplitude=1,visco=1e-3,normeG,rect=0;

/**
We define here an obstacle, borrowed from aberny's sandbox (cgsBool.c)*/
double cubeF (double x, double y, double z, coord center, double size)
{

  /**
  Our cube is defined as the intersection of 4 orthogonal planes. We
  define the first two planes, $P1_{Plus}$ and $P1_{Minus}$.
  
  We then define P1 has: 
  $$ P1 = P1_{Plus}\cap -P1_{Minus}$$*/

  double P1_Plus = x - size/2. + center.x;
  double P1_Minus = x + size/2. + center.x;
  double P1 = max (P1_Plus, -P1_Minus);

  /**
  We apply the same process to otain P2. */

  double P2_Plus = y - size/2. + center.y;
  double P2_Minus = y + size/2. + center.y;
  double P2 = max (P2_Plus, -P2_Minus);

  /**
  The cube is finally given by:

  $$P1 \cap P2 $$*/

  double c=  max (P1, P2);
  return c;
}


int main(int argc, char * argv[]) {
  if (argc > 1)
    rect=atof(argv[1]);
  L0 = 1.;
  H0=0.5;
  angle=30.*M_PI/180;
  amplitude=0.001;
  init_grid(128);
  
  const face vector mdpdx[] = {amplitude*sin(angle),-amplitude*cos(angle)};
  normeG=sqrt(sq(amplitude*sin(angle))+sq(amplitude*cos(angle)));
  a = mdpdx;
  mu = muv;
  rho = rhov;
  alpha = alphav;

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
  //p[bottom] = neumann(cos(angle));
  //pf[bottom] = neumann(cos(angle));


  DT=0.01;
  run(); 
}

double geometry (double x, double y, double z)  {
  return cubeF (x, y, z, center, 0.1); 
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

  
  /**/
  char name[80];
  sprintf(name, "./init_sol_at_%g",t);
  output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);
  /**/
  
}

/**
total density
*/

#define drho(f) ((f) + RHOF*(1. - (f)))
event properties(i++) {
  trash({ alphav });
  scalar eta[];

  scalar fa[];
  foreach()
    fa[] = (4.*f[] +
    2.*(f[-1, 0] + f[1, 0] + f[0, -1] + f[0, 1]) +
    f[1, 1] + f[-1, 1] + f[1, -1] + f[-1, -1]) / 16.;
  boundary({ fa });

  foreach() {
    eta[] = mug;
    if (f[] > 0.) {
      eta[] = visco;
    }
  }
  boundary({ eta });
  
  foreach_face() {
    double fm = (fa[] + fa[-1,0]) / 2.;
    muv.x[] = fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug;
    alphav.x[] = 1./drho(fm);
  }
  foreach()
    rhov[] = drho(fa[]);
  boundary({ muv, alphav, rhov });
}



scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}

event conv (t += 0.1; i <= 100000) {
  update_perf();
  fprintf(stderr,"%d %g %ld\n",i,t,perf.tnc);
  fflush(stderr);

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

    sprintf(name, "./fin_sol_at_%g",t);
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
  sprintf (name, "sol_at_%g.png",t);
  output_ppm(f, file = name, min = 0, max = 4, spread = 2, n = 512, linear = true,
  	  box = { { 0, 0 }, { 1, 1 } });

  sprintf(name, "./sol_at_%g",t);
  output_vtu((scalar *) {f,p,champs_solide}, (vector *) {u}, name);

}


/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++)
  adapt_wavelet ((scalar*){f,u}, (double[]){0.01,3e-3,3e-3}, LEVEL, 6);

event stop (t = 20){

  /**/
  char name[80];
  sprintf(name, "./sol_at_%g",t);
  output_vtu((scalar *) {f,p}, (vector *) {u}, name);
  /**/

  return 1; 
}
