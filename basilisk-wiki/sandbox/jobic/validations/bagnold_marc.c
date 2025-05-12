
#include "grid/quadtree.h"
#include "granular_without_vof.h"
//#include "navier-stokes/double-projection.h"
//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

/**
Not parallel yet (with MPI) 
*/
double error_rms_vit(void);

#define mug  1e-5
#define LEVEL 7

double Lx=14.,angle=0.43,oldErreur=HUGE;
long int itmax=1000000;

double computeSolEx(double posy) {
    double H=Lx;
    double U0=I0*sqrt(cos(angle))*(tan(angle)-mus)/(mus+dmu - tan(angle));
    return ((pow(H,1.5) - pow(H - posy, 1.5))*2./3.*U0) ;
}

int main(int argc, char * argv[]) {
  L0 = Lx; 
  size(Lx);
  origin(0,0,0);
  lambdaR=1e-4;
  N=1<<5;
  
  const face vector mdpdx[] = {sin(angle),-cos(angle)};
  a = mdpdx;
  
  u.t[top] = neumann(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  p[top] = dirichlet(0);

  periodic (right);
  
  DT=0.1;
  run(); 
}

event init (t = 0) {
  
  foreach() {
    u.x[] = computeSolEx(y);
    u.y[] = 0.;
    p[] = cos(angle)*(1-y);
  }

  /**/
  char name[80];
  sprintf(name, "./init_sol_at_%g",t);
  output_vtu((scalar *) {p}, (vector *) {u}, name);
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

/**
 old value of the velocity is saved
*/
scalar un[];
event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event conv (t+=1) {
  double du = change (u.x, un);
  double err=error_rms_vit(),errDu=fabs(err-oldErreur)/oldErreur;
  oldErreur=err;

  fprintf(stdout,"%g %g %g %g\n",t,du,err,errDu);
  fflush(stdout);
  if (i > 0 && du < 1.0e-5) {
    scalar solEx[],myI[],dvdx[];
    vector uu[];
    
    foreach() {
      double D2 = 0.;
      foreach_dimension() {
        double dxx = u.x[1,0] - u.x[-1,0];
        double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
        D2 += sq(dxx) + sq(dxy);
      }
      D2 = sqrt(2*D2)/(2.*Delta); // this D2 is sqrt(2) D2
      myI[] = D2/sqrt(p[]); //Adim, for dim, add dg*
      dvdx[]=D2; // this D2 is sqrt(2) D2
    
      solEx[]=computeSolEx(y);
    }
    
    char name[80];
    sprintf(name, "./sol_at_%g",t);
    output_vtu((scalar *) {p,solEx,myI,dvdx}, (vector *) {u}, name);

    return 1; /* stop */
  }
}


/**
We adapt according to the error on the velocity field. 
*/
/*
event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){3e-3,3e-3}, LEVEL, 5);
*/

event stop (i = itmax){

  /**/
  char name[80];
  sprintf(name, "./sol_FIN_%g",t);
  output_vtu((scalar *) {p}, (vector *) {u}, name);
  /**/

  return 1; 
}

/**
Not parallel yet (with MPI) 
*/
double error_rms_vit(void) {
  double rms=0;
  int cpt=0;
  
  foreach(reduction(+:rms) reduction(+:cpt)) {
    double val,valEx;
    val=u.x[];
    valEx=computeSolEx(y);
    rms+=(val-valEx)*(val-valEx);
    //printf("%g %d %g %g\n",t,cpt,val,valEx);
    cpt++;
  }
  return sqrt(rms)/cpt;
}


#if 0
event movie (t += 0.05) {
    static FILE * fp1 = popen ("ppm2mpeg > level.mp4.ppm", "w");
    scalar l[];
    foreach()
    l[] = level;
    boundary ({l});
    output_ppm (l, fp1, min = 0, max = LEVEL,
                n = 512, box = {{0,0},{Lx,Lx}});
    
    foreach()
    l[] = sqrt(sq(u.x[]) + sq(u.y[]));
    boundary ({l});
    
    static FILE * fp2 = popen ("ppm2mpeg > velo.mp4.ppm", "w");
    output_ppm (l, fp2, min = 0, max = 2., linear = true,
                n = 512, box = {{0,0},{Lx,Lx}});
}
#endif

