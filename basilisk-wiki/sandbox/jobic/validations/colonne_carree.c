
#include "grid/quadtree.h"
#include "granular_without_vof.h"
//#include "navier-stokes/double-projection.h"
#include "utils.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"

#define LEVEL 8

double Lx=28.,milieu,fixedPressure=28,oldErreur=HUGE;

double error_rms_vit(void);

double computeSolEx(double pval, double posx, double pmilieu) {
  double val, xx=fabs(posx), x0=mus*pval, dmuP=dmu*pval;
  if (xx <  x0) {
    val=-I0*sqrt(pval)*((x0-pmilieu)
                       +dmuP*log(dmuP/(dmuP-(pmilieu-x0))));
  } else {
    val=-I0*sqrt(pval)*((xx-pmilieu)
                       +dmuP*log((dmuP-(xx-x0))/(dmuP-(pmilieu-x0))));
  }
  return val;
}

double computeShearEx(double pval, double posx) {
  double val, xx=fabs(posx), x0=mus*pval, dmuP=dmu*pval;
  if (xx <  x0) 
    val = I0*sqrt(pval)*(x0-x0)/(dmuP-(Lx/2-x0));
  else 
    val = I0*sqrt(pval)*(xx-x0)/(dmuP-(Lx/2-x0));
  return val;
}

int main(int argc, char * argv[]) {
  L0 = Lx; 
  size(Lx);
  milieu=Lx/2.;
  origin(-milieu,0,0);
  lambdaR=1e-4;
  N=1<<LEVEL;
  
  const face vector mdpdx[] = {0.,-1.};
  a = mdpdx;
  
  u.n[right] = dirichlet(0);
  u.t[right] = dirichlet(0);
  
  u.n[left] = dirichlet(0);
  u.t[left] = dirichlet(0);

  periodic (top);
  //stokes=true;
  
  DT=0.001;
  run(); 
}

event init (t = 0) {
  
  foreach() {
    u.x[] = 0.;
    //u.y[] = computeSolEx(fixedPressure,x,14.)*(rand()/RAND_MAX/0.2 + 0.9);
    //u.y[] = -2;
    u.y[] = computeSolEx(fixedPressure,x,14.);
    p[] = fixedPressure;
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
/*
void mg_print(mgstats mg)
{
	if (mg.i > 0 && mg.resa > 0.)
		fprintf(stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
		exp(log(mg.resb / mg.resa) / mg.i));
}
*/
/**
 old value of the velocity is saved
*/
scalar un[];
event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event conv (i++) {
  double du = change (u.x, un);
  double err=error_rms_vit(),errDu=fabs(err-oldErreur)/oldErreur;
  oldErreur=err;
  stats s = statsf (u.y);

  fprintf(stdout,"%g %g %g %g %g\n",t,du,err,errDu,s.min);
  fflush(stdout);
  if (t > 10 && du < 1.0e-5) {
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
    
      solEx[]=computeSolEx(fixedPressure,x,14.);
    }
    
    char name[80];
    sprintf(name, "./sol_at_%g",t);
    output_vtu((scalar *) {p,solEx,myI,dvdx}, (vector *) {u}, name);

    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double x = -Lx/2; x < Lx/2; x += 1./pow(2.,LEVEL)) {
      double SolEx=computeSolEx(fixedPressure,x,14.);
      fprintf (fp, "%g %g %g %g %g %g %g %g \n", x, interpolate (u.y, x, L0/2), interpolate (shear, x, L0/2),
               SolEx, interpolate (p, x, L0/2),
               interpolate (eta, x, L0/2), interpolate (foo, x, L0/2),computeShearEx(fixedPressure,x) );
    }
    fclose (fp);

    return 1; /* stop */
  }
}

/**
Not parallel yet (with MPI) 
*/
double error_rms_vit(void) {
  double rms=0;
  int cpt=0;
  
  foreach(reduction(+:rms) reduction(+:cpt)) {
    double val,valEx;
    val=u.y[];
    valEx=computeSolEx(fixedPressure,x,14.);
    rms+=(val-valEx)*(val-valEx);
    //printf("%g %d %g %g\n",t,cpt,val,valEx);
    cpt++;
  }
  return sqrt(rms)/cpt;
}


event stop (t = 20){

  /**/
  char name[80];
  sprintf(name, "./sol_FIN_%g",t);
  output_vtu((scalar *) {p}, (vector *) {u}, name);
  /**/

  return 1; 
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

