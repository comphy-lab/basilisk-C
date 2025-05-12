/**
# The buoyant plume

The aim of this code was to get to use the Navier-Stokes solver.
Lots of it derives from [Antoon van Hooft](../Antoonvh/README)'s work.

The buoyancy is due to the gradient of concentration of a chemical, equations are the same than with the traditional Boussinesq approach.

The Péclet number in the heat equation is now $\Pi$, an equivalent for solute flows.
*/



#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"


/**  */
double Pi= 2.e4;

double be = 2.5e-3;
int maxlevel = 11;

scalar A[];
scalar * tracers = {A};
face vector diff_tra[];


face vector av[];

/** Arguable BC */
A[bottom]   = dirichlet(fabs(x)<0.01);
A[top]      = dirichlet(0.);
u.t[bottom] = dirichlet(0.);


/**  */
int main(){
  origin(-0.5,0.);
  a = av;

  display_control (maxlevel,6,12);
  display_control (Pi,10.,100000.);
  run();
}

event init(i=0){
  DT = 0.01;
  refine((y<0.04 && fabs(x) < 0.02 ) && level <= maxlevel);
  foreach(){
    A[]  = 0;//(fabs(x)<0.1)*(y<0.1);
    u.y[]= 0;//(fabs(x)<0.1)*y*(1-y);
  }
  boundary({A});
  
}
/**  */
event acceleration(i++){
  foreach_face(y)
    av.y[] = (A[] + A[0,-1])/2.;
}

//event whereWeAre(t+=0.1)
//  fprintf(stderr,"%d %g \n",i,t);


event properties (i++){
  foreach_face()
    diff_tra.x[] = 1./Pi;
}

event tracer_diffusion (i++){
  diffusion(A,dt,diff_tra);
}


/**  */
event adapt(i++){
  adapt_wavelet ((scalar *){A,u}, (double[]){be,be,be},maxlevel,3);
}

/**  */
void black_body (double cmap[NCMAP][3])
{
  /* black body color map from:
   * http://www.kennethmoreland.com/color-advice/
   */
  static double basemap[33][3] = {
    {0.0,0.0,0.0},
    {0.0857913205762,0.0309874526184,0.0173328711915},
    {0.133174636606,0.0588688899571,0.0346802666087},
    {0.180001956037,0.0730689545154,0.0515393237212},
    {0.22981556179,0.0840603593119,0.0647813713857},
    {0.281397607223,0.093912584278,0.075408501413},
    {0.334521638801,0.102639499627,0.0842454688083},
    {0.388957802186,0.110254429637,0.0927990674821},
    {0.444611925648,0.116732501721,0.101402659637},
    {0.501422312285,0.122025816585,0.110058408122},
    {0.559331322331,0.126067584009,0.118767796491},
    {0.618285970576,0.128767919785,0.127531801155},
    {0.678237857955,0.130007052818,0.136351016263},
    {0.712849583079,0.181721849923,0.13081678256},
    {0.743632057947,0.232649759358,0.120991817028},
    {0.774324938583,0.279315911516,0.108089917959},
    {0.804936242903,0.323627020047,0.0907961686083},
    {0.835473266757,0.366524681419,0.0662363460741},
    {0.865942668698,0.408541395043,0.026029485466},
    {0.876634426153,0.46401951695,0.0173065426095},
    {0.883455346031,0.518983528803,0.0149628730405},
    {0.88905246237,0.572164381169,0.013499801006},
    {0.893375939063,0.624108797455,0.0130334871745},
    {0.89637036663,0.675180034619,0.013680092215},
    {0.897973818846,0.725630730259,0.015555776796},
    {0.898116710502,0.775642817733,0.0187767015864},
    {0.896720396485,0.825350944866,0.023459027255},
    {0.927670131094,0.859991226192,0.319086199143},
    {0.956158602738,0.893933112845,0.503316730316},
    {0.97827065392,0.92856476667,0.671307024002},
    {0.993196411712,0.963913323002,0.83560909192},
    {1.0,1.0,1.0},
  };
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(31 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

#if MOVIE
scalar b[];

event movie(t+=0.01){
  b = A;
  boundary ({b});
  scalar db[];
  foreach(){
    db[] = 0.;
    foreach_dimension()
      db[] += sq((b[1] - b[-1])/(2*Delta));
    if (db[] > 0.)
      db[] = log (sqrt (db[]) + 1.);
  }
  boundary({db});
  output_ppm(db, file = "plume.mp4", n = (1 << (maxlevel)), map = black_body,
	     linear = true, box = {{-0.2, 0.}, {.2, 1}}, min = 0., max = 6.);
  output_ppm(A,file = "A.mp4",box = {{-0.2,0.},{0.2,1.}},min=0,max = 1, linear = true);
  output_ppm(u.y, file = "uy.mp4",min=-1.,max = 1,linear = true);
}
#endif // MOVIE

scalar b[];

event picture(t=7){
  b = A;
  boundary ({b});
  scalar db[];
  foreach(){
    db[] = 0.;
    foreach_dimension()
      db[] += sq((b[1] - b[-1])/(2*Delta));
    if (db[] > 0.)
      db[] = log (sqrt (db[]) + 1.);
  }
  boundary({db});
  output_ppm(db, file = "plume.png", n = (1 << (maxlevel)), map = black_body,
	     linear = true, box = {{-0.2, 0.}, {.2, 1}}, min = 0., max = 6.);
}
/** 
## Output

![Picture](plume_buoyant/plume.png){width="50%"}

## References:

to be done 
*/