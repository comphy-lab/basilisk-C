/**
#Stable solidification of an initially sinusoidal interface

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_2 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations can be found on the mini_cell test
case.


We output only the interface at different times during the simulation.


~~~gnuplot speed
plot 'log' u 1:2 w p t 'max', 'log' u 1:3 w p t 'min'
#f(x)  = A*exp(-x/B)-C
#f2(x) = A2*exp(-x/B2)-C2
#fit f(x)  'log' u ($1):($2) via A,B,C
#fit f2(x)  'log' u ($1):($3) via A2,B2,C2
#ftitle(a,b,c) = sprintf('%.3f*(exp)^{t/%.2f}-%.3f',a,b,c)
#plot 'test' every 10:10 u 1:2 w p pt 64 ps 1.25 t 'max', \
#     'test' every 10:10 u 1:3 w p pt 64 ps 1.25 t 'min', \
#  f(x) w l dt 2 lw 3 t ftitle(A,B,C), \
#  f2(x) w l dt 2 lw 3 t ftitle(A2,B2,C2)
~~~
*/
#define BICUBIC 1
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0
#define Pi 3.14159265358979323846

#include "embed.h"
#include "../double_embed-tree.h"
#include "../advection_A.h"
#include "diffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../level_set.h"
#include "../LS_advection.h"
#include "../LS_curvature.h"


#define T_eq          0.
#define TL_inf        1.
#define TS_inf       -1.

#define tstart 0.

int MINLEVEL, MAXLEVEL; 
double H0;
double latent_heat;

#define DT_MAX  1.

#define T_eq         0.


#define plane(x, y, n) (y  + 0.2*sin(2.*n*Pi*x))


scalar TL[], TS[], dist[];
scalar * tracers = {TL};
scalar * tracers2 = {TS};

vector vpc[];
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];

#if Gibbs_Thomson
double  epsK = 0.0005, epsV = 0.0005;
#else
double  epsK = 0.000, epsV = 0.000;
#endif
scalar curve[];

#define GT_aniso 0
#if GT_aniso
int aniso = 6;
#else
int aniso = 1;
#endif


double lambda[2];

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB


  
mgstats mg1,mg2;

TL[embed] = dirichlet(Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, aniso));
TL[top]   = dirichlet(TL_inf); 

TS[embed]  = dirichlet(Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, aniso));
TS[bottom] = dirichlet(TS_inf); 

int j;
int k_loop = 0;
/**
The domain is 4 units long, centered vertically. */



int main() {
  periodic(right);

  L0 = 1.;
  origin (-0.5*L0, -0.5*L0);

  
/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
  latent_heat  = 1;
  MAXLEVEL  = 7 ;
  MINLEVEL  = 4 ;

  H0 = 0.5*L0; 
  N = 1 << MAXLEVEL;
  init_grid (1 << MAXLEVEL);
  run();
}

event init(t=0){

  TOLERANCE = 1.e-7;
  DT = 0.2*L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach(){
     dist[] = plane(x,y,2);
  }
  boundary ({dist});
  restriction({dist});

  // LS_reinit(dist, RK2 = 1);

  vertex scalar dist_n[];
  cell2node(dist,dist_n);

  fractions (dist_n, cs, fs);
  fractions_cleanup(cs,fs);
  boundary({cs,fs});
  restriction({cs,fs});

  // curvature(cs,curve);
  // boundary({curve});
  curvature_LS(dist, curve);

  foreach() {
    TL[] = TL_inf;
    TS[] = TS_inf;
  }

  foreach_face(){
    vpc.x[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});
}


event tracer_diffusion(i++){
  int kk;
  mgstats mg1;
  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
  stats s;
  for (kk=1;kk<=4;kk++){
    if(i%2==0){
      boundary({TL});
      mg1 = diffusion(TL, dt, D = muv , theta = cs);
      s = statsf(TL);
    }
    else{
      boundary({TS});
      mg1 = diffusion(TS, dt, D = muv, theta = cs);
      s = statsf(TS);
    }
    // fprintf(stderr, "%g %g %g\n", s.min, s.stddev, s.max);
    if(mg1.resa > TOLERANCE) {
/**
If the calculation crashes (it often does if the Poisson solver does not
converge) we save the last state of the variables
*/
      scalar ffsx[], ffsy[];
      foreach(){
        ffsx[] = fs.x[];
        ffsy[] = fs.y[];
      }
      boundary({ffsx,ffsy});
      dump();
      exit(1);
    }
  }
}

event LS_advection(i++,last){
  if(i%2 == 1){
    invertcs(cs,fs);
    double lambda1 = lambda[0], lambda2 = lambda[1];
    
    advection_LS(
    dist,
    latent_heat,
    cs,fs,
    TS,TL,
    vpc,
    lambda1,lambda2,
    epsK,epsV,aniso,
    curve,
    &k_loop,
    k_limit = 0,
    overshoot = 0.3,
    deltat = 0.45*L0 / (1 << grid->maxdepth),
    itredist = 20,
    itrecons = 300,
    tolrecons = 1.e-5,
    s_clean = 1.e-10,
    NB_width);

    invertcs(cs,fs);
  }
}

#if DOUBLE_EMBED
event double_calculation(i++,last){
  invertcs(cs,fs);
}
#endif

event movies ( i++,last;i<400)
{
  if(i%8 == 1) {

    double y_max = -L0;
    boundary((scalar *){vpc});
    restriction((scalar *){vpc});
    stats s2 = statsf (vpc.y);

    vector h[];
    heights (cs, h);
    boundary((scalar *){h});
    foreach(reduction(max:y_max)){
      if(interfacial(point, cs)){
        double yy = y+Delta*height(h.y[]);
        y_max = max(y_max,yy);
      }     
    }
    fprintf (stderr, "###%.9f %.9f %.9f %g\n",
      t, s2.max, s2.min, y_max);

    if(i%20==1) {
      output_facets (cs, stdout);
    }
  }
}

#if 1
event adapt (i++, last) {
  if(i%2 == 1 ){


    foreach_cell(){
      cs2[] = 1.-cs[];
    }
    foreach_face(){
      fs2.x[]      = 1.-fs.x[];
    }

    boundary({cs,cs2,fs,fs2});
    fractions_cleanup(cs,fs,smin = 1.e-14);
    fractions_cleanup(cs2,fs2,smin = 1.e-14);
    restriction({cs,cs2,fs,fs2});
    int n=0;
    stats s2 = statsf(curve);
    fprintf(stderr, "%g %g %g\n",t, s2.min, s2.max);
    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    restriction({visu});

    adapt_wavelet ({cs,visu,curve},
      (double[]){1.e-3,1.e-3,1.e-3},MAXLEVEL, MINLEVEL);
    foreach(reduction(+:n)){
      n++;
    }
    phase_change_velocity_LS_embed (cs, fs ,TL, TS, vpc, latent_heat,
      lambda,epsK, epsV, aniso);
    // curvature(cs,curve);
    // boundary({curve});
    curvature_LS(dist, curve);


    fprintf(stderr, "##nb cells %d\n", n);
  }
}
#endif


/**
~~~gnuplot Evolution of the interface
set size ratio 0.5
f(x) = 0
plot 'out' w l t 'dt = 0.5*Delta', f(x)
#plot 'out' w l t 'dt = 0.5*Delta', f(x), 'valeurs' u 1:2 w p pt 7 ps 1 lc -1
~~~
*/