/**
#Solidification of an undercooled liquid

We want to study the dendritic growth of a solid in an undercooled liquid. We
simulate the diffusion of two tracers separated by an embedded
boundary. The interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
$$

where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The boundaries of the domain are set at $T_L = -1$. The ice particle is initially
at $T_S = 0$.

The temperature on the interface is derived from the Gibbs-Thomson relation:
$$
 T_{\Gamma} = T_{m} - \epsilon_{\kappa} * \kappa - \epsilon_{v} v_{pc}
$$

~~~gnuplot Interface
set key top right
set xrange [-2:2]
set size ratio -1
set xlabel 'x'
set ylabel 'y'
set object 1 circle front at 0.,0. size 1.1 fillcolor rgb "black" lw 1 dt 2
plot 'out' w l lw 0.8 t 'Interface' 
~~~



*/
#define ANISO 1
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define Pi 3.14159265358979323846


#include "../../ghigo/src/myembed.h"
#include "../double_embed-tree.h"

#include "../advection_A.h"
#include "../../ghigo/src/mydiffusion.h"

#include "fractions.h"
#include "curvature.h"


double fac1(coord n, double eps4){
  double theta = atan2(n.y, n.x);
  return 1.+eps4*(8./3.*sq(sq(sin(3*(theta-Pi/2.)))) - 1. );
}
#include "view.h"
#include "../level_set.h"
#include "../LS_curvature.h"
#include "../LS_advection.h"


#define T_eq         0.
#define TL_inf       -0.8
#define TS_inf       -0.8

/**
Setup of the numerical parameters
*/
int MAXLEVEL = 9;
int MINLEVEL = 4;
double H0;

/**
Setup of the physical parameters + level_set variables
*/
scalar TL[], TS[], dist[];
vector vpc[],vpcf[];

scalar * tracers    = {TL};
scalar * tracers2 = {TS};

scalar * level_set  = {dist};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];
double DT2;


double  latent_heat = 1.;
double  lambda[2]; // thermal capacity of each material
#if Gibbs_Thomson // parameters for the Gibbs-Thomson's equation
double  epsK = 0.001, epsV = 0.001;
#else
double  epsK = 0.000, epsV = 0.000;
#endif

#define GT_aniso 1
/**
Initial interface is a circle, we want to study the effect of the anisotropy on
the curvature constant in the Gibbs-Thomson equation, i.e. we want to introduce
dendritic growth purely due to anisotropy in the GT equation.
*/
#if GT_aniso
double eps4 = 0.4;
#else
double eps4 = 0.;
#endif

int nb_cell_NB;
double  NB_width ;    // length of the NB

int itrecons;

scalar curve[];


#define Pi 3.14159265358979323846

double geometry(double x, double y, coord center, double Radius) {

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = ( sqrt(R2)*(1.+0.2*cos(6*theta)) - Radius);

  return s;
}


double TL_init(double val,double V, double t){
  if(val>V*t)
    return TL_inf+exp(-V*(val-V*t));
    // return TL_inf;
  else
    return 0.;
}

#include "../alex_functions2.h"

/**
$k_{loop}$ is a variable used when the interface arrives in a new cell. Because
we are using embedded boundaries, the solver needs to converge a bit more when
the interface is moved to a new cell. Therefore, the Poisson solver is used
$40*k_{loop} +1$ times.
*/
int k_loop = 0;



int main() {
  L0 = 4.;
  CFL = 0.3;
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  double theta = atan2(0.5,0.5);

  origin(-L0/2., -L0/2.);
  init_grid (1 << 8);
  run();
}


event init(t=0){

  lambda[0]  = 1.;
  lambda[1]  = 1.;
  // DT2        = 0.5;  // diffusion time scale
  DT2        = sq(L0/(1<<MAXLEVEL))/lambda[0];
  nb_cell_NB = 1 << 2 ; // number of cell in the 
                        // narrow band 
  itrecons = 50;
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  coord center = {0.,0.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = 0.1;

  vertex scalar distn[];

  for (int k = 0; k < 5; k++){
    foreach() {
      dist[] = geometry(x,y,center,size);
    }
    boundary({dist});
    restriction({dist});
    cell2node(dist,distn);

    fractions (distn, cs, fs);

    boundary({cs,fs});
    restriction({cs,fs});
    foreach() {
      TL[] = TL_inf;
      TS[] = 0.;
    }
    boundary({dist,TL,TS});
    restriction({dist,TL,TS});
    double lambda1 = lambda[0], lambda2 = lambda[1]; 
    LS_speed(
    dist,latent_heat,cs,fs,TS,TL,T_eq,
    vpc,vpcf,lambda1,lambda2,
    epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
    itrecons = 30,tolrecons = 1.e-12, NB_width);
    event("adapt");
    
  }

  foreach() {
    dist[] = geometry(x,y,center,size);
  }
  boundary ({dist});
  restriction({dist});

  cell2node(dist,distn);

  fractions (distn, cs, fs);

  boundary({cs,fs});
  restriction({cs,fs});

  foreach_face(){
    vpc.x[] = 0.;
  }
  boundary((scalar *){vpc});

#if CURVE_LS
  curvature_LS(dist, curve);
#else
  curvature(cs,curve);
  boundary({curve});
#endif

  foreach() {
    TL[] = TL_inf;
    TS[] = 0.;
  }
  boundary({TL,TS});
  restriction({TL,TS});

  myprop(muv,fs,lambda[0]);
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);

}

event velocity(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = itrecons,tolrecons = 1.e-4, NB_width
  );
  double DT3 = timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
  // dump();
  // exit(1);
  // fprintf(stderr, "## %g %g %g\n", t, DT2, DT3);
}

event tracer_diffusion(i++){
  advection_LS(
  dist,
  cs,fs,
  TS,TL,
  vpcf,
  itredist = 5,
  s_clean = 1.e-10,
  NB_width,
  curve
  );
  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});


  boundary({TL});
  myprop(muv,fs,lambda[0]); // MANDATORY, the interface has moved !!!!!
  mgT = diffusion (TL,dt,muv,theta = cm);
  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[1]);
  boundary({TS});
  mgT = diffusion (TS,dt,muv,theta = cm);
  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[0]);
}

/**
This event is the core the of the hybrid level-set/embedded boundary.
*/


event snapshot(t+=3.e-3,last){
  scalar normvpc[];
  foreach(){
    coord n;
    if(cs[]*(1-cs[])!=0.){
      foreach_dimension()
        n.x = vpc.x[];  
      normvpc[] = norm2(n);
    }
    else{
      normvpc[] = nodata;
    }
  }

  boundary({normvpc});
  fprintf(stderr, "##OUTPUT\n" );
  dump();
}

event interface2(t+=1.e-3,last; t<3.6e-2){
  output_facets(cs,stdout);
}


/**
Mesh adaptation. still some problems with coarsening inside the embedded
boundary.
*/
#if 1
event adapt (i++, last) {
  foreach_cell(){
    cs2[] = 1.-cs[];
  }
  foreach_face(){
      fs2.x[]      = 1.-fs.x[];
  }

  boundary({cs,cs2,fs,fs2});
  fractions_cleanup(cs,fs,smin = 1.e-10);
  fractions_cleanup(cs2,fs2,smin = 1.e-10);
  restriction({cs,cs2,fs,fs2});
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});

  adapt_wavelet ({cs,visu,vpcf.x, vpcf.y},
    (double[]){1.e-3,1.e-2,1.e-2,1.e-2},MAXLEVEL, MINLEVEL);

/* Small on-the-fly outputs*/
  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,2.);
  }
  fprintf(stderr, "##time %g dt %g nbcells %d solid volume %g\n", t, dt, n ,
   solid_volume);
}
#endif
/**

![Animation of the approximated temperature field](cube/visu.mp4)(loop)


*/

