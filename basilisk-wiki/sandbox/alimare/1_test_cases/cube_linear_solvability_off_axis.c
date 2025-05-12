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

~~~gnuplot Evolution of the interface (zoom)
set key top right
set size ratio -1
plot 'out' w l lw 0.5 t 'Interface' 
~~~

~~~gnuplot Speed of the top dendrite
set size ratio 0.5
set key top right
set yrange [0:0.1]
f(x) = 0.017
set xlabel 'time'
set ylabel 'velocity'
set ytics scale 3
set mytics 4
unset x2tics
unset y2tics
set grid lw 0.9
set pointsize 1.5 
plot 'log' u 1:2 pt 6 ps 0.4 t 'Tip velocity', \
     '' u 1:3 pt 6 ps 0.4 t 'Tip velocity deduced', \
  f(x) t 'Theoretical value' lt rgb 'red' lw 1.5 
~~~



*/
#define ANISO 1
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define LS_CURVE 1

double fac1(coord n, double eps4){
  double theta = atan2(n.y, n.x);
  return 1.-15.*eps4*cos(4.*(theta+M_PI/4.));
}

#include "../../ghigo/src/myembed.h"
#include "../embed_extrapolate_3.h"

#include "../double_embed-tree.h"

#include "../advection_A.h"
#include "../../ghigo/src/mydiffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../level_set.h"
#include "../LS_curvature.h"
#include "../LS_advection.h"


#define T_eq         0.
#define TL_inf       -0.55
#define TS_inf       -0.55

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
double  epsK = 0.5, epsV = 0.000;
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
double eps4 = 0.05;
#else
double eps4 = 0.;
#endif

int nb_cell_NB;
double  NB_width ;    // length of the NB

int itrecons;

scalar curve[];


#define Pi 3.14159265358979323846

double geometry(double x, double y, coord center, double Radius) 
{
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  return ( sqrt(R2) - Radius);
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

double max_output = 0.;
double y_max_pre = 10.;


int main() {
  L0 = 800.;
  CFL = 0.5;
  TOLERANCE = 1.e-6;
  NITERMIN = 4;
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  TL[top]    = dirichlet(TL_inf); 
  TL[bottom] = dirichlet(TL_inf); 
  TL[left]   = dirichlet(TL_inf); 
  TL[right]  = dirichlet(TL_inf); 
  
  TS[top]    = dirichlet(TS_inf); 
  TS[bottom] = dirichlet(TS_inf); 
  TS[left]   = dirichlet(TS_inf); 
  TS[right]  = dirichlet(TS_inf); 

  origin(-L0/2., -L0/2.);
  init_grid (1 << 6);
  run();
}


event init(t=0){

  lambda[0]  = 1.;
  lambda[1]  = 1.;
  // DT2        = 0.5;  // diffusion time scale
  DT2        = 0.1*sq(L0/(1<<MAXLEVEL))/lambda[0]; // temperature equilibrium
  nb_cell_NB = 1 << 3 ; // number of cell in the 
                        // narrow band 
  itrecons = 20;
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  coord center = {0.,0.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = 10.1;
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
    itrecons = 10,tolrecons = 1.e-12,NB_width);
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

event stability(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1];
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = itrecons,tolrecons = 1.e-4,
  NB_width);
  double DT3 = 0.5*timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
  // dump();
  // exit(1);
  // fprintf(stderr, "## %g %g %g\n", t, DT2, DT3);
}

event tracer_diffusion(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1];
  advection_LS(
  dist,
  latent_heat,
  cs,fs,
  TS,TL,
  T_eq,
  vpc,vpcf,
  lambda1,lambda2,
  epsK,epsV,eps4,
  curve,
  &k_loop,
  deltat = 0.45*L0 / (1 << grid->maxdepth),
  itredist = 2*nb_cell_NB,
  tolredist = 3.e-3,
  itrecons = itrecons,
  tolrecons = 1.e-4,
  s_clean = 1.e-10,
  NB_width,
  0);
  foreach_face(){ // for GHIGO
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


event interface2(t+=180,last; t<12000){
  output_facets(cs,stdout);
}

event tip_velocity(i++,last){
  // event to output the velocity on the top tip of the dendrite and its
 // position
  double y_max=0, velomax = 0.;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs) && y >0){
      double yy = y+Delta*height(h.y[]);
      if(yy < 1.e10){
        y_max = max(y_max,yy);
        velomax = vpcf.y[];
      }
    }     
  }
  double d0 = epsK;
  fprintf(stderr, "%g %g %g\n", t/sq(d0), (y_max-y_max_pre)/dt,
   velomax*d0);
  y_max_pre = y_max;
}

event movie(t+=200,last){
  squares("TL");
  save("TL.mp4");
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
    (double[]){1.e-3,1.e-3,1.e-3,1.e-3},MAXLEVEL, MINLEVEL);

/* Small on-the-fly outputs*/
  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,2.);
  }
  fprintf(stderr, "##nbcells %d solid volume %g time %g\n", n , solid_volume, t);
}
#endif
/**

![Animation of the approximated temperature field](cube/visu.mp4)(loop)


TODO: 
Add biblio!

*/

