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

These figures can be produced using the *crystal_growth_gridX.tst* cases defined in my Makefile.

~~~gnuplot Evolution of the interface Maxlevel 6
set key top right
set size ratio -1
set xrange [-2:2]
unset border
unset xlabel
unset xtics
unset ytics
unset margin
unset ylabel
plot 'out_growth_grid6' w l lw 3 lc rgb "black" t ''
~~~

~~~gnuplot Evolution of the interface Maxlevel 7
set key top right
set size ratio -1
set xrange [-2:2]
unset border
unset xlabel
unset xtics
unset ytics
unset margin
unset ylabel
plot 'out_growth_grid7' w l lw 3 lc rgb "black" t ''
~~~

~~~gnuplot Evolution of the interface Maxlevel 8
set key top right
set size ratio -1
set xrange [-2:2]
unset border
unset xlabel
unset xtics
unset ytics
unset margin
unset ylabel
plot 'out_growth_grid8' w l lw 3 lc rgb "black" t ''
~~~

*/
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define CURVE_LS 1
#define DEBUG 0
#define Pi 3.14159265358979323846
#define ANISO 0

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
#define TL_inf       -0.5
#define TS_inf       0.

/**
Setup of the numerical parameters
*/
#ifdef GRIDLEVEL
int MAXLEVEL = GRIDLEVEL; 
int MINLEVEL = GRIDLEVEL-2;
#else
int MAXLEVEL = 7; 
int MINLEVEL = 5;
#endif
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
#ifdef EPSK
double epsK = EPSK;
#else
double epsK = 2.e-3;
#endif
#ifdef EPSV
double epsV = EPSV;
#else
double epsV = -2.e-3;
#endif

#else // Gibbs_Thomson
double  epsK = 0.000, epsV = 0.000;
#endif // Gibbs_Thomson

#define GT_aniso 0
/**
We take into account an [eps4tropic](../phase_change_velocity.h) effect.
We define 2 functions to take into account the eps4tropy in the change velocity
where we set the following condition :
$$
\epsilon = \overline{\epsilon} \left( 1. - \alpha \cos(n \theta +
\theta_0)\right)
$$
where $\theta$ is the angle between the interface and the x-axis, $\alpha =0.5$
and $\theta_0 = 0$ are hardcoded. This will be used in the Gibbs-Thomson
formulation.
*/
#if GT_aniso
double eps4 = 0.1;
#else
double eps4 = 0.;
#endif

int itrecons;

int nb_cell_NB;
double  NB_width ;    // length of the NB

scalar curve[];


#define Pi 3.14159265358979323846


/**
The initial geometry is a crystal seed:

$$
r\left(1+ - 0.25 *cos(4\theta) \right) - R
$$
*/

double geometry(double x, double y, coord center, double Radius) {

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = ( sqrt(R2)*(1.-0.3*cos(4*theta)) - Radius);

  return s;
}

double TL_init(double val,double V, double t){
  if(val>V*t)
    return -1+exp(-V*(val-V*t));
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
  TOLERANCE = 1.e-8;
  NITERMIN  = 4;
  L0 = 4.;
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  origin(-L0/2., -L0/2.);
  init_grid (1 << MAXLEVEL);
  run();
}


event init(t=0){

  lambda[0]  = 1.;
  lambda[1]  = 1.;
  DT2        = 5.e-4;  // diffusion time scale
  // DT2        = 0.3*sq(L0/(1<<MAXLEVEL))/lambda[0];
  nb_cell_NB = 1 << 3 ; // number of cell in the 
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

  foreach()
    foreach_dimension()
      vpc.x[] = 0.;
  boundary((scalar *){vpc});

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

#if CURVE_LS
  curvature_LS(dist, curve);
#else
  curvature(cs,curve);
  boundary({curve});
#endif
    double lambda1 = lambda[0], lambda2 = lambda[1]; 
    LS_speed(
    dist,latent_heat,cs,fs,TS,TL,T_eq,
    vpc,vpcf,lambda1,lambda2,
    epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
    itrecons = 30,tolrecons = 1.e-12);
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
  itrecons = 60,tolrecons = 1.e-12
  );
  double DT3 = timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
  // fprintf(stderr, "## %g %g %g\n", t, DT2, DT3);
}

event tracer_diffusion(i++){
  fprintf(stderr, "## DIFFUSION\n" );
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
event LS_advection(i++,last){
  fprintf(stderr, "## LS_ADVECTION\n" );
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
  itredist = 10,
  tolredist = 3.e-3,
  itrecons = 60,
  tolrecons = 1.e-12,
  s_clean = 1.e-10,
  NB_width);

  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});  
}

event interface2(t+=0.1,last; t<0.8){
  output_facets (cs, stdout);
}



event final(t=end){
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  dump();
}

/**
Mesh adaptation. still some problems with coarsening inside the embedded
boundary.
*/
event adapt (i++, last) {
  fprintf(stderr, "##ADAPT\n");
  scalar normvpc[];

  foreach(){
    if(fabs(dist[]< NB_width)){
      coord temp;
      foreach_dimension()
        temp.x = vpcf.x[];

      normvpc[] = norm2(temp);
    }
    else{
      normvpc[] = 0.;
    }
  }

  boundary({normvpc});
  restriction({normvpc});

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

  adapt_wavelet ({cs,visu,normvpc},
    (double[]){1.e-3,1.e-3,1.e-2},MAXLEVEL, MINLEVEL);

/* Small on-the-fly outputs*/
  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,2.);
  }
  fprintf(stderr, "##time %g dt %g nbcells %d solid volume %g\n", t, dt, n ,
   solid_volume);
  stats s3 = statsf(visu);
  if((s3.min< -10) ||  (s3.max> 10)){
    fprintf(stderr, "#temperature %g %g\n",  s3.min, s3.max);  
    dump();
    exit(1);
  }

}

#if 0
event movies (t+=0.01,last)
{
  fprintf(stderr, "##MOVIES\n");
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
// 
  char filename [100];
  sprintf(filename,"temperature_grid_%d.mp4",MAXLEVEL);
  view (width = 800, height = 800);
  draw_vof("cs");
  squares("visu", min = -0.2 , max = 0.01);
  save (filename);

//   draw_vof("cs");
//   squares("vpc.x", min = 1. , max = 3.);
//   save("vpcx.mp4");
}
#endif
/**

*/

