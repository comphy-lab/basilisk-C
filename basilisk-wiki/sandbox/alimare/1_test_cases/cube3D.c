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
plot 'out' w l lw 3 t 'Interface' 
~~~


*/
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define CURVE_LS 0
#define ANISO 1

#include "grid/octree.h"
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
#define TL_inf       -1.
#define TS_inf       -1.

/**
Setup of the numerical parameters
*/
int MAXLEVEL = 7;
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
double  epsK = 0.001, epsV = 0.;
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

scalar curve[];


#define Pi 3.14159265358979323846

double geometry(double x, double y, double z, coord center, double Radius,
 double theta0,
  double A) 
{

  double R2  =  sq(x - center.x) + sq (y - center.y) + sq (z - center.z);
  return sqrt(R2) - Radius;

  // double theta = atan2 (y-center.y, x-center.x);
  // double s = ( sqrt(R2)*(1.-A*cos(4*theta+ theta0)) - Radius);
  // return s - Radius;
}

double TL_init(double val,double V, double t){
  if(val>V*t)
    return TL_inf+exp(-V*(val-V*t));
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
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  TL[top]    = dirichlet(TL_inf); 
  TL[bottom] = dirichlet(TL_inf); 
  TL[front]  = dirichlet(TL_inf); 
  TL[back]   = dirichlet(TL_inf); 
  TL[left]   = dirichlet(TL_inf); 
  TL[right]  = dirichlet(TL_inf); 
  
  TS[top]    = dirichlet(TS_inf); 
  TS[bottom] = dirichlet(TS_inf); 
  TS[front]  = dirichlet(TS_inf); 
  TS[back]   = dirichlet(TS_inf); 
  TS[left]   = dirichlet(TS_inf); 
  TS[right]  = dirichlet(TS_inf); 

  origin(-L0/2., -L0/2., -L0/2.);
  init_grid (1 << (MAXLEVEL-4));
  run();
}


event init(t=0){

  lambda[0]  = 1.;
  lambda[1]  = 1.;
  // DT2        = 0.9*sq(L0/(1<<MAXLEVEL))/lambda[0];  // diffusion time scale
  DT2        = 1.e-4;  // diffusion time scale
  // DT2        = 1.e-5;  // diffusion time scale
  nb_cell_NB = 1 << 3 ; // number of cell in the 
                        // narrow band 

  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  // coord center1 = {0.,0.};
  // coord center = {0.,0.,0.};

  double spacing = 15.;
  coord center1 = {-L0/spacing,L0/spacing,0.};
  coord center2 = {L0/spacing,-L0/spacing,0.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = 1.5*L0/(1 << (MAXLEVEL -3));
  // double theta = 0.;
  double theta0 = 0., theta1 = 60.*Pi/180.;
  double A = 0.3;
  vertex scalar distn[];

  for (int k = 0; k < 5; k++){
    foreach() {
    dist[] = min(geometry(x,y,z,center1,size,theta0,A),
      geometry(x,y,z,center2,size,theta1,A));
      // dist[] = geometry(x,y,z,center,size,theta,A);
    }
    boundary({dist});
    restriction({dist});
    cell2node(dist,distn);

    fractions (distn, cs, fs);

    boundary({cs,fs});
    restriction({cs,fs});
    foreach() {
      TL[] = TL_init(dist[], 60, 0.);
    }
    boundary({dist});
    restriction({dist});

    refine (level < MAXLEVEL && fabs(dist[]) < NB_width);
    adapt_wavelet ({cs,TL},
      (double[]){1.e-3,3.e-3},MAXLEVEL, MINLEVEL);
  }

  foreach() {
    // dist[] = geometry(x,y,z,center,size,theta,A);
    dist[] = min(geometry(x,y,z,center1,size,theta0,A),
      geometry(x,y,z,center2,size,theta1,A));
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
    TL[] = TL_init(dist[], 60, 0.);
    TS[] = 0.;
  }
  boundary({TL,TS});
  restriction({TL,TS});

  myprop(muv,fs,lambda[0]);
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  u.r[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
  uf.r[embed] = dirichlet (0.);
}

event velocity(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = 10,tolrecons = 1.e-12
  );
  double DT3 = timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
  // fprintf(stderr, "## %g %g %g\n", t, DT2, DT3);
}

event tracer_diffusion(i++){
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
  itrecons = 30,
  tolrecons = 1.e-12,
  s_clean = 1.e-10,
  NB_width);

  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});

  // stats s = statsf(vpc.z);
  // norm n = normf(vpc.z);
  // fprintf(stderr, "##verification symetrie vpcz,  avg %g ,min-max %g\n",
  //  n.avg,s.min-s.max);

}

// event backup(t+=6.e-4,last){
//   char name[80];
//   sprintf (name, "snapshot%.4f", t);
//   if(t>0){
//     dump(name);
//     // exit(1);
//   }
// }

event interface2(t+=5.e-4,last){
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});

  int count = 0;
  view(quat = {0.,0.,0.,1.});
  char name[80];
  sprintf (name, "temperature-plane%d-time%.5f.png", ++count,t);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf));
  save(name);

  view(quat = {-0.707107,-0,-0,0.707107});
  sprintf (name, "temperature-plane%d-time%.5f.png", ++count,t);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf), n = {0,1,0});
  save(name);
  
  view (quat = {0.322738,0.470285,0.689538,0.446324});
  sprintf (name, "temperature-plane%d-time%.5f.png", ++count,t);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf), n = {0.5,0.5,0});
  save(name);
}

event movies (t+=1.e-4; t<3.)
{
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  stats s3 = statsf(visu);
  fprintf(stderr, "#temperature %g %g\n",  s3.min, s3.max);  

  char name[80];
  int count = 0;
  view(quat = {0.,0.,0.,1.});
  sprintf (name, "temperature-plane%d-time.mp4", ++count);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf));
  save(name);

  view(quat = {-0.707107,-0,-0,0.707107});
  sprintf (name, "temperature-plane%d-time.mp4", ++count);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf), n = {0,1,0});
  save(name);
  view (quat = {0.322738,0.470285,0.689538,0.446324});
  sprintf (name, "temperature-plane%d-time.mp4", ++count);
  draw_vof("cs");
  squares("visu", min = TL_inf , max = 0.01*fabs(TL_inf), n = {0.5,0.5,0});
  save(name);
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

  adapt_wavelet ({cs,visu,vpcf.x, vpcf.y, vpcf.z},
    (double[]){1.e-3,3.e-3,1.e-2,1.e-2,1.e-2},MAXLEVEL, MINLEVEL);

/* Small on-the-fly outputs*/
  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,3.);
  }
  fprintf(stderr, "%g %d %g\n",t, n , solid_volume);
}
#endif
/**

![Animation of the approximated temperature field](cube/visu.mp4)(loop)


*/

