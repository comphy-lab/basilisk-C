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
set object 1 circle front at 0.,0. size 1.75 fillcolor rgb "black" lw 1 dt 2
plot 'out' w l lw 0.8 t 'Interface' 
~~~

*/
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define CURVE_LS 1

#include "grid/octree.h"
#include "../../ghigo/src/myembed.h"
#include "../double_embed-tree.h"

#include "../../ghigo/src/mydiffusion.h"
#include "../advection_A.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"

#define dtLS 1
#include "../LS_diffusion.h"

#define TL_inf       -0.8
#define TS_inf       -0.8

/**
Setup of the numerical parameters
*/
int MAXLEVEL = 8;
int MINLEVEL = 4;
double H0;

/**
Initial interface is a circle, we want to study the effect of the anisotropy on
the curvature constant in the Gibbs-Thomson equation, i.e. we want to introduce
dendritic growth purely due to anisotropy in the GT equation.
*/

double geometry(double x, double y, double z, coord center, double Radius) {

  double R2  =  sq(x - center.x) + sq (y - center.y) + sq (z - center.z) ;
  double s = ( sqrt(R2) - Radius);

  return s;
}


double TL_init(double val,double V, double t){
    return TL_inf+exp(-V*(val-V*t));
}


int main() {
  L0 = 4.;
  CFL = 0.3;
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  origin(-L0/2., -L0/2., -L0/2.);
  init_grid (1 << 6);
  run();
}


event init(t=0){

  lambda[0]  = 1.;
  lambda[1]  = 1.;

  epsK = 0.001;
  epsV = 0.;
  eps4 = 0.4;

  DT        = 0.1*sq(L0/(1<<MAXLEVEL))/lambda[0];
  DT_LS = 0.45*(L0)/(1 << MAXLEVEL);
  nb_cell_NB = 1 << 3 ; // number of cell in the 
                        // narrow band 
  itrecons = 30;
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  coord center = {0.,0.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = 0.1;

  vertex scalar distn[];

  for (int k = 0; k < 4; k++){
    foreach() {
      dist[] = geometry(x,y,z,center,size);
    }
    boundary({dist});
    restriction({dist});
    LS2fractions(dist,cs,fs);

    foreach() {
      TL[] = TL_init(dist[],40,0.);
      TS[] = 0.;
    }
    boundary({dist,TL,TS});
    restriction({dist,TL,TS});
    double lambda1 = lambda[0], lambda2 = lambda[1]; 
    LS_speed(
    dist,latent_heat,cs,fs,TS,TL,T_eq,
    vpc,vpcf,lambda1,lambda2,
    epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
    itrecons = 30,tolrecons = 1.e-3,NB_width);
    event("adapt");
    
  }

  foreach() {
    dist[] = geometry(x,y,z,center,size);
  }
  boundary ({dist});
  restriction({dist});
  LS2fractions(dist,cs,fs);

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
    if(cs[]>0.)
      TL[] = TL_init(dist[],40,0.);
    else
      TL[] = 0.;
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

event movies (t+=1.e-4,last)
{
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});

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

  boundary({cs,cs2});
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
    (double[]){1.e-3,1.e-2,1.e-2},MAXLEVEL, MINLEVEL);

/* Small on-the-fly outputs*/
  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,3.);
  }
  fprintf(stderr, "##time %g dt %g nbcells %d solid volume %g\n", t, dt, n ,
   solid_volume);
}
#endif
/**

![Animation of the approximated temperature field](cube/visu.mp4)(loop)


*/

