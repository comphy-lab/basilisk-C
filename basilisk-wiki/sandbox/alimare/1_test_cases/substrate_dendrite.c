/**
#Solidification of an undercooled liquid

We want to study the dendritic growth of a solid in an undercooled liquid on a
cold plate. We simulate the diffusion of two tracers separated by an embedded
boundary. The interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
$$

where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The boundaries of the domain are set at $T_L = -1$. The ice particle is initially
at $T_S = 0$. And the substrate's influence is modelled with a neumann boundary
condition

The temperature on the interface is derived from the Gibbs-Thomson relation:
$$
 T_{\Gamma} = T_{m} - \epsilon_{\kappa} * \kappa - \epsilon_{v} v_{pc}
$$

~~~gnuplot Interface every dt=0.01
set term svg size 1000,400
set key top right
set size ratio -1
set xrange [-2:2]
plot 'out' w l lw 3 lc rgb "black" t ''
~~~

~~~gnuplot Height of the top dendrite
unset xrange
set size square
plot 'log' u 1:2 w l t 'Height'
~~~

~~~gnuplot Radius of the ice on the cold plate
unset xrange
unset size
plot 'log' u 1:3 w l t 'Radius'
~~~

*/

#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define CURVE_LS 1

#include "../../ghigo/src/myembed.h"
#include "../double_embed-tree.h"

#include "../../ghigo/src/mydiffusion.h"
#include "../advection_A.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"

#define dtLS 1
#include "../LS_diffusion.h"


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

#define GT_aniso 0

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
  DT = 1.;
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TL[bottom] = neumann(-1.);
  TS[bottom] = neumann(-1.);
  TL.third = true;
  TS.third = true;
 
  origin(-L0/2., -L0/2.);
  init_grid (1 << MAXLEVEL);
  run();
}


event init(t=0){

  epsK = 2.e-3;
  epsV = 2.e-3;
  eps4 = 0.;
  itrecons = 50;
  tolrecons = 1.e-10;
  DT_LS = 0.45*(L0)/( 1 << MAXLEVEL);
  
  nb_cell_NB = 1 << 3 ; // number of cell in the 
                        // narrow band 
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);
  coord center = {0.,-L0/2.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = 0.1;
  for (int k = 0; k < 5; k++){
    foreach() {
      dist[] = geometry(x,y,center,size);
    }
    boundary({dist});
    restriction({dist});
    
    LS2fractions(dist,cs,fs);

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
    itrecons = 30,tolrecons = 1.e-12,NB_width);
    event("adapt");
   
  }

  foreach() {
    dist[] = geometry(x,y,center,size);
  }
  boundary ({dist});
  restriction({dist});

  LS2fractions(dist,cs,fs);

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

/**
This event is the core the of the hybrid level-set/embedded boundary.
*/

event interface2(t+=0.01,last; t<0.13){
  output_facets (cs, stdout);
}

event output(i++,last){
  double y_max=0,x_max= 0.;
  double y_ref = -L0/2.;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max) reduction(max:x_max)){
    if(interfacial(point, cs)){
      double yy = y-y_ref+Delta*height(h.y[]);
      double xx = fabs(x + Delta*height(h.x[]));
      if(yy < 1.e10)y_max = max(y_max,yy);
      if(xx < 1.e10)x_max = max(x_max,xx);
    }     
  }
  fprintf(stderr, "%g %g %g\n", t+dt, y_max, x_max);
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

