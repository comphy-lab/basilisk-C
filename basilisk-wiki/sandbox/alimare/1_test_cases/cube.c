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

#define TL_inf       -1.
#define TS_inf       -1.

/**
Setup of the numerical parameters
*/
int MAXLEVEL = 6; 
int MINLEVEL = 4;
double H0;

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

int main() {
  TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));
  TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, eps4));

  TL[top]    = dirichlet(TL_inf); 
  TL[left]   = dirichlet(TL_inf); 
  TL[bottom] = dirichlet(TL_inf); 
  TL[right]  = dirichlet(TL_inf); 

  TS[top]    = dirichlet(TS_inf); 
  TS[bottom] = dirichlet(TS_inf); 
  TS[left]   = dirichlet(TS_inf); 
  TS[right]  = dirichlet(TS_inf); 

  origin(-L0/2., -L0/2.);
  init_grid (1 << MAXLEVEL);
  run();
}


event init(t=0){

  DT         = 5.*sq(L0/(1<<MAXLEVEL));  // Delta
  epsK       = 0.0001 ; 
  epsV       = 0.;
  eps4       = 0.;
  nb_cell_NB = 1 << 3 ; // number of cell in the 
  // narrow band 
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0]  = 0.1;
  lambda[1]  = 0.1;

  coord center1 = {0.,0.};
  // coord center1 = {-L0/20.,L0/20.};
  // coord center2 = {L0/20.,-L0/20.};

/**
And because, why not, we put two crystal seeds, to see them compete during
crystal growth.
*/
  double size = L0/14.99;

  foreach() {
    dist[] = geometry(x,y,center1,size);
  }

  boundary ({dist});
  restriction({dist});

  LS2fractions(dist,cs,fs);

  curvature(cs,curve);
  boundary({curve});

  foreach() {
    TL[] = TL_init(dist[], 60, 0.);
    TS[] = 0.;
  }
  boundary({TL,TS});
  restriction({TL,TS});

  myprop(muv,fs,lambda[0]);
#if GHIGO // mandatory for GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
#endif
}

event interface2(t+=0.1,last){
  output_facets (cs, stdout);
}

event movies (t+=0.01; t<1.)
{
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  stats s3 = statsf(visu);
  fprintf(stderr, "#temperature %g %g\n",  s3.min, s3.max);  
  view (width = 800, height = 800);
// 
  char filename [100];
  sprintf(filename,"temperature_epsK%g_epsV%g_eps4%g.mp4",epsK,epsV,eps4);

  draw_vof("cs");
  squares("visu", min = -1. , max = 0.01);
  save (filename);
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

  boundary({cs,cs2});
  fractions_cleanup(cs,fs,smin = 1.e-10);
  fractions_cleanup(cs2,fs2,smin = 1.e-10);
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

  adapt_wavelet ({cs,visu},
    (double[]){1.e-3,1.e-4},MAXLEVEL, MINLEVEL);
  foreach(reduction(+:n)){
    n++;
  }
  fprintf(stderr, "##nb cells %d\n", n);
}
#endif
/**

![Animation of the approximated temperature field](cube/visu.mp4)(loop)


*/

