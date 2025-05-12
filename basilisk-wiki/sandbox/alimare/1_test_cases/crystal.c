/**
#Melt of a solid particle

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The boundaries of the domain are set at $T_L = 1$. The ice particle is initially
at $T_S = -1$.

The temperature on the interface is derived from the Gibbs-Thomson relation,
withtout any anisotropic effect :
$$
 T_{\Gamma} = T_{m} - \epsilon_{\kappa} * \kappa - \epsilon_{v} v_{pc}
$$

![Animation of the temperature field](crystal/temperature.mp4)(loop)
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


#define TL_inf       1.
#define TS_inf       -1.

/**
Setup of the numerical parameters
*/
int MAXLEVEL, MINLEVEL;
double  H0;

/**
Initial geometry definition. Here the interface equation is :

$$
r\left(1+ 0.2 *cos(8\theta) \right) - R
$$
where $r = \sqrt{x^2 + y^2}$, $\theta = \arctan(x,y)$ and $R = \frac{L_0}{5}$

Notice that the initial dist[] field is not really a distance, it is modified
after a few iterations of the LS_reinit() function.
*/
double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.5;
  center.y = 0.5;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = -( sqrt(R2)*(1.+0.2*cos(6*theta)) - Radius);
  // double s = -( sqrt(R2)*(1.+0.*cos(6*theta)) - Radius);


  return s;
}

double TL_init(double val, double factor){
    return T_eq+(TL_inf-T_eq)*tanh(val*L0/factor);
}

/**
$k_{loop}$ is a variable used when the interface arrives in a new cell. Because
we are using embedded boundaries, the solver needs to converge a bit more when
the interface is moved to a new cell. Therefore, the Poisson solver is used
$40*k_{loop} +1$ times.
*/
int k_loop = 0;

int main() {
  CFL = 0.5;
  MAXLEVEL  = 7;
  MINLEVEL = 4;
  int N     = 1 << MAXLEVEL;
  init_grid (N);
  TOLERANCE = 1.e-6;
  NITERMIN = 4;
  TL.third = true;
  TS.third = true;  
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

  run();
}


event init(t=0){
  DT   = sq(L0/(1 << MAXLEVEL));  // Delta
  epsK = 0.001 ; epsV = 0.001;
  eps4 = 0.;
  DT_LS = 0.45*(L0)/(1 << MAXLEVEL);
  itrecons = 60;


  nb_cell_NB = 8 ;               // number of cell in the 
                                      // narrow band 
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach() {
    dist[] = clamp(-geometry(x,y,L0/4.),-1.5*NB_width,1.5*NB_width);
  }

  boundary ({dist});
  restriction({dist});

  LS2fractions(dist,cs,fs);
  
  curvature(cs,curve);

  boundary({curve});
  foreach() {
    TL[] = cs[] > 0. ? TL_init(dist[],0.05) : 0.;
    TS[] = TS_inf;
  }
  boundary({TL,TS});
  restriction({TL,TS});

#if GHIGO
  foreach()
    csm1[]= cs[];

  boundary({csm1});
  restriction({csm1});
#endif // GHIGO

  myprop(muv,fs,lambda[0]);
#if GHIGO // mandatory for GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
#endif

}

event snapshot(t+=0.012){
  output_facets (cs, stdout);
}

event movies ( t+=0.0024,last;t<0.12)
{
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  view (tx = -0.489668, ty = -0.453663);
  squares("visu", min = -1 , max = 1.);
  draw_vof("cs");
  save("temperature.mp4");
}

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
  if(i%10==0.)fprintf(stderr, "%g %g %g\n",t, s2.min, s2.max);
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});

  adapt_wavelet ({cs,visu, vpcf.x, vpcf.y},
    (double[]){1.e-3,1.e-3,5.e-3,5.e-3},MAXLEVEL, MINLEVEL);
  foreach(reduction(+:n)){
    n++;
  }
  curvature(cs,curve);
  boundary({curve});
  if(i%10==0) fprintf(stderr, "##nb cells %d\n", n);
}
#endif

/**

~~~gnuplot Curvature
set key left
set yrange [-30:30]
set xrange [0:0.095]
plot 'log' u 1:2 w l t 'curvature min',  \
     'log' u 1:3 w l  t 'curvature max'
~~~

~~~gnuplot Evolution of the interface (zoom)
set size ratio -1
set yrange [*:*]
set xrange [*:*]
set key top right
unset xlabel
unset xtics
unset ytics
unset border
unset margin
unset ylabel
plot 'out' w l lw 2 t 'interface'
~~~
*/