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



#if GHIGO 
#include "../../ghigo/src/myembed.h"
#include "../../ghigo/src/mydiffusion.h"
#else
#include "embed.h"
#include "diffusion.h"
#endif
#include "advection.h"
#include "tracer.h"
#include "curvature.h"

#include "view.h"

#define T_eq         0.
#define TL_inf       1.


/**
Setup of the numerical parameters
*/
int MAXLEVEL, MINLEVEL;
double  H0;


/**
Setup of the physical parameters + level_set variables
*/
scalar TL[];
face vector muv[];

scalar * tracers = {TL};

TL[top]    = dirichlet(TL_inf); 
TL[bottom] = dirichlet(TL_inf); 
TL[left]   = dirichlet(TL_inf); 
TL[right]  = dirichlet(TL_inf); 

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
  run();
}


event init(t=0){
  DT        = L0/(1<<MAXLEVEL);  // Delta

  vertex scalar dist[];
  foreach_vertex() {
    dist[] = -geometry(x,y,L0/4.);
  }

  boundary ({dist});
  restriction({dist});
  fractions (dist, cs, fs);
  fractions_cleanup(cs,fs,smin = 1.e-10);
  boundary({cs,fs});
  restriction({cs,fs});

#if GHIGO
  csm1 = cs;
#endif // GHIGO
  
  foreach() {
    TL[] = TL_inf;
  }
  boundary({TL});
  restriction({TL});

  TL[embed] = dirichlet(T_eq);
#if GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
#endif
}

event properties(i++){
  foreach_face()
    muv.x[] = fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion (i++,t < 2.) {
  boundary({TL});
  diffusion (TL,dt,muv,theta = cm);
  fprintf(stderr, "%g %g\n", dt, DT);
}

event final(t=2.)
  fprintf(stderr, "FIN\n");
