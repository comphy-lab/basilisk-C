/**
# Validation of gradient calculation

We set up an initial field  :
$$
T(x,t) = \left\{\begin{array}{cc}
T_{eq} + e^{y0-y} & , y>y_0\\
T_{eq}        & , y \leq y_0\\
\end{array}\right.
$$

We check that the gradient on the interface is 1.
*/

#define DOUBLE_EMBED 1
#define LevelSet     1

#include "embed.h"
#include "../double_embed-tree.h"

#include "advection.h"

#include "fractions.h"
#include "curvature.h"

#include "../level_set.h"
#include "view.h"
#include "../phase_change_velocity.h"

#define T_eq          0.
#define TL_inf        1.
#define TS_inf        0.

int MINLEVEL, MAXLEVEL; 
double latent_heat;

#define T_eq         0.

#define plane( y,H0) (y -H0)
#define Stefan 1.

double T_init(double y,double y0){
  if(y>y0)
    return T_eq + 1.-exp(y0-y);
  else
    return T_eq ;
  
}

scalar TL[], TS[], dist[];
vector vpc[];
face vector vpcf[];


scalar * tracers   = {TL};
scalar * tracers2  = {TS};
scalar * level_set = {dist};
scalar grad1[], grad2[];
scalar curve[];

TL[embed] = dirichlet(T_eq);
TS[embed] = dirichlet(T_eq);

double lambda[2];

double  epsK = 0., epsV = 0.;
int aniso = 1;

// u.n[embed] = dirichlet(vpcfn(point,vpc));
// u.t[embed] = dirichlet(0.);

int j;

int main() {
  periodic(left);

  TL.third = true;
  TS.third = true;

  for (j=0;j<=2;j++){

/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 1/Stefan;
    MAXLEVEL = MINLEVEL = 5+j ;
    N = 1 << MAXLEVEL;
    init_grid (N);
    lambda[0] = 1.;
    lambda[1] = 1.;
    double y0 = 0.1;

/**
not mandatory
*/
    foreach(){
      dist[] = plane(y,y0);
    }
    boundary ({dist});
    restriction({dist});

    vertex scalar dist_n[];
    foreach_vertex(){
      dist_n[] = plane(y,y0);
    }
    boundary({dist_n});
    restriction({dist_n});

    fractions (dist_n, cs, fs);
    fractions_cleanup(cs,fs);

    boundary({cs,fs});
    restriction({cs,fs});

    curvature(cs,curve);

    boundary({curve});
    foreach() {
      TL[] = T_init(y,y0);
      TS[] = T_init(y,y0);
    }

    foreach_face(){
      vpc.x[] = 0.;
    }

    boundary({TL,TS});
    restriction({TL,TS});
    phase_change_velocity_LS_embed (cs, fs ,TL, TS, T_eq, vpc, latent_heat, 
      lambda,epsK=0, epsV =0, aniso =1);
    double myprint = 0;
    foreach_face(y){
      if(vpc.y[] !=0. && myprint == 0){
        myprint  = vpc.y[];
      }
    }
    fprintf(stderr, "%d %g\n", N,fabs(-myprint-Stefan));
    squares("TL");
    save("TL.png");
    if(fabs(myprint+Stefan) > 1.e-3 ){
      fprintf(stderr, "ERROR, %g\n", myprint-Stefan);
      exit(1);
    }
  }
}

/**
~~~gnuplot
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b
set logscale
plot 'log' u 1:2 t 'error on speed' , exp(f(log(x))) t ftitle(a,b)
~~~
*/