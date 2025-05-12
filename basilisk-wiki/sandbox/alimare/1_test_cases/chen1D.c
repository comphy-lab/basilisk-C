/**

#Solidifying planar interface 

Exact solution.
Speed should be constant throughout the simulation.
Temperature profile :

$$
T(x,t) = \left\{\begin{array}{cc}
-1 + e^{V(x-Vt)}& , x > Vt\\
0               & , x \leq Vt\\
\end{array}\right.
$$

Plot of the $L_1$-error.

*/

#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0
#define BIQUADRATIC 1
#define GHIGO 1

#if GHIGO
#include "../../ghigo/src/myembed.h"
#else
#include "embed.h"
#endif
#include "../double_embed-tree.h"

#include "../advection_A.h"
#if GHIGO
#include "../../ghigo/src/mydiffusion.h"
#else
#include "diffusion.h"
#endif

#include "fractions.h"
#include "curvature.h"

#include "../level_set.h"
#include "../LS_advection.h"
#include "view.h"


#define T_eq          0.
double TL_inf(double x,double V, double t){
  return -1+exp(-V*(x-V*t));
}
#define TS_inf        0.


int MINLEVEL, MAXLEVEL; 
double H0;
double latent_heat;


/**
Setup of the physical parameters + level_set variables
*/
scalar TL[], TS[], dist[];
vector vpc[];
face vector vpcf[];


scalar * tracers   = {TL};
scalar * tracers2  = {TS};
scalar * level_set = {dist};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];
double DT2;

double  epsK = 0., epsV = 0.;

scalar curve[];


double lambda[2];

#define GT_aniso 0
#if GT_aniso
int aniso = 6;
#else
int aniso = 1;
#endif


int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

double s_clean = 1.e-10; // used for fraction cleaning
  
mgstats mg1,mg2;

#define V  1.
TL[embed] = dirichlet(T_eq);
TS[embed] = dirichlet(T_eq);

TL[top]    = dirichlet(TL_inf(y,V,t-1.e-6)); 
TS[bottom] = dirichlet(TS_inf); 


int j;
int k_loop = 0;


/**
*/

#define geo(x, V,t) (x-V*t)


char filename [100];
FILE * fp1,  * fp2;


/**
*/

void myprop(face vector muv, face vector fs, 
  double lambda){
  foreach_face()
    muv.x[] = lambda*fs.x[];
  boundary((scalar *) {muv});
}

void writefile(mgstats mgd){
  if((mgd.resa > TOLERANCE) ) {
/**
If the calculation crashes (it often does if the Poisson solver does not
converge) we save the last state of the variables
*/
    scalar ffsx[], ffsy[];
    scalar point_f[];
    vector x_interp[];
    foreach(){
      ffsx[] = fs.x[];
      ffsy[] = fs.y[];
      if(interfacial(point,cs)){
        coord n       = facet_normal( point, cs ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (cs[], n);
        line_length_center (n, alpha, &p);
        x_interp.x[] = p.x;
        x_interp.y[] = p.y;
        coord segment[2];
        point_f[] = facets (n, alpha, segment);
      }
      else{
        x_interp.x[] = nodata;
        x_interp.y[] = nodata;
        point_f[]    = 0;
      }
    }
    boundary({ffsx,ffsy});
    vertex scalar distn[];
    cell2node(dist,distn);
    dump();
    fprintf(stderr, "#CRASH#");
    exit(1);
  }
}


int main() {

  L0 = 1.;
  CFL = 0.5;
  origin (-0.5*L0, -0.5*L0);
  periodic(left);  
  j = 1;
  for (j=0;j<=2;j++){
    DT2 = 0.0001/powf(2,j);
    fprintf(stderr, "##DT %g\n", DT2);
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");
    snprintf(filename, 100,  "out%d", j);
    fp2 = fopen (filename,"w");
/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 1;
    MAXLEVEL = 6+j, MINLEVEL = 5 ;
    TL.third = true;
    TS.third = true; 
    H0 = 0.5*L0; 
    N = 1 << MAXLEVEL;
    init_grid (N);
    run(); 
    fclose(fp1);
    fclose(fp2);
  }
  exit(1);
}

/**
*/


event init(t=0){

  TOLERANCE = 1.e-9;
  NITERMIN = 4;

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach(){
      dist[] = geo(y,V,t-1.e-6);
  }
  boundary ({dist});
  restriction({dist});


  vertex scalar dist_n[];
  cell2node(dist,dist_n);

  fractions (dist_n, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  curvature(cs,curve);
  boundary({curve});

  foreach(){
    if(y>V*(t-1.e-6)){
      TL[] = TL_inf(y,V,t-1.e-6);
    }
    else{
      TL[]=  0.;
    }
    TS[] = 0.;
    boundary({TS,TL});
  }

  boundary({TL,TS});
  restriction({TL,TS});
  int n =0;
  foreach(reduction(+:n)){
      n++;
    }
  myprop(muv,fs,lambda[0]);

  double y_max=0;;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs) && y >0){
      double yy = y+Delta*height(h.y[]);
      if(yy < 1.e10)y_max = max(y_max,yy);
    }     
  }
  fprintf(fp1, "%g %g\n", t+dt, y_max);

#if GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
#endif
}

/**
This event is the core the of the hybrid level-set/embedded boundary.
*/


/**
*/
event velocity(i++){
  
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  advection_LS(
  dist,
  latent_heat,
  cs,fs,
  TS,TL,
  T_eq,
  vpc,vpcf,
  lambda1,lambda2,
  epsK,epsV,aniso,
  curve,
  &k_loop,
  deltat = 0.45*L0 / (1 << grid->maxdepth),
  itredist = 10,
  tolredist = 3.e-3,
  itrecons = 60,
  tolrecons = 1.e-6,
  s_clean = 1.e-10,
  NB_width);

  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});
}
/**
*/


event tracer_diffusion(i++){
  boundary({TL});

// #if GHIGO // adding a reaction  in interfacial cells ?
//   scalar r[],Teq[];
//   foreach(){
//     // here we suppose that Tm = 0
//     Teq[] = T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, aniso);
//   }
//   boundary({Teq});
//   restriction({Teq});
    
//   foreach() {

//     double xc = x, yc = y, zc = z;
//     double grad = 0;
//     if (cs[] > 0. && cs[] < 1.) {
//       coord n = facet_normal (point, cs, fs), p;
//       double alpha = plane_alpha (cs[], n);
//       plane_center (n, alpha, cs[], &p); // Different from line_area_center
//       xc += p.x*Delta, yc += p.y*Delta, zc += p.z*Delta;
//       double temp = Teq[];
//       double c    = 0.;
//       double area = plane_area_center (n, alpha, &p);
//       grad += dirichlet_gradient(point, TL, cs , n, p, 
//         temp, &c)*n.x*area;
//     }
//     r[] = grad;
//   }
//   mgT = diffusion (TL,dt,muv,theta = cm, r = r);
// #else
  mgT = diffusion (TL,dt,muv,theta = cm);
// #endif

  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[1]);
  boundary({TS});
  mgT = diffusion (TS,dt,muv,theta = cm);
  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[0]);
}





event movies ( i++,last;t<0.01)
{

  double y_max=0;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs) && y >0){
      double yy = y+Delta*height(h.y[]);
      if(yy < 1.e10)y_max = max(y_max,yy);
    }     
  }
  fprintf(fp1, "%g %g\n", t, y_max);
}

event final_error(t=end){

  scalar e[];
  foreach() {
    if (cs[] == 0.)
      e[] = nodata;
    else {
      e[] = TL[] - TL_inf( y,V,t-1.e-6);
    }
  }
  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g \n",N, n.avg, n.max);
  fflush (ferr);
}

#if 1
event adapt (i++, last) {

  foreach_cell(){
    cs2[] = 1.-cs[];
  }
  foreach_face(){
      fs2.x[]      = 1.-fs.x[];
  }

  boundary({cs,cs2,fs,fs2});
  fractions_cleanup(cs,fs,smin = 1.e-14);
  fractions_cleanup(cs2,fs2,smin = 1.e-14);
  restriction({cs,cs2,fs,fs2});
  int n=0;
  stats s2 = statsf(curve);
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});

  adapt_wavelet ({cs,curve,visu},
    (double[]){1.e-3,1.e-3,5.e-4},MAXLEVEL, MINLEVEL);
  foreach(reduction(+:n)){
    n++;
  }

  curvature(cs,curve);
  boundary({curve});

}
#endif

/**



~~~gnuplot Position of the interface
f(x) = x
set xrange [0:0.01]
plot 'log0' u 1:2 w p t '32x32',\
     'log1' u 1:2 w p t '64x64',\
     'log2' u 1:2 w p t '128x128',\
     f(x) w l 
~~~

~~~gnuplot L1-error
set xrange [32:512:2]
set logscale
plot 'log' u 1:2 w l t 'avg', '' u 1:3 w l t 'max'
~~~


~~~bib

@article{chen_simple_1997,
  title = {A simple level set method for solving {Stefan} problems},
  volume = {135},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Chen, S and Merriman, B and Osher, S and Smereka, P},
  year = {1997},
  pages = {8--29}
}

~~~
*/
