/**

#Solidifying planar interface 

Exact solution.
Speed should be constant throughout the simulation.
Temperature profile :

$$
T(x,t) = \left\{\begin{array}{cc}
-1 + e^{-V(x-Vt)}& , x > Vt\\
0               & , x \leq Vt\\
\end{array}\right.
$$

Plot of the $L_1$-error.

~~~gnuplot Position of the interface
f(x) = x
set xrange [0:0.2]
set key left
plot 'log0' u 1:2 w l t '32x32   VOF height',\
     'log1' u 1:2 w l t '64x64   VOF height',\
     'log2' u 1:2 w l t '128x128 VOF height',\
     'log0' u 1:3 w p pt 7 lc 'blue' t '32x32   LS',\
     'log1' u 1:3 every 2 w p pt 7 lc 'green' t '64x64   LS',\
     'log2' u 1:3 every 4 w p pt 7 lc 'red' t '128x128 LS'
~~~

~~~gnuplot Average error convergence CFL = 0.5
unset xrange
unset yrange

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)

f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b

f2(x) = a2 + b2*x
fit f2(x) 'log' u (log($1)):(log($4)) via a2,b2

f3(x) = a3 + b3*x
fit f3(x) 'log' u (log($1)):(log($6)) via a3,b3

set ylabel 'Average error'
set xrange [16:256]
set yrange [*:*]
set xtics 16,2,256
set format y "%.1e"
set logscale
plot '' u 1:2 pt 7 lc 'blue' t 'all cells', exp(f(log(x))) t ftitle(a,b) lc 'blue', \
     '' u 1:4 pt 7 lc 'green' t 'partial cells', exp(f2(log(x))) t ftitle(a2,b2) lc 'green', \
     '' u 1:6 pt 7 lc 'red' t 'full cells', exp(f3(log(x))) t ftitle(a3,b3) lc 'red'
~~~

~~~gnuplot Initial and final temperature plots
set term svg size 1000,1000
set key right
unset logscale
unset xrange
set ylabel 'Temperature profile' 
set xrange[0:0.5]
plot 'out0' u 1:2 w p pt 7 lc 'blue' t 'Final temperature 32x32  ',\
     'out1' u 1:2 w p pt 7 lc 'green' t '64x64  ',\
     'out2' u 1:2 w p pt 7 lc 'red' t '128x128', \
     '' u 1:3 w l lt -1 t 'Reference final',\
     'init0' u 1:2 w p pt 7 lc 'blue' t 'Initial temperature 32x32 ',\
     'init1' u 1:2 w p pt 7 lc 'green' t '64x64  ',\
     'init2' u 1:2 w p pt 7 lc 'red' t '128x128',\
     '' u 1:3 w l lt -1 t 'Reference initial'
~~~

~~~gnuplot Final error plots
set term svg size 600,600
set key left
set logscale y
set ylabel 'Error profile'
set xrange [0.2:0.5]
plot 'out0' u 1:4 w p pt 7 lc 'blue' t '32x32',\
     'out1' u 1:4 w p pt 7 lc 'green' t '64x64',\
     'out2' u 1:4 w p pt 7 lc 'red' t '128x128'
~~~


*/

#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0
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
#include "../LS_speed.h"
// #include "../simple_discretization.h"
#include "../LS_reinit.h"

#include "view.h"


#define T_eq          0.
double TL_inf(double x,double V, double t){
  if(x>V*t)
    return -1+exp(-V*(x-V*t));
  else
    return 0.;
}

double TL_interf(double x,double V, double t){
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


scalar * tracers   = {TL};
scalar * tracers2  = {TS};
scalar * level_set = {dist};
face vector muv[];
scalar grad1[], grad2[];
double DT2;

double  epsK = 0., epsV = 0.;



double lambda[2];

#define GT_aniso 0
int aniso = 1;

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

double s_clean = 1.e-10; // used for fraction cleaning

TL[embed] = dirichlet(T_eq);
TS[embed] = dirichlet(T_eq);

/**
The inteface moves at constant speed 1.
*/
#define V  1.
#define t0 1.e-6
TL[top]    = dirichlet(TL_inf(y,V,t-t0)); 
TS[bottom] = dirichlet(TS_inf); 


int j;
int k_loop = 0;


/**
*/

#define geo(x, V,t) (x-V*t)


char filename [100];
FILE * fp1, * fp2, * fp3;


/**
*/

void myprop(face vector muv, face vector fs, 
  double lambda){
  foreach_face()
  muv.x[] = lambda*fs.x[];
  boundary((scalar *) {muv});
}


int main() {
  L0 = 1.;
  CFL = 0.5;
  origin (-0.5*L0, -0.5*L0);
  periodic(left);  
  j = 1;
  for (j=0;j<=2;j++){
    MAXLEVEL = 5+j, MINLEVEL = 5 ;
    double dtref = 1.e-2; // < 0.5*L0/(1<< MAXLEVEL);
    DT2 = dtref/powf(4,j); 
    DT  = DT2;
    fprintf(stderr, "##DT %g\n", DT2);
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");
    snprintf(filename, 100,  "out%d", j);
    fp2 = fopen (filename,"w");
    snprintf(filename, 100,  "init%d", j);
    fp3 = fopen (filename,"w");
/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 1;
    TL.third = true;
    TS.third = true; 
    H0 = 0.5*L0; 
    N = 1 << MAXLEVEL;
    init_grid (N);
    run(); 
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
  }
}

/**
*/

event init(t=0){

  TOLERANCE = 1.e-9;
  NITERMIN = 4;

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);
  // NB_width = 1.;

  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach(){
    dist[] = geo(y,V,t-t0);
  }
  boundary ({dist});
  restriction({dist});


  vertex scalar distn[];
  foreach_vertex(){
    if((point.j-1) > (1 << grid-> depth)){
      distn[] = dist[] +Delta/2.;
    }
    else{
      distn[] = dist[] - Delta/2.; 
    }
  }
  boundary({distn});
  restriction({distn});

  fractions (distn, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  foreach(){
    if(cs[]>0.){
      TL[] = TL_interf(y,V,t);
    }
    else{
      TL[]=  nodata;
    }
    TS[] = 0.;
    boundary({TS,TL});
  }

  boundary({TL,TS});
  restriction({TL,TS});

  foreach(){
    if(x==Delta/2. && TL[] != nodata){
      fprintf(fp3, "%.8g %.8g %.8g %.8g\n", y, TL[], TL_interf(y,V,t-t0),
       fabs(TL[]-TL_interf(y,V,t-t0)));}
  }
  int n =0;
  foreach(reduction(+:n)){
    n++;
  }
  myprop(muv,fs,lambda[0]);

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

event velocity(i++){
  tnext= t + DT2;
  dt = DT2;
}

event tracer_diffusion(i++){
  boundary({TL});
  myprop(muv,fs,lambda[0]);
  diffusion (TL,dt,muv);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[1]);
  boundary({TS});
  diffusion (TS,dt,muv);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[0]);
}

/**
Note about the events order:\\
For this test case, the best order of accuracy is obtained when the LS_advection
is done after de the diffusion (it's the opposite for the case where the phase
change velocity is not calculated but imposed see [here]
(#mydirichlet_expand_1D.c))
*/
event LS_advection(i++){
/**
Previous state of the metric is saved.
*/
  foreach(){
    csm1[]   = cs[];
  }

  boundary({csm1});
  restriction({csm1});

  vector vpc[];
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, T_eq, vpc, latent_heat, 
      lambda,epsK=0, epsV =0, aniso =1);/**
We copy this value in vpcr.
*/
  vector vpcr[];
  foreach(){
    foreach_dimension(){
      if(interfacial(point,cs))vpcr.x[] = vpc.x[];
      else vpcr.x[] = 0.;
    }
  }

  boundary((scalar * ){vpcr});
  restriction((scalar * ){vpcr});

  scalar * speed_recons  = {vpcr.x,vpcr.y};

/**
We reconstruct a cell-centered vpc field in the vicinity of the interface where
vpcr is the solution to the bilinear interpolation of vpc on face centroids.
*/
  int err = 0;

  recons_speed(dist, dt = 0.5 * L0/(1<<MAXLEVEL), LS_speed = speed_recons,
   tolerance = 1.e-12, &err, 
   nb_iter = 60, 
   cs, fs);

  RK3(dist,vpcr,dt, NB_width);

  boundary ({dist});
  restriction({dist});
  LS_reinit(dist, it_max = 10);

  vertex scalar distn[];
  foreach_vertex(){
    if((point.j-1) > (1 << grid-> depth)){
      distn[] = dist[] +Delta/2.;
    }
    else{
      distn[] = dist[] - Delta/2.; 
    }
  }
  boundary({distn});
  restriction({distn});
  fractions (distn, cs, fs);
  fractions_cleanup(cs,fs,smin = s_clean);

  boundary({cs,fs});
  restriction({cs,fs});


  foreach() {
    if (cs[] > 0. && csm1[] <= 0.) { // Emerged cells
      fprintf(stderr, "%g\n", y);
      // Barycenter and normal of the embedded fragment
      coord b, n;
      embed_geometry (point, &b, &n);

      // Boundary condition at time t
      bool dirichlet = true;
      double ab = (TL.boundary[embed] (point, point, TL, &dirichlet));
      assert (dirichlet);

      TL[] = dirichlet_extrapolate (point, TL, n, b, ab, dirichlet);
    }
    if(cs[] <=0. && csm1[]>0.){
      TL[] = nodata;
    }
  }
  /**
  We update the boundary condition for *TL*. */
  boundary({TL});
  restriction({TL});

  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});
}


event movies ( i++,last;t<=0.2)
{
  double y_max=0,y_LS = 0.;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max) reduction(max:y_LS)){
    if(interfacial(point, cs) && y >0){
      double yy = y+Delta*height(h.y[]);
      if(yy < 1.e10)y_max = max(y_max,yy);
      y_LS = max(y_LS,y-dist[]);
    }     
  }

/**
We've already done the advection of the interface, so the time is $t+dt$.
*/
  fprintf(fp1, "%g %g %g\n", t+dt, y_max,y_LS);
}

event final_error(t=end,last){
  scalar e[], ep[], ef[];
    foreach() {
    if (cs[] == 0.)
      e[] = nodata;
    else {
    if(cs[] <1.){
        e[] = TL[] - TL_interf( y,V,t-t0+dt);
      }
      else{
        e[] = TL[] - TL_inf( y,V,t-t0+dt);
      }
      ep[] = cs[] < 1. ? e[] : nodata;
      ef[] = cs[] >= 1. ? e[] : nodata;
    }
  }
  foreach(){
    if(x==Delta/2. && cs[] > 0){
      fprintf(fp2, "%.8g %.8g %.8g %.8g\n", y, TL[], TL_inf(y,V,t-t0),
        fabs(TL[]-TL_inf(y,V,t-t0)));
    }
  }
  norm n = normf (e), np = normf (ep), nf = normf (ef);
  fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %g\n",
   N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, dt);
  fflush (ferr);

  if(j==2){
    squares("e");
    save("e.png");
  }
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
}
#endif

/**


Here we first check that the interface is correctly advected, error on the
position of the interface should follow f(x) = x




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
