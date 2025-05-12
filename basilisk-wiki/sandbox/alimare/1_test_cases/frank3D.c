/**
#Frank's Spheres

This theory of this test case has been studied originally by 
[Frank](#frank1950radially). We start with a solid sphere at melting
temperature surrounded by an undercooled liquid. This case is stable. 

##Analytical Solutions

A class of self similar solutions has been developed by Frank, in one, two and
three dimensions. For $d =1, 2, 3$ the solid interior of a slab, cylinder or
sphere of radius,
$$
R(t) = S t^{1/2}
$$
parameterized by S and the temperature field is 0 if $s <$ and:
$$
T(x,t) = T(s) = T_\infty * (1- \frac{F_d(s)}{F_d(S)})
$$
if s> S, and:
$$
F_2(s) = E_1(\frac{s^2}{4})
$$
We use the GNU Scientific Library for comparison with theoretical values. $E_1$
is readily available in the GSL, using `gsl_sf_expint_E1` function which is :
$$
E_1(x) = \int_1^\infty \frac{e^{-xt}}{t}dt = \int_x^\infty \frac{e^{-t}}{t}dt
$$
*/

/**
Plots:
~~~gnuplot Test
set term pdf enhanced font ",18"
set key top left
set xlabel "t"
set ylabel "Radius"
plot 'log0' u 1:4 pt 5 ps 0.4 lc 'black'  t '32x32',\
     'log1' u 1:4 every 2 pt 6 ps 0.4 lc 'red'    t '64x64',\
     'log2' u 1:4 every 4 pt 7 ps 0.4 lc 'blue' t '128x128',\
      '' u 1:5 w l lt -1 lw 1.5 lc rgb 'red' t 'theo'
~~~

~~~gnuplot Average error convergence CFL = 0.5
ftitle(a,b) = sprintf("%.3fx^\{%4.2f\}", exp(a), -b)

f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b

f2(x) = a2 + b2*x
fit f3(x) 'log' u (log($1)):(log($4)) via a2,b2

f3(x) = a3 + b3*x
fit f3(x) 'log' u (log($1)):(log($6)) via a3,b3

set ylabel 'Average error'
set xrange [16:256]
set yrange [*:*]
set xtics 16,2,256
set format y "%.1e"
set logscale
set key top right
plot '' u 1:2 pt 7 lc 'black' t 'all cells'  , exp(f(log(x)))  lt -1 lw 1.5 t ftitle(a,b)   , \
     '' u 1:4 pt 7 lc 'red' t 'partial cells', exp(f3(log(x))) lt 2 lw 1.5 t ftitle(a2,b2) , \
     '' u 1:6 pt 7 lc 'blue' t 'full cells'  , exp(f3(log(x))) lt 4 lw 1.5 t ftitle(a3,b3)
~~~

![Final error field](frank_spheres/error.png)

*/


#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define CURVE_LS 0
#define DEBUG 0
#define ANISO 0

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
#include "../LS_advection.h"

#define T_eq          0.
#define TL_inf       -1./2
#define TS_inf        0.

#define tstart 0.

int MINLEVEL, MAXLEVEL; 
double latent_heat;

#define DT_MAX  1.



/**
Setup of the physical parameters + level_set variables
*/
scalar TL[], TS[], dist[];
vector vpc[],vpcf[];

scalar * tracers   = {TL};
scalar * tracers2  = {TS};
scalar * level_set = {dist};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];
double DT2;

double  epsK = 0.000, epsV = 0.000;
scalar curve[];


double lambda[2];

double eps4 = 0.;


int     nb_cell_NB =  1 << 2 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

double s_clean = 1.e-10; // used for fraction cleaning

#include "../basic_geom.h"
  
mgstats mg1,mg2;

int j;
int k_loop = 0;

/**
Theoretical functions and parameters for this test case.
Only working in 2D for now.
*/

#include <gsl/gsl_sf_erf.h>
#pragma autolink -lgsl -lgslcblas

double F3(double s){
  return 1/s*exp(-sq(s)/4.)-0.5*(sqrt(Pi)*gsl_sf_erfc(s/2.));
}


double Theo(double r, double t, double undercooling, double F3_S){
  double s = r/pow(t,1./2);
  return undercooling * (1 - F3(s) / F3_S);
}


TL[embed] = dirichlet(T_eq);
TS[embed] = dirichlet(T_eq);

#define radius(x, y, z) (sqrt(sq(x)+sq(y) + sq(z)))
TL[top]    = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));
TL[front]  = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));
TL[back]   = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));
TL[bottom] = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));
TL[left]   = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));
TL[right]  = dirichlet(Theo(radius(x,y,z) , 1., undercooling, F3(S)));

TS[top]    = dirichlet(TS_inf); 
TS[front]  = dirichlet(TS_inf); 
TS[back]   = dirichlet(TS_inf); 
TS[bottom] = dirichlet(TS_inf); 
TS[left]   = dirichlet(TS_inf); 
TS[right]  = dirichlet(TS_inf); 
/**
*/
/**
These are the parameters from the Almgren paper.
*/
double undercooling = -1/2., S = 1.56;
double t_fin = 0.4; // duration of the calculation, final time is therefore 2.
double t0 = 1.;
char filename [100];
FILE * fp1,  * fp2;

#include "../alex_functions2.h"

int main() {

  L0 = 10.;
  CFL = 0.5;
  origin (-0.5*L0, -0.5*L0, -0.5*L0);
  
  j = 1;
  for (j=1;j<=1;j++){
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");
    snprintf(filename, 100,  "out%d", j);
    fp2 = fopen (filename,"w");
/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 1;
    MAXLEVEL = 5+j, MINLEVEL = 4 ;
    TL.third = true;
    TS.third = true; 
    N = 1 << MAXLEVEL;
    init_grid (N);
    DT2 = 0.2*sq(L0/(1<< MAXLEVEL)); // dt = sq(dx)
    fprintf(stderr, "### DT2 %g \n", DT2);
    run(); 
    fclose(fp1);
    fclose(fp2);
  }
}

/**
*/


event init(t=0){
  dt = DT2;
  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;
  double R_init = 1.56;
  coord center = {0.,0.,0.};
    double F3_S  = F3(S);

  for (int k = 0; k < 2; k++){
  foreach(){
      dist[] = sphere(x,y,z,center,R_init);
  }
  boundary ({dist});
  restriction({dist});


  vertex scalar dist_n[];
  cell2node(dist,dist_n);

  fractions (dist_n, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  foreach(){
    double ray =  radius(x,y,z);
    if(cs[]>0.){
      TL[] = Theo(ray , 1., undercooling, F3_S);
    }
    else{
      TL[]=  0.;
    }
    TS[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
    LS_speed(
    dist,latent_heat,cs,fs,TS,TL,T_eq,
    vpc,vpcf,lambda1,lambda2,
    epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
    itrecons = 30,tolrecons = 1.e-12);
    event("adapt");
  } 


  foreach(){
    double ray =  radius(x,y,z);
    if(cs[]>0.){
      TL[] = Theo(ray , 1., undercooling, F3_S);
    }
    else{
      TL[]=  0.;
    }
    TS[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});


#if LS_curve
  curvature_LS(dist,curve);
#else
  curvature(cs,curve);
  boundary({curve});
#endif




  myprop(muv,fs,lambda[0]);

#if GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  u.r[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
  uf.r[embed] = dirichlet (0.);
#endif

  scalar myrad[];
  foreach(){
    if(cs[]*(1-cs[]) !=0.){
      coord n       = facet_normal( point, cs ,fs) , p;
      normalize(&n);
      double alpha  = plane_alpha (cs[], n);
      plane_area_center (n, alpha, &p);
      myrad[] = radius(x+Delta*p.x, y+Delta*p.y, z+Delta*p.z);
    }
    else
      myrad[] = nodata;
  }

  dump();
  exit(1);
}


event velocity(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = 40,tolrecons = 1.e-5
  );
  double DT3 = timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
}

/**
*/


event tracer_diffusion(i++){
  boundary({TL});
  myprop(muv,fs,lambda[0]);
  mgT = diffusion (TL,dt,muv,theta = cm);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[1]);
  boundary({TS});
  mgT = diffusion (TS,dt,muv,theta = cm);
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
  itredist = 20,
  tolredist = 3.e-3,
  itrecons = 40,
  tolrecons = 1.e-5,
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
event intefaces(i++,last){
  draw_vof("cs");
  save("cs.mp4");
}


event movies (i++,last;t<t_fin)
{
// the radius we output is just the maximum of the radius

  boundary({TL,TS});

  double y_max=0,y_max2 = 0.;
  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs) && y >0){
      double yy = y+Delta*height(h.y[]);
      if(yy < 1.e10)y_max = max(y_max,yy);
      if(x == Delta/2){
        y_max2 = y - dist[];
      }
    }     
  }
  // y_ref = sqrt(R^2 - x^2), x=Delta/2
  // double ref = sqrt(sq(S)*(t+t0+dt)-sq(L0/(2.*(1 << MAXLEVEL))));
  fprintf(fp1, "%g %g %g\n", sqrt(t), y_max, y_max2);
}

event final_error(t=end, last){
    scalar e[], ep[], ef[];

/**
Here we output the error given by the theory of [Frank](#frank1950radially).
*/
  double F3_S  = F3(S);
  foreach() {
    double ray =  radius(x,y,z);
    if (cs[] == 0.)
      e[] = nodata;
    else {
      e[] = TL[] - Theo(ray , t +t0+dt, undercooling, F3_S);
      ep[] = cs[] < 1. ? e[] : nodata;
      ef[] = cs[] >= 1. ? e[] : nodata;
    }
  }
  norm n = normf (e), np = normf (ep), nf = normf (ef);
  fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %g\n",
   N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, t);
  fflush (ferr);
  // if(MAXLEVEL == 5){
    scalar loge[];
    foreach(){
      if (cs[] == 0.)
        loge[] = nodata;
      else
        loge[] = log(fabs(e[]))/log(10.);
    }
    view(fov = 8.7, width = 1000, height = 1000);
    cells();
    squares("loge", min = -4, max = -2.5);
    draw_vof("cs");
    save("error.png");


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
    dump();
  // }
  fprintf(stderr, "### 1.e-2i2ale %d\n", i);
}

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

  boundary({cs,cs2,fs,fs2});
  fractions_cleanup(cs,fs,smin = 1.e-14);
  fractions_cleanup(cs2,fs2,smin = 1.e-14);
  restriction({cs,cs2,fs,fs2});
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});

  adapt_wavelet ({cs,visu,normvpc},
    (double[]){1.e-3,1.e-2,2.e-2},MAXLEVEL, MINLEVEL);


  double solid_volume = 0.;
  int n=0;
  foreach(reduction(+:n) reduction(+:solid_volume)){
    n++;
    solid_volume += (1.-cs[])*powf(Delta,3.);
  }
  fprintf(stderr, "##time %g nbcells %d solid volume %g\n", 
    t+1., n ,solid_volume);
#if LS_curve
  curvature_LS(dist,curve);
#else
  curvature(cs,curve);
  boundary({curve});
#endif
}
#endif

/**


~~~bib

@article{frank1950radially,
  title={Radially symmetric phase growth controlled by diffusion},
  author={Frank, Frederick Charles},
  journal={Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences},
  volume={201},
  number={1067},
  pages={586--599},
  year={1950},
  publisher={The Royal Society London}
}

@Article{Almgren1993,
  author        = {R. Almgren},
  title         = {Variational Algorithms and Pattern Formation in Dendritic Solidification},
  year          = {1993},
  volume        = {106},
  pages         = {337-354},
  issn          = {0021-9991},
  doi           = {10.1006/jcph.1993.1112},
}

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
