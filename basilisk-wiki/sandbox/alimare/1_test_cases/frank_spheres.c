/**
#Frank's Spheres

This theory of this test case has been studied originally by 
[Frank](#frank1950radially). We start with a solid circle (or sphere in 3D) at melting
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
fit f2(x) 'log' u (log($1)):(log($4)) via a2,b2

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
     '' u 1:4 pt 7 lc 'red' t 'partial cells', exp(f2(log(x))) lt 2 lw 1.5 t ftitle(a2,b2) , \
     '' u 1:6 pt 7 lc 'blue' t 'full cells'  , exp(f3(log(x))) lt 4 lw 1.5 t ftitle(a3,b3)
~~~


~~~gnuplot Evolution of the interface
unset logscale
unset format y
set xrange [0:2.5]
set yrange [0:2.5]
set xtics 0,0.5,2.5
unset ylabel
plot 'out0' w l lc 'black' lw 1.5 t '32x32',\
     'out1' w l lc 'red'   lw 1.5 t '64x64',\
     'out2' w l lc 'blue'  lw 1.5 t '128x128'
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


int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
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

#include <gsl/gsl_sf_expint.h>
#pragma autolink -lgsl -lgslcblas

#if dimension == 2
double F2(double s) {return gsl_sf_expint_E1(pow(s,2.)/4.);}
#endif

double Theo(double r, double t, double undercooling, double F2_S){
  if(dimension == 2){;
    double s = r/pow(t,1./2);
    return undercooling * (1 - F2(s) / F2_S);
  }
  exit(1); // dimension == 3 not coded
}


TL[embed] = dirichlet(T_eq);
TS[embed] = dirichlet(T_eq);

#define radius(x, y) (sqrt(sq(x)+sq(y)))
TL[top]    = dirichlet(Theo(radius(x,y) , 1., undercooling, F2(S)));
TL[bottom] = dirichlet(Theo(radius(x,y) , 1., undercooling, F2(S)));
TL[left]   = dirichlet(Theo(radius(x,y) , 1., undercooling, F2(S)));
TL[right]  = dirichlet(Theo(radius(x,y) , 1., undercooling, F2(S)));

TS[top]    = dirichlet(TS_inf); 
TS[bottom] = dirichlet(TS_inf); 
TS[left]   = dirichlet(TS_inf); 
TS[right]  = dirichlet(TS_inf); 
/**
*/



/**
These are the parameters from the Almgren paper.
*/
double undercooling = -1/2., S = 1.56;
double t_fin = 1.; // duration of the calculation, final time is therefore 2.
double t0 = 1.;
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

  L0 = 16.;
  CFL = 0.5;
  origin (-0.5*L0, -0.5*L0);
  
  j = 1;
  for (j=0;j<=2;j++){
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
  coord center = {0.00,0.00};
  foreach(){
      dist[] = circle(x,y,center,R_init);
  }
  boundary ({dist});
  restriction({dist});


  vertex scalar dist_n[];
  cell2node(dist,dist_n);

  fractions (dist_n, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

#if LS_curve
  curvature_LS(dist,curve);
#else
  curvature(cs,curve);
  boundary({curve});
#endif


  double F2_S  = F2(S);
  scalar error[],theo[];
  foreach(){
    double ray =  radius(x,y);
    if(cs[]>0.){
      TL[] = Theo(ray , 1., undercooling, F2_S);
    }
    else{
      TL[]=  0.;
    }
    TS[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});
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


event velocity(i++){
  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,eps4,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = 40,tolrecons = 1.e-12
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
  itrecons = 40,
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
event intefaces(t+=0.1,last){
  output_facets(cs,fp2);
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
  double ref = sqrt(sq(S)*(t+t0+dt)-sq(L0/(2.*(1 << MAXLEVEL))));
  fprintf(fp1, "%g %g %g %g %g\n", t, y_max, S*pow(t+t0+dt,1./2),y_max2,ref);
}



event final_error(t=end, last){
    scalar e[], ep[], ef[];

/**
Here we output the error given by the theory of [Frank](#frank1950radially).
*/
  double F2_S  = F2(S);
  foreach() {
    double ray =  radius(x,y);
    if (cs[] == 0.)
      e[] = nodata;
    else {
      e[] = TL[] - Theo(ray , t +t0+dt, undercooling, F2_S);
      ep[] = cs[] < 1. ? e[] : nodata;
      ef[] = cs[] >= 1. ? e[] : nodata;
    }
  }
  norm n = normf (e), np = normf (ep), nf = normf (ef);
  fprintf (stderr, "%d %.3g %.3g %.3g %.3g %.3g %.3g %g\n",
   N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max, t);
  fflush (ferr);
  if(MAXLEVEL == 7){
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
    // dump();
  }
  fprintf(stderr, "### ite finale %d\n", i);
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
  int n=0;
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});

  squares("normvpc");
  draw_vof("cs");
  save("speed.mp4");

  adapt_wavelet ({cs,visu,normvpc},
    (double[]){1.e-3,1.e-3,1.e-2},MAXLEVEL, MINLEVEL);
  foreach(reduction(+:n)){
    n++;
  }
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
