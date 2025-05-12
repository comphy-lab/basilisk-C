/**
#Mullins Sekerka instability

This theory of this test case has been studied originally by [Mullins and
Sekerka](#Mullins1964). Then it was used as a validation test case for several
numerical method (see the work of [Almgren](#Almgren1993) and [Chen et al.](#chen_simple_1997))

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_2 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations can be found on the mini_cell test
case.

Here we plot the same figure Fig. 4 of [Chen et al.](#chen_simple_1997)

~~~gnuplot Evolution of the interface
unset xlabel
unset ylabel
set term svg size 1200,400
plot 'out0' w l lc 'blue'   t '32x32',\
     'out1' w l lc 'red'    t '64x64',\
     'out2' w l lc 'black'  t '128x128'
~~~


~~~gnuplot Interface position
set term svg size 900,900
set key left

ftitle(a,b,c) = sprintf("%.3fexp^{%4.2ft}-%.2f", exp(a), b,c)
f(x) = a*exp(b*x)-c
f2(x) = a2*exp(b2*x)-c2
f3(x) = a3*exp(b3*x)-c3
fit f(x) 'log0' u 1:2 via a,b,c
fit f2(x) 'log1' u 1:2 via a2,b2,c2
fit f3(x) 'log2' u 1:2 via a3,b3,c3

set logscale xy
plot 'log0' u 1:2  w p pt 7 lc 'blue'  t '32x32', \
  'log1' u 1:2  w p pt 7 lc 'red'  t '64x64', \
  'log2' u 1:2  w p pt 7 lc 'black'  t '128x128', \
  f(x) t ftitle(a,b,c) lc 'blue',\
  f2(x) t ftitle(a2,b2,c2) lc 'red',\
  f3(x) t ftitle(a3,b3,c3) lc 'black'
~~~


*/

#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 0
#define Pi 3.14159265358979323846
#define QUADRATIC 1
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
#define TS_inf        0.

#define tstart 0.

int MINLEVEL, MAXLEVEL; 
double latent_heat;

#define DT_MAX  1.

#define T_eq         0.


#define plane(x, y, A, n, L0) (y - A*cos(2.*n*Pi*x*L0))

double T0(double x,double y, double V, double t){
  return -1+exp(-V*(y-V*t));
}

double TL_inf(double x, double y, 
  double V, double t, double A, int n, double L0){

  double ystar = y - A*cos(2.*n*Pi*x*L0);
  if(ystar>V*t)
    return -1+exp(-V*(ystar-V*t));
  else
    return 0.;
}


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



double  epsK, epsV = 0.;


double lambda[2];
int aniso = 1;


int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

double s_clean = 1.e-10; // used for fraction cleaning

  
mgstats mg1,mg2;

scalar curve[];
TL[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, aniso));
TS[embed] = dirichlet(T_eq + Temp_GT(point, epsK, epsV, vpc, curve, fs, cs, aniso));



#define V  1.
#define t0 0.
double A; 
double DomainSize;
int n;


TL[top]    = dirichlet(TL_inf(x,y,V,t-t0,A,2,DomainSize)); 
TS[bottom] = dirichlet(TS_inf); 


int j;
int k_loop = 0;
double T_final = 0.1;
int nb_img = 30;

/**
Variable for varying parameters
*/
char filename [100];
FILE * fp1,  * fp2;


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
  periodic(left);

  DomainSize = L0 = 1.;
  CFL = 0.5;
  origin (-0.5*L0, -0.5*L0);
  // dist[top]    = neumann(-1.);  
  // dist[bottom] = neumann(-1.); 
  j = 1;
  for (j=0;j<=2;j++){
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");
    snprintf(filename, 100,  "out%d", j);
    fp2 = fopen (filename,"w");
    // A = 0.005*(j+1);
    A = 0.01;
#if Gibbs_Thomson
    epsK = 1.e-3;
#else
    epsK = 0.;
#endif
/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 1.;
    MAXLEVEL = 5+j, MINLEVEL = 5 ;
    DT2  = sq(L0/(1<< MAXLEVEL));
    N = 1 << MAXLEVEL;
    TL.third = true;
    TS.third = true;  

    init_grid (N);
    run(); 
    fclose(fp1);
    fclose(fp2);
  }
}

event init(t=0){
  NB_width = nb_cell_NB*L0/(1<<MAXLEVEL);
  n= 2;
  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach(){
      dist[] = clamp(plane(x,y,A,n,DomainSize),-1.2*NB_width, 1.2*NB_width);
  }
  boundary ({dist});
  restriction({dist});

  vertex scalar distn[];
  cell2node(dist,distn);

  fractions (distn, cs, fs);
  fractions_cleanup(cs,fs,smin = s_clean);
  boundary({cs,fs});
  restriction({cs,fs});

  foreach(){
    if(cs[]>0.){
      TL[] = TL_inf(x,y,V,t,A,n,DomainSize);
    }
    else{
      TL[]=  0.;
    }
    TS[] = 0.;
    boundary({TS,TL});
  }

  foreach_face(){
    vpc.x[] = 0.;
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



event velocity(i++){

  double lambda1 = lambda[0], lambda2 = lambda[1]; 
  LS_speed(
  dist,latent_heat,cs,fs,TS,TL,T_eq,
  vpc,vpcf,lambda1,lambda2,
  epsK,epsV,aniso,deltat=0.45*L0/(1<<MAXLEVEL),
  itrecons = 60,tolrecons = 1.e-12
  );
  double DT3 = timestep_LS(vpcf,DT2,dist,NB_width);
  tnext = t+DT3;
  dt = DT3;
  fprintf(stderr, "## %g %g %g\n", t, DT2, DT3);
}

event tracer_diffusion(i++){
  boundary({TL});
  myprop(muv,fs,lambda[0]); // MANDATORY, the interface has moved !!!!!
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
  epsK,epsV,aniso,
  curve,
  &k_loop,
  deltat = 0.45*L0 / (1 << grid->maxdepth),
  itredist = 10,
  tolredist = 3.e-3,
  itrecons = 60,
  tolrecons = 1.e-12,
  s_clean = 1.e-10,
  NB_width);

  foreach_face(){
    uf.x[] = 0.;
  }
  boundary((scalar *){uf});
  restriction((scalar *){uf});  
}

event interface_output(i++,last){
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
  fprintf(fp1, "%g %g %g\n", t+dt+t0, y_max-t-dt-A,y_LS-t-dt-A);
}

event interface2(t+=0.01,last){
  output_facets (cs, fp2);
}

event movies ( t+=0.005,last;t < T_final)
{
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  squares("visu", min = -1 , max = 0.);
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
    (double[]){1.e-3,2.e-4},MAXLEVEL, MINLEVEL);
  foreach(reduction(+:n)){
    n++;
  }
}
#endif



/**


~~~bib

@Article{Mullins1964,
  author        = {Mullins, William W and Sekerka, RF},
  title         = {Stability of a planar interface during solidification of a dilute binary alloy},
  journal       = {Journal of applied physics},
  year          = {1964},
  volume        = {35},
  number        = {2},
  pages         = {444--451},
  publisher     = {AIP},
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
