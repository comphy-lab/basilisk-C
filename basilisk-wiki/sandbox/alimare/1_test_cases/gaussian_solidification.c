/**
#Solidification of an interface slightly perturbed by a Gaussian function

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations can be found on the mini_cell test
case.

![Animation of cs*u.x + (1-cs)*u2.x.](gaussian_solidification/visu.mp4)(loop)

![Animation of cs*u.x + (1-cs)*u2.x.](gaussian_solidification/v_pcy.mp4)(loop)

We output only the interface at different times during the simulation.

~~~gnuplot Evolution curvature
set term pngcairo size 800,600 enhanced font 'Helvetica ,18'
set output 'curve_gaussian.png'
set xlabel 't'
set ylabel 'Curvature'
plot 'log' u 1:7 w l lw 2 t 'max', \
    'log' u 1:6 w l lw 2 t 'min'
~~~

~~~gnuplot Evolution of the interface (zoom)
set output 'v_pc_gaussian.png'
set xlabel 't'
set ylabel 'phase change speed'
plot 'log' u 1:8 w l lw 2 t 'max', \
    'log' u 1:9 w l lw 2 t 'min'

~~~

~~~gnuplot Evolution of the interface (zoom)
set term pngcairo size 900,600 enhanced font 'Helvetica ,36'
set output 'interfaces_gaussian.png'
set xtics offset 0,graph 0.05
set xtics -2,1,2
set ytics -1,0.5,0
unset xlabel
unset ylabel
set yrange [-1:0]
plot 'out' w l lw 4 t 'Interface'
~~~

As one can see, the Gaussian bump, slightly melts in the middle. This is due to
the Gibbs-Thomson formulation that we use :

$$
T_m =  T_{eq} - \epsilon_\kappa \kappa - \epsilon_v v_{pc}
$$
where $T_m$ is the local equilibrium temperature on the interface, $T_{eq}$ is
the phase change temperature (here 0), $\kappa$ is the local curvature of the
interface, $v_{pc}$ the phase change velocity and $\epsilon_\kappa, \epsilon_v$
are constants.

In the middle of the domain, the local temperature on the interface is such that
it is above the melting temperature $T_m$, therefore the ice melts.
*/

#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1

#include "embed.h"
#include "advection.h"
#include "diffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../level_set.h"

#define T_eq          0.
#define TL_inf        1.
#define TS_inf       -1.

int MIN_LEVEL, MAXLEVEL; 
double H0;
double latent_heat;
char filename [100];
FILE * fp1;

#define DT_MAX  1.

#define T_eq         0.


#define plane(x, y, H) (y - exp(-8.*(x)*(x))/L0 - H)
// #define plane(x, y, H) (y - H)


scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};

vector v_pc[];
scalar * LS_speed   = {v_pc.x,v_pc.y};
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];

#if Gibbs_Thomson
double  epsK = 0.005, epsV = 0.005;
#else
double  epsK = 0.000, epsV = 0.000;
#endif
scalar curve[];



double lambda[2];

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  
mgstats mg1,mg2;

TL[embed] = dirichlet(T_eq + epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq + epsK*curve[]-epsV*sqrt(v_pc.x[]*v_pc.x[]+v_pc.y[]*v_pc.y[]));
TS[bottom] = dirichlet(TS_inf); 

int j;
int k_loop = 0;
/**
The domain is 4 units long, centered vertically. */

double timestep_LS (const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}


int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  j = 1;
  for (j=1;j<=1;j++){

/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 100./j;
    MAXLEVEL = MIN_LEVEL = 7;

    H0 = -0.2*L0; 
    N = 1 << MAXLEVEL;
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");

    init_grid (N);
    run();
    // fclose(fp1); 
  }
}

event init(t=0){

  TOLERANCE = 2.e-7;
  DT = L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;

  foreach(){
      dist[] = clamp(plane(x,y,H0),-1.1*NB_width, 1.1*NB_width);
  }
  boundary ({dist});
  restriction({dist});

  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  curvature(cs,curve);
  boundary({curve});
  stats s2    = statsf(curve);

  foreach() {
    TL[] = TL_inf;
    TS[] = T_eq;
  }

  foreach_face(){
    v_pc.x[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});
}


event properties(i++){

  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  int kk;
  for (kk=1;kk<=(10*k_loop+1);kk++){
    if(i%2==0){
      mg1 = diffusion(TL, L0/(1 << MAXLEVEL), D = muv);
    }
    else{
      mg2 = diffusion(TS, L0/(1 << MAXLEVEL), D = muv );
    }
  }
}


event LS_advection(i++,last){
  if(i%2 == 1 && i> 100){
    double L_H       = latent_heat;  

    scalar cs0[];


    foreach(){
      cs0[]   = cs[];
      cs[]    = 1.-cs[];
    }
    foreach_face(){
      fs.x[]  = 1.-fs.x[];
    }

    boundary({cs,fs,cs0});
    restriction({cs,fs});

    phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width,
      nb_cell_NB,lambda,epsK, epsV);

    recons_speed(dist, 0.5*DT, nb_cell_NB, NB_width, LS_speed);

    dt = timestep_LS (v_pc, DT);
    
    boundary((scalar *){v_pc});

    advection_LS (level_set, v_pc, dt);
    boundary ({dist});
    restriction({dist});

    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

    curvature(cs,curve);
    boundary({curve});
    

    foreach(){
      cs[]      = 1.-cs[];
    }
    foreach_face(){
      fs.x[]      = 1.-fs.x[];
    }

    boundary({cs,fs});
    restriction({cs,fs});

    k_loop = 0;
    foreach(){
      if(cs0[] != 1. && cs[] ==1.)k_loop = 1;
    }
  }
}


event LS_reinitialization(i++,last){
  if(i>0 && i%2==1){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      1);
  }
}

#if DOUBLE_EMBED
event double_calculation(i++,last){
// modify the fs , cs, copy the outer fields in the partially covered cells

  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face(){
    fs.x[]      = 1.-fs.x[];
  }

  boundary({cs,fs});
  restriction({cs,fs});
}
#endif

event movies ( i++,last;t<150.)
{
  stats s2    = statsf(curve);
  stats s3    = statsf(v_pc.y);
  if(i%40 == 1) {

    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    
    view (fov = 20.5678, quat = {0,0,0,1}, tx = 0.0308985, 
      ty = 0.0325629, bg = {1,1,1}, width = 800, height = 800, samples = 1);    
    // draw_vof("cs");
    // squares("curve", min =-5., max = 5.);
    // save ("curve.mp4");

    draw_vof("cs");
    squares("visu", min =-1., max = 1.);
    save ("visu.mp4");

    // boundary((scalar *){v_pc});
    // draw_vof("cs");
    // squares("v_pc.y", min =s3.min, max = s3.max);
    // save ("v_pcy.mp4");



    fprintf(stderr, "## %g %g %d\n", mg1.resa, mg2.resa, k_loop);
    
    if(i%600==1) {
      output_facets (cs, stdout);
    }
  }
  if(i%2 ==1){
    Point p     = locate(-1.51*L0/(1<<MIN_LEVEL),-3.51*L0/(1<<MIN_LEVEL));

    double cap  = capteur(p, TL);
    double cap3 = capteur(p, TS);
    double cap2 = capteur(p, cs);
    double T_point = cap2*cap + (1.-cap2)*cap3;
   fprintf (stderr, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
    t, cap, cap2, cap3, T_point, s2.min, s2.max, s3.min, s3.max);
   }
 
}