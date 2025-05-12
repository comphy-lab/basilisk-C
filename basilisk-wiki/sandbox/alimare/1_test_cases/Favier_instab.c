/**
# Rayleigh Benard instability with a Stefan problem.

This test case is similar to the [Favier et al.](#Favier2019) case. 
Still some work to do regarding linear stability analysis.

![Temperature field](Favier_instab/temperature.mp4)
*/

#define BICUBIC 1
#define DOUBLE_EMBED  1
#define LevelSet      1
#define Gibbs_Thomson 0
#define Pi 3.14159265358979323846


#include "embed.h"
#include "../double_embed-tree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"

#include "../level_set.h"
#include "../LS_advection.h"

#define T_eq         0.05
#define TL_inf       1.
#define TS_inf       0.

#define Ra 256
#define St 10
#define epsilon 4.e-3
#define alpha 2500
#define h0 0

#define MAXLEVEL 7
#define MINLEVEL 3


scalar TL[];
scalar * tracers = {TL};

scalar TS[];
scalar * tracers2 = {TS};

scalar dist[];
scalar * level_set = {dist};

vector vpc[];
scalar * LS_speed  = {vpc.x,vpc.y};

face vector muv[];

stats s3;

double  latent_heat = 10.;
double  lambda[2];   

// no Gibbs Thomson effect.
double  epsK = 0.000, epsV = 0.000;
int aniso = 1;

int k_loop = 0;

scalar curve[];

int     nb_cell_NB;
double  NB_width ;    // length of the NB

mgstats mgd;  
TL[bottom]    = dirichlet(TL_inf);
TL[embed]     = dirichlet(T_eq);

TS[top]       = dirichlet(TS_inf);
TS[embed]     = dirichlet(T_eq);

u.n[embed]    = dirichlet(0.);
u.t[embed]    = dirichlet(0.);

double nu = 0.001;
face vector muc[], av[];

double TLfield (double x, double y, double Thetam, double H0){
  y+=0.5;
  H0+=0.5;
  if (y>H0)
    return Thetam;
  else
    return 1+(Thetam-1)*y/H0;
}

double TSfield (double x, double y, double Thetam, double H0){
  y+=0.5;
  H0+=0.5;
  if (y>H0)
    return Thetam*(y-1)/(H0-1);
  else
    return Thetam;
}

int main() {
  L0 = 1;
  origin (-L0/2,-L0/2);
  int N = 1 << MAXLEVEL;
  init_grid (N);
  mu = muc;
  a = av;
  run();
}

event init (t = 0) {
  DT = 0.8*(L0)/( 1 << MAXLEVEL);
  periodic(right);
  nb_cell_NB = 1 << 3 ;  //
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);
  foreach()
    dist[] = -y;
  boundary ({dist});
  restriction({dist});
  LS_reinit(dist);
  vertex scalar dist_n[];
  cell2node(dist,dist_n);
  fractions (dist_n, cs, fs);
  fractions_cleanup(cs,fs);
  boundary({cs,fs});
  restriction({cs,fs});
  foreach_face(){
    vpc.x[] = 0.;
  }
  boundary((scalar *){vpc,av});
  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach() {
    TL[] = TLfield(x,y,T_eq,0.);
    TS[] = TSfield(x,y,T_eq,0.);
  }
  boundary({TL,TS});
  restriction({TL,TS});
}

event stability (i++) {
/**
Sometimes, because the embedded boundary is moving, some cells have uf.x[] != 0.
&& fm.x[] == 0. Therefore, we do an ugly quickfix for cell that are incoherent
with the embedded boundary stability conditions. ()
*/
  foreach_face()
    if(fm.x[]==0. && uf.x[] !=0.) uf.x[] = 0.;
  boundary((scalar *){uf});
}


event properties (i++){
  foreach_face()
    muc.x[] = fm.x[]*nu;
  boundary((scalar*){muc});
}

event tracer_diffusion (i++) {
  mgd = diffusion (TL,dt,muc,theta = cm);
  invertcs(cs,fs);
  event("properties");
  mgd = diffusion (TS,dt,muc,theta = cm);
  invertcs(cs,fs);
  event("properties");

}

event acceleration (i++) {
  foreach_face(x)
    av.x[] = 0;
  foreach_face(y)
    av.y[] = cs[]*1000*(TL[]+TL[-1])/2.;
  boundary((scalar *) {av});

}

event LS_advection(i++,last){  
  double lambda1 = lambda[0], lambda2 = lambda[1];
  
  advection_LS(
  dist,
  latent_heat,
  cs,fs,
  TS,TL,
  T_eq,
  vpc,
  lambda1,lambda2,
  epsK,epsV,aniso,
  curve,
  &k_loop,
  k_limit = 0,
  overshoot = 1.,
  deltat = 0.8*L0 / (1 << grid->maxdepth),
  itredist = 5,
  tolredist = 2.e-2,
  itrecons = 30,
  tolrecons = 1.e-5,
  s_clean = 1.e-10,
  NB_width);

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
  fractions_cleanup(cs,fs,smin = 1.e-10);
  fractions_cleanup(cs2,fs2,smin = 1.e-10);
  restriction({cs,cs2,fs,fs2});
  int n=0;
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  adapt_wavelet ({cs,visu,u},
    (double[]){0.005,0.001,0.005,0.005},MAXLEVEL, MINLEVEL); // cs, visu, u.x,
  // u.y
  foreach(reduction(+:n)){
    n++;
  }
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, T_eq, vpc, latent_heat,
    lambda, epsK, epsV, aniso);
  fprintf(stderr, "##Time %g Nbcells %d\n",t, n);
}
#endif

event movie1 (t+=0.005,last){
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  draw_vof("cs");
  squares("visu", min = 0. , max = 1.);
  save("temperature.mp4");

  
  // draw_vof ("cs", "fs", fc = {0.9, 0.9, 0.9});
  // squares ("TL", min = 0, max = 1);
  // save ("TL.mp4");

  Point p = locate(0.,-0.1);

  double val  = capteur(p, TL);
  double val2 = capteur(p, u.y);
  fprintf(stderr, "##ite %d temps %g valeur temperature %g vitesse %g \n", i,t,  val,val2);
}

// event movie2 (i = 2; i++){
//   scalar Theta[];
//   foreach()
//     Theta[] = u.y[];
//   boundary({Theta});
//   restriction({Theta});
//   output_ppm(Theta,linear = true,spread = 2, file ="uy.mp4",n=200);
// }

event final (t = 0.8)
{
  // char name[80];
  // sprintf (name, "mu-%g.png", t);
  // output_ppm (TL, file = name, n = 200, linear = true, spread = 2);
  draw_vof ("cs", "fs", fc = {0.9, 0.9, 0.9});
  squares ("TL", min = 0, max = 1);
  cells();
  save ("TL.png");
  dump();
}

/**

~~~bib
@Article{Favier2019,
  author        = {B. Favier and J. Purseed and L. Duchemin},
  title         = {Rayleigh–Bénard convection with a melting boundary},
  year          = {2019},
  volume        = {858},
  pages         = {437-473},
  issn          = {0022-1120},
  doi           = {10.1017/jfm.2018.773},
}

~~~

*/