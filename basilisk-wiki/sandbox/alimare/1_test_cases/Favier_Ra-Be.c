/**
# Rayleigh-Benard convection with a melting boundary

Test case taken from [Favier et al., 2019](#Favier2019). Mimics condition where
a liquid is heated from below and its solid is above it and melts. As the upper
boundary (interface solid/liquid) melts the apparent Rayleigh number changes and
the flow characteristics change.

A Navier-Stokes solver is used in the liquid and the heat equation is solved in
the solid.


Global Rayleigh is $10^5$.

![Animation of the temperature](Favier_Ra-Be/temperature.mp4)
*/


/**
#Figures


~~~gnuplot Average Height
set key bottom
set xlabel "Time"
set ylabel "Average height"
plot 'out' u 1:4 w l t 'Height'
~~~

~~~gnuplot Kinetic energy versus effective Rayleigh number
set key right
set xlabel "Ra_{e}"
set xrange[0:2000]
set xtics('1000' 1000,'1707.76' 1707.76,'2000' 2000)
set logscale y
set yrange[0.0000000000001:1.]
set ylabel 'Kinetic energy density'
plot 'out' u 2:3 w l t 'Energy density'
~~~ 
*/

/**
Required for using hybrid level-set/embedded boundary method
*/
#define DOUBLE_EMBED 1
#define Gibbs_Thomson 0
#define Pi 3.14159265358979323846
#define QUADRATIC 1
#define GHIGO 1
#define NS_emerged 1

#include "ghigo/src/myembed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#include "curvature.h"
#include "view.h"
// #define dtLS 1
#include "alimare/LS_diffusion.h"

/**
## Miscalleneous functions for this test case

Function to initalize TL
Convention:
-cell is completely fluid => frac == 1
-cell is completely solid => frac == 0
*/

#define ratio 8.            // Ratio of length to width in the domain
double TLfield (double y, double frac, double Thetam, double H0){
  y += ratio/2.; // change of reference for having y = 0 in the bottom rather
  // than y = -4
  if (frac > 0.)
    return 1.+(Thetam-1.)*y/H0;
  else
    return Thetam;
}

/**
Function to initialize TS
*/

double TSfield (double y, double frac, double Thetam, double H0){
  y += ratio/2.;
  if (frac < 1.)
    return Thetam*(y-1.)/(H0-1.);
  else
    return Thetam;
}

/**
Function to calculate the kinetic energy density
*/

double energy_density(){
  double se = 0.;
  double ke = 0.;
  foreach(reduction(+:ke))
    if (cs[]>1.e-6){
      ke+=dv()*(sq(u.x[])+sq(u.y[]));
    }
  foreach(reduction(+:se)){
    if (cs[]>1.e-6)
      se+=dv();
  }
  if(se==0.){
    fprintf(stderr, "#ERROR no liquid cells\n");
    exit(1);
  }
  return ke/se;
}

/**
Function to calculate the average height of the interface
*/

double output_height (struct OutputFacets p){
  scalar c = p.c;
  face vector s = p.s;
  if (!p.fp) p.fp = stdout;
  if (!s.x.i) s.x.i = -1;
  double hsum = 0.;
  double xsum = 0.;
  foreach(reduction(+:hsum))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha1 = plane_alpha (c[], n);
      coord segment[2];
      if (facets (n, alpha1, segment) == 2){
        hsum += y+(segment[0].y*Delta+segment[1].y*Delta)/2.;
      }
    }
  foreach(reduction(+:xsum))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha1 = plane_alpha (c[], n);
      coord segment[2];
      if (facets (n, alpha1, segment) == 2){
        xsum += 1.;
      }
    }
  return hsum/xsum+ratio/2.;
  fflush (p.fp);
}


/**
## Parameters for the simulation.
*/


#define TL_inf       1.     // T1 in the paper
#define TS_inf       0.     // T0 in the paper

double Ra;             // Global Rayleigh number
double h0;             // Initial height of the interface

int MAXLEVEL;
#define MINLEVEL 3



mgstats mgd;  

TL[bottom]  = dirichlet(TL_inf);
TL[embed]   = dirichlet(T_eq);

TS[top]     = dirichlet(TS_inf);
TS[embed]   = dirichlet(T_eq);

double nu = 1.;             // Viscocity
face vector muc[], av[];


int main(int argc, char * argv[]) {

  if (argc > 1)
    Ra = atoi (argv[1]);
  else
    Ra = 1.e5; // default value

  if (argc > 2)
    MAXLEVEL = atoi (argv[2]);
  else
    MAXLEVEL = 8; // default value
  if(argc > 3)
    h0 = atof(argv[3]);
  else
    h0 = 0.05; // default value

  L0 = ratio;  
  origin (-L0/2,-L0/2);
  int N = 1 << MAXLEVEL;
  init_grid (N);
/**
We use mask() to create rectangle domain
*/
  mask(y > 1.-ratio/2. ? top : none); // 
  mu = muc;
  a = av;
  run();
}

event init (t = 0) {
  DT_LS = 0.45*(L0)/( 1 << MAXLEVEL);
  DT = 1.;
  CFL = 0.1;

  T_eq   =      0.3;    // Theta m in the paper
  latent_heat = 10;
  epsK = 0.;
  epsV = 0.;
  eps4 = 0.;
  itrecons = 50;
  tolrecons = 1.e-12;

  periodic(right);
  nb_cell_NB = 6 ; // maybe 8 ... 6 is sufficient 
  NB_width   = nb_cell_NB * L0 / (1 << MAXLEVEL);


/**
### Initialization of the interface and temperature field
*/
  foreach()
    dist[] = h0-ratio/2.-y;

  fprintf(stderr, "## init %g %g %g %g\n",Ra, ratio/2.,h0, L0/(1<<MAXLEVEL));

  boundary ({dist});
  restriction({dist});
  LS_reinit(dist);

  LS2fractions(dist,cs,fs);

  foreach_face(){
    vpc.x[] = 0.;
  }
  boundary((scalar *){vpc,av});
  lambda[0] = 1.;
  lambda[1] = 1.;
  foreach() {
    TL[] = TLfield(y,cs[],T_eq,h0);
    TS[] = TSfield(y,cs[],T_eq,h0);
  }
  boundary({TL,TS});
  restriction({TL,TS});

#if GHIGO // mandatory for GHIGO
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  // Added to account for moving boundary algorithm
  
  uf.n[embed] = dirichlet (0.);
  uf.t[embed] = dirichlet (0.);
#endif
}

event stability (i++) {
/**
Sometimes, because the embedded boundary is moving, some cells have uf.x[] != 0.
&& fm.x[] == 0. Therefore, we do an ugly quickfix for cell that are incoherent
with the embedded boundary stability conditions.
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

/**
The top condition works not well if i>1. So we do an ugly quickfix for
reforcement of mask and top condition.
**/

event renforcemask(i++){
  mask(y > 1.-ratio/2. ? top : none);
  foreach()
    if(y > 1.-ratio/2.)
      TS[] = TS_inf;
}

/**
Boussinesq term.
*/
event acceleration (i++) {
  foreach_face(x)
    av.x[] = 0;
  foreach_face(y)
    av.y[] = cs[]*Ra*(TL[]+TL[-1])/2.;
  boundary((scalar *) {av});

}

/**
FIXME :
Still some issues with mesh adaptation, especially an assert() in
refine_embed_linear().
*/
#if 1 
event adapt (i++, last) {

  foreach_cell(){
    cs2[] = 1.-cs[];
  }
  foreach_face(){
      fs2.x[] = 1.-fs.x[];
  }

  boundary({cs,cs2,fs,fs2});
  fractions_cleanup(cs,fs,smin = 1.e-10);
  fractions_cleanup(cs2,fs2,smin = 1.e-10);
  restriction({cs,cs2,fs,fs2});
  boundary({TL,TS});
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
  }
  boundary({visu});
  restriction({visu});
  adapt_wavelet ({cs,visu,u},
    (double[]){0.005,1.e-4,0.005,0.005},MAXLEVEL, MINLEVEL);
}
#endif


/**
## Outputs
*/
event movie1 (t+=0.005,last;t <9.){
  fprintf(stderr, "t = %g\n",t );

  view (fov = 2.75, ty = 0.44, width = 1600, height = 200);
  scalar visu[];
  foreach(){
    visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ; // linear combination of both
    // temperatures for visualization purposes
  }
  boundary({visu});
  restriction({visu});
  draw_vof("cs");
  squares("visu", min = 0. , max = 1.);
  save("temperature.mp4");  
}
/**
We create an event to output effective Rayleigh number, everge height of interface and kinetic energy density 
**/


event Calculation(i++){
  double Rae = Ra*(1-T_eq)*pow(output_height(cs,stdout),3);
  fprintf(stdout, "%g %g %g %g \n", t, Rae, energy_density(), output_height(cs,stdout));
  fprintf(stderr, "t = %g \n", t);
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
  __markedentry = {[limare:6]},
  doi           = {10.1017/jfm.2018.773},
}
~~~
*/