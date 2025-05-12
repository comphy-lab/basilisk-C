/**
we want to simulate the rise of a lower density fluid bubble at the surface of water with air over.    
rho_water = 1000.  
rho_bubble = 700.  
rho_air = 1.21   
gravity 9.81   

[rising lens](http://bipbip.ida.upmc.fr/~ray/three_phases/mov_rising_lens.mpg)
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "navier-stokes/conserving3.h"
#include "tension.h"

#define radius 0.2
#define Dre 2.7e-3
#define Vre 0.257
#define Mu_re2 20e-3
#define Mu_re3 20e-3
#define Rho2_re 1000.
#define Rho3_re 700. 
#define Rho1_re 1.21
#define Re (Rho2_re*Vre*Dre/Mu_re2)
#define SIGMA21 32e-3
#define SIGMA31 40e-3
#define SIGMA32 56e-3
#define gad 1.

#define vrho2 1.
#define vrho3 (Rho3_re/Rho2_re)
#define vrho1 (Rho1_re/Rho2_re)
#define vmu2 (2.*radius/Re)
#define vmu3 (2.*radius/Re*(Mu_re3/Mu_re2))
#define vmu1 (2.*radius/Re*(1.81e-5/Mu_re2))

#define ycc 0.
#define hpool (5.*radius*2.)
#define xcc0 0.
#define xcc2 (xcc0 + hpool)
#define hgout (0.25*radius*2.)
#define xcc (xcc2 - 1.5*radius)

#define circle(x,y) (sq(x - xcc) + sq(y - ycc))
#define drop(x,y) (sq(radius) - circle(x,y))
#define pool(x) (xcc2-x)
#define pool2(x,y) (min(xcc2-x,-drop(x,y)))

#define epi radius/8.
#define rfdrop(x,y) ((sq(radius-epi) < circle(x,y)) && (sq(radius+epi) > circle(x,y)))
#define rfpool(x) (((xcc2-epi) < x) && ((xcc2+epi) > x))

int maxlevel = 7;
int minlevel = 1;
double uem = 0.01;
double fem = 0.01;

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.t[right] = neumann(0.);  
p[right]   = neumann(0.);  
u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0.);
p[bottom]   = neumann(0.);
u.n[top] = dirichlet(0.);
u.t[top] = neumann(0.); 
p[top]   = neumann(0.);

int main (int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uem = atof (argv[2]);

  init_grid (64);
  origin (0., 0.);
  size (4.);

  rho1 = vrho1; 
  rho2 = vrho2;
  rho3 = vrho3;
  mu1 = vmu1;
  mu2 = vmu2; 
  mu3 = vmu3;
  f2.sigma = (SIGMA21 + SIGMA32 - SIGMA31)/2.;
  f3.sigma = (SIGMA31 + SIGMA32 - SIGMA21)/2.;
  f1.sigma = (SIGMA31 + SIGMA21 - SIGMA32)/2.;
  run();
}

event init (t = 0) {

  scalar f4[], f5[], f6[];

  refine ((rfpool(x) || rfdrop(x,y)) && (level < maxlevel));

  fraction (f4, pool2(x,y));
  fraction (f6, drop(x,y));
  foreach() {

    f2[] = f4[];
    f3[] = f6[];
    f1[] = clamp(1.-f2[]-f3[],0.,1.); 
    rhov[] = rho(f2[],f3[]);
  }
  boundary ((scalar *){f1,f2,f3,rhov});
}
/** 
ouput gfs files */   
/**
event snapshot (t = 0.; t += 0.020; t <= 10.000) {   
  char name[80];   
  sprintf (name, "lensaxi-%5.3f.gfs", t);   
  scalar vel[];   
  foreach() {    
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));   
  }   
  boundary ((scalar *){vel});    
  output_gfs (file = name, t = t, list = {u.x,u.y,vel,p,f2,f3,f1});   
}*/    

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= gad;
}

event adapt (i++) {
  scalar ff2[],ff3[],r3[];
  foreach() {
    ff2[] =  (2.*f2[] - 1.);
    ff3[] =  (2.*f3[] - 1.);
    r3[]= (2.*f3[]*f2[]*f1[] - 1.);
  }
  boundary ({ff2,ff3,r3});    
  adapt_wavelet ({r3,ff2,ff3,u,f2,f3}, (double[]){0.01,0.1,0.1,uem,uem,fem,fem}, maxlevel, minlevel);
  event ("properties");
}