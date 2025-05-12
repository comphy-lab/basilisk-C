/**     
We try to simulate an oil lens on water */   

#include "axi.h"
#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

/**
for beginning we choose two fluids with same viscosity and density with air over
*/

#define radius 0.2
#define Dre 2.7e-3
#define Vre 0.257
#define Mu_re2 20e-3
#define Mu_re3 20e-3
#define Rho2_re 1000.
#define Rho3_re 1000. 
#define Rho1_re 1.21
#define Re (Rho2_re*Vre*Dre/Mu_re2)

/**
and choose surface tension (case two)
*/

#define SIGMA2 64e-3
#define SIGMA3 80e-3
#define SIGMA32 112e-3

/**
We adimension with gravity even if we don't use it, we will use it
after below.*/

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
#define xcc (xcc2)

#define circle(x,y) (sq(x - xcc) + sq(y - ycc))
#define drop(x,y) (sq(radius) - circle(x,y))
#define pool(x) (xcc2-x)
#define pool2(x,y) (min(xcc2-x,-drop(x,y)))

#define epi radius/8.
#define rfdrop(x,y) ((sq(radius-epi) < circle(x,y)) && (sq(radius+epi) > circle(x,y)))
#define rfpool(x) (((xcc2-epi) < x) && ((xcc2+epi) > x))

int maxlevel = 8;
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
  f2.sigma = (SIGMA2 + SIGMA32 - SIGMA3)/2.;
  f3.sigma = (SIGMA3 + SIGMA32 - SIGMA2)/2.;
  f1.sigma = (SIGMA3 + SIGMA2 - SIGMA32)/2.;
  run();
}

event init (t = 0) {

  scalar f4[], f6[];

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

event snapshot (t = 0.; t += 1.0; t <= 5.0) {
  char name[80];
  sprintf (name, "lensaxi-%3.1f", t);
  dump(name);
}

event movie (t =0.; t += 0.020; t <= 5.000) {

  clear();
  view (fov = 6.989, quat = {0,0,-0.707,0.707}, tx = 1e-6, ty = -0.5, width = 1200, height = 600);
  draw_vof ("f1", lw = 2);
  draw_vof ("f2", lw = 2);
  draw_vof ("f3", lw = 2);
  mirror ({1}) {
    draw_vof ("f1", lw = 2);
    draw_vof ("f2", lw = 2);
    draw_vof ("f3", lw = 2);
  }  
  cells();
  save ("mov_lens.mp4");
  
  clear();
  view (fov = 1.54689,  tx = 0.0617166, ty = -0.502029);
  draw_vof ("f1", lw = 2, lc = {1,0,0});
  draw_vof ("f2", lw = 2, lc = {0,1,0});
  draw_vof ("f3", lw = 2, lc = {0,1,1});
  cells();
  save ("zoom.mp4");
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

/**
![Animation of the interfaces and adaptive mesh.](oil_lens/mov_lens.mp4)

![Zoom on the triple point.](oil_lens/zoom.mp4)
*/