/**
#Three-phasic splash of a water drop over an oil film on a deep water pool.
*/

/**
we want to simulate a three-phasic splash under the experimental conditions carried out by Guy-Jean Michon at the Institut d'Alembert.
*/

#include "axi.h"
#include "./centered.h"
#include "./three-phase.h"
#include "./conserving3p.h"
#include "tension.h"
#include "view.h"


#define radius 0.2
#define Dre 2.8e-3
#define Vre 1.935
#define Mu_re2 1.0e-3
#define Mu_re3 2.67e-2
#define Rho2_re 1000.
#define Rho3_re 900. 
#define Rho1_re 1.21
#define Re (Rho2_re*Vre*Dre/Mu_re2)
#define Sig_re2 72e-3
#define Sig_re3 22e-3
#define Sig_re32 50e-3
#define We2 (Rho2_re*Vre*Vre*Dre/Sig_re2)
#define We3 (Rho2_re*Vre*Vre*Dre/Sig_re3)
#define We32 (Rho2_re*Vre*Vre*Dre/Sig_re32)
#define SIGMA2 (2.*radius/We2)/2.
#define SIGMA3 (2.*radius/We3)/2.
#define SIGMA32 (2.*radius/We32)/2.

#define Lref (Dre/(2.*radius))   
#define gad (Lref*9.81/(Vre*Vre))

#define vrho2 (1.)
#define vrho3 (Rho3_re/Rho2_re)
#define vrho1 (Rho1_re/Rho2_re)
#define vmu2 (2.*radius/Re)
#define vmu3 (2.*radius/Re*(Mu_re3/Mu_re2))
#define vmu1 (2.*radius/Re*(1.81e-5/Mu_re2))

#define ycc 0.
#define hpool (4.*radius*2.)
#define hlame (0.707*radius*2.*1.)
#define xcc0 0.
#define xcc2 (xcc0 + hpool)
#define xcc3 (xcc0 + hpool + hlame)
#define hgout (0.25*radius*2.)
#define xcc (xcc3 + hgout + radius)

#define circle(x,y) (sq(x - xcc) + sq(y - ycc))
#define drop(x,y) (sq(radius) - circle(x,y))
#define pool(x) (xcc2 -x)
#define lame(x) (xcc3 -x)*(x-xcc2)

#define epi radius/8.
#define epi2 radius/8.
#define rfdrop(x,y) ((sq(radius-epi) < circle(x,y)) && (sq(radius+epi) > circle(x,y)))
#define rfpool(x) (((xcc2-epi) < x) && ((xcc2+epi) > x))
#define rflame(x) (((xcc3-epi2) < x) && ((xcc3+epi2) > x))

int maxlevel = 9;
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

  
  scalar f4[], f5[], f6[];

  refine ((rflame(x) || rfdrop(x,y) || rfpool(x)) && (level < maxlevel));
  
  fraction (f4, pool(x));
  fraction (f5, lame(x));
  fraction (f6, drop(x,y));
  foreach() {
    f2[] = f4[] + f6[];
    f3[] = f5[];
    f1[] = clamp(1.-f2[]-f3[],0.,1.);
    rhov[] = rho(f2[],f3[]);
    u.x[] = -f6[];
  }
  boundary ((scalar *){f1,f2,f3,rhov});
}

event movie (t =0.; t += 0.02; t <= 3.00) {
  clear();
  view (fov = 9.93351, quat = {0,0,-0.707,0.707}, tx = 0.300417, ty = -0.432814, bg = {1,1,1}, width = 1024, height = 768, samples = 1);
  draw_vof ("f2",filled=1, fc = {0.,0.7,0.7});
  draw_vof ("f3",filled=1, fc = {0.7,0.7,0.});
  cells();
  save("splash_mesh.mp4");
  
  clear();
  view (fov = 9.93351, quat = {0,0,-0.707,0.707}, tx = 0.300417, ty = -0.432814, bg = {1,1,1}, width = 1024, height = 768, samples = 1);
  draw_vof ("f2",filled=1, fc = {0.,0.7,0.7});
  draw_vof ("f3",filled=1, fc = {0.7,0.7,0.});
  save("splash.mp4");
}

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
  }
  boundary ({ff2,ff3});    
  adapt_wavelet ({ff2,ff3,u,f2,f3}, (double[]){0.01,0.01,uem,uem,fem,fem}, maxlevel, minlevel);  
  event ("properties");
}

/**
At level 9 the oil film is not refined correctly
*/

/**
![Animation of the three phases - level 9](splash3ph/splash.mp4)

![Three phases with mesh - level 9](splash3ph/splash_mesh.mp4)
*/

/**
At level 12 the refinement is correct
*/

/**
![Animation of three phases - level 12](final_splash_color.mp4)

![Three phases with mesh - level 12](final_splash_mesh.mp4)
*/