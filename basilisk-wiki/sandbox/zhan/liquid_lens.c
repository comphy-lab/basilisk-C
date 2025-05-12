/**
# Liquid lens test case

In this test, Phase 3 is sandwiched by phase 1 (top) and phase 2 (bottom).
Under the effect of the surface tension, we will have analytical solution for the length of the liquid lens.

* 2D
$$
D = 2\sqrt{\frac{A}{\frac{1}{\sin\theta_1}\left(\frac{\theta_1}{\sin\theta_1}-\cos\theta_1\right)+\frac{1}{\sin\theta_2}\left(\frac{\theta_2}{\sin\theta_2}-\cos\theta_2\right)}}
$$

* AXI
$$
D = \sqrt[3]{\frac{V}{\pi\left( \frac{  \left(1-\cos\theta_1\right)^2\left(2+\cos\theta_1\right)}{24\sin^3\theta_1}  +     \frac{  \left(1-\cos\theta_2\right)^2\left(2+\cos\theta_2\right)}{24\sin^3\theta_2}  \right)}}.
$$

![The three phase system can be simulated with three tracers of VOF method.](liquid_lens/movie.mp4)

 */

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "norm-3p.h"
#include "view.h"


#define L 		2.5 	//size of the square domain
#define R 		.5		//radius of the circle (phase 3)
#define OH 		.25	// Ohnesorge number	
#define MAXLEVEL 	6
#define MINLEVEL 	3

vertex scalar phi[]; 	     //distance function used to define the fractions

FILE * fp;

int main(){
  size (L);
  origin (-L/2.0, -L/2.0);
  N = 1 << MAXLEVEL;
  init_grid (N);
  double sigma12 	= 1.0;
  double sigma23 	= 1.0;
  double sigma31 	= 1.0;
  f1.sigma 		= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 		= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 		= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 			= 1.0;
  rho2 			= 1.0;
  rho3 			= 1.0;
  mu1 			= OH * 1.0;
  mu2			= OH * 1.0;
  mu3 			= OH * 1.0;
  run();
  return 0; 
}

event init (t = 0){
  if(!restore("dumpfilekk")){
    refine(sq(x) + sq(y) + sq(z) < sq(1.1 * R) && sq(x) + sq(y) + sq(z) > sq(0.9 * R) && level < MAXLEVEL);
    refine((y < 0.1 * L/2.0 && y > -0.1 * L/2.0 && x > R && x < -R) && level < MAXLEVEL);
    foreach_vertex()
      phi[] = sq(R) - ((sq(x) + sq(y) + sq(z)));
    boundary({phi});
    restriction({phi});
    fractions(phi, f3);
    foreach(){
      if (y <= 0){
	f2[] = 1. - f3[];
	f1[] = 0.;
      }
      else{
	f1[] = 1. - f3[];
	f2[] = 0.;
      }
    }
    boundary({f1, f2, f3});
    restriction({f1, f2, f3});
  }
  else{}
}

event adapt(i++){
  adapt_wavelet({f1, f2, f3, u.x, u.y}, (double[]) {1.0e-9, 1.0e-9, 1.0e-9, 1.0e-3, 1.0e-3}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event movie(t+=0.1){
  clear();
  box();
  cells();
  draw_vof("f1");
  draw_vof("f2", filled = 1, fc = {0.9, 0.4, 0.2});  
  draw_vof("f2");
  draw_vof("f3", filled = 1, fc = {0.2, 0.4, 0.9});
  draw_vof("f3");
  save("movie.mp4");
}

event logfile (t=end){
  coord p1, p2;
  locate_triple_point(&p1, f1, f2, f3, -L/2., 0., -L/2., L/2., 1.1);
  locate_triple_point(&p2, f1, f2, f3, 0., L/2., -L/2., L/2., 1.1);
  output_cells(stdout);
  output_facets(f1, stderr);
  output_facets(f2, stderr);
  output_facets(f3, stderr);
  char filename[80];
  sprintf(filename, "triple_point");
  fp = fopen(filename, "w");
  fprintf(fp, "%g %g\n%g %g\n", p1.x, p1.y, p2.x, p2.y);
}

event end(t = 5.0){}
