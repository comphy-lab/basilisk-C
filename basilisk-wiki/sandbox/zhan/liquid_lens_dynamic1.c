/**
# Dynamic test 1 for liquid lens test case

![Time evolution of the interface.](liquid_lens_dynamic1/movie.mp4)


 */

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

#define L 			2.5 	//size of the square domain
#define R 			0.5		//radius of the circle (phase 3)
#define OH 			0.1 	// Ohnesorge number	
#define MAXLEVEL 	6
#define MINLEVEL 	3
#define M1 			50./pi
#define M2 			5./4.

vertex scalar phi[]; 				//distance function used to define the fractions

int main(){
  size (L);
#if AXI
#else
  origin (-L/2.0, -L/2.0);
#endif
  N = 1 << MAXLEVEL;
  init_grid (N);
  double sigma12 = 1.0;
  double sigma23 = 0.51;
  double sigma31 = 0.51;
  f1.sigma 	= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 	= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 	= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 		= 1.0;
  rho2 		= 1.0;
  rho3 		= 1.0;
  mu1 		= OH * 1.0;
  mu2		= OH * 1.0;
  mu3 		= OH * 1.0;
  run();
  // movie generation is broken by the return below. This is a bug of Basilisk.
  //  return 0; 
}

event init (t = 0){
  refine(((x+5-y) < 1) && ((x+5-y)>-1) && ((x-5-y) > -1) && ((x-5-y) < 1) && level < MAXLEVEL);
  refine((x < 1.1* L/2.0 && x > 0.9 * L/2.0) && level < MAXLEVEL);
  foreach_vertex(){
    if(x < 0)
      phi[] = x+1.25-y;
    else
      phi[] = x-1.25-y;
  }
  fractions(phi,f1);
  foreach_vertex(){
    if(x<0 && y>0)
      phi[] = 	M1*x + M2 - y;
    if(x>=0 && y>0)
      phi[] = - 	M1*x + M2 - y;
    if(x<0 && y<=0)
      phi[] = y + M1*x + M2 ;
    if(x>=0 && y<=0)
      phi[] = y - M1*x + M2 ;

  }
  fractions(phi, f3);
  foreach(){
    if(f3[] + f1[] >= 1.)
      f1[] = 1. - f3[];
  }
  foreach(){
    f2[] = 1.-f1[]-f3[];
  }
  boundary({f1, f2, f3});
  restriction({f1, f2, f3});
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


event end(t = 30.){}
