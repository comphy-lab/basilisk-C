/**
# Dynamic test 2 for liquid lens test case

![Time evolution of the interface.](liquid_lens_dynamic2/movie.mp4)
*/

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

#define L 			2.5 	//size of the square domain
#define LL 			2	//radius of the circle (phase 3)
#define R 			0.047325
#define OH 			0.1 	// Ohnesorge number	
#define MAXLEVEL        6
#define MINLEVEL 	3
#define M1 			50./pi
#define M2 			5.

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
  double sigma23 = 20.;
  double sigma31 = 20.;
  f1.sigma 	= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 	= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 	= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 		= 1.0;
  rho2 		= 1.0;
  rho3 		=  1.0;
  mu1 		= OH*1.0;
  mu2		= OH*1.0;
  mu3 		= OH*1.0;
  run();
  // movie generation is broken by the return below. This is a bug of Basilisk.
  //  return 0; 
}

event init (t = 0){
  if(!restore("dumpfile")){
    foreach_vertex(){
      if(x < - LL/2.)
	phi[] = sq(R) - sq(x + LL/2.) - sq(y);
      else if(x > LL/2.)
	phi[] = sq(R) - sq(x - LL/2.) - sq(y);			
      else
	phi[] = R - fabs(y);
    }
    fractions(phi, f3);
    foreach(){
      if(y> 0 ){
	f1[] = 1. - f3[];
	f2[] = 0.;
      }
      else{
	f2[] = 1. - f3[];
	f1[] = 0.;
      }
    }
    boundary({f1, f2, f3});
    restriction({f1, f2, f3});
  }
  else{}
}


event adapt(i++){
  adapt_wavelet({f1, f2, f3, u.x, u.y}, (double[]) {1.0e-9, 1.0e-9, 1.0e-9, 1.0e-2, 1.0e-2}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}


event movie(t+=0.01){
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


event end(t = 3.){}
