/**
 
# Liquid lens with graivity
The test case is token from [Zhang et al, 2016](#zhang_diffuse_2016)

![Time evolution of the interface for Bo=3.](liquid_lens_gravity/movie.mp4)


*/

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "reduced-3p.h"
#include "view.h"

#define L 			5.0 	//size of the square domain
#define R 			0.5		// Radius of the bubbles
#define OH 			0.25 	// Ohnesorge number	
#define BO          3.0	// Bond number
#define MAXLEVEL 	6
#define MINLEVEL 	3

vertex scalar phi[]; 				//distance function used to define the fractions

int main(){
  size (L);
#if AXI
#else
  origin (-L/2.0, -L/2.0);
#endif
  N = 1 << MAXLEVEL;
  init_grid (N);
  double sigma12 	= 1.;
  double sigma23 	= 1.;
  double sigma31 	= 1.;
  f1.sigma 		= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 		= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 		= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 			= 1.0;
  rho2 			= 10.;
  rho3 			= 5.;
  mu1 			= OH*1.0;
  mu2				= OH*1.0;
  mu3 			= OH*1.0;
  G.y = - BO;
  run();
  // movie generation is broken by the return below. This is a bug of Basilisk.
  //  return 0; 
}

event init (t = 0){
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

event end(t=15.){}

/**
## References

~~~bib

@article{zhang_diffuse_2016,
	title = {Diffuse interface simulation of ternary fluids in contact with solid},
	volume = {309},
	issn = {00219991},
	url = {https://linkinghub.elsevier.com/retrieve/pii/S0021999115008797},
	doi = {10.1016/j.jcp.2015.12.054},
	language = {en},
	urldate = {2021-07-10},
	journal = {Journal of Computational Physics},
	author = {Zhang, Chun-Yu and Ding, Hang and Gao, Peng and Wu, Yan-Ling},
	month = mar,
	year = {2016},
	pages = {37--51},
	file = {Zhang 等。 - 2016 - Diffuse interface simulation of ternary fluids in .pdf:/home/xue/Zotero/storage/ZXJ6LHSS/Zhang 等。 - 2016 - Diffuse interface simulation of ternary fluids in .pdf:application/pdf},
}

~~~

*/
