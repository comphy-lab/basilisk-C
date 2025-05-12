/**
# Two drops on solid 
The test case is token from [Zhang et al, 2016](#zhang_diffuse_2016)

![Time evolution of the interface with all the contact angle 90 degrees.](two_drop_on_solid/movie.mp4)


*/


#include "navier-stokes/centered.h"
#include "contact.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

#define L 			4.0 	//size of the square domain
#define R 			1.0		// Radius of the bubbles
#define OH 			0.1 	// Ohnesorge number	
#define MAXLEVEL 	6
#define MINLEVEL 	3

vertex scalar phi[]; 				//distance function used to define the fractions

vector h1[];
vector h2[];
vector h3[];

double theta1 = 90.;
double theta2 = 90;
double theta3 = 90;

h1.t[bottom] = contact_angle (x<0? (180.-theta1)*pi/180.:(180.-theta3)*pi/180. );
h2.t[bottom] = contact_angle (fabs(x)<R/2.? theta2*pi/180.:theta1*pi/180. );
h3.t[bottom] = contact_angle (fabs(x)<R/2.? (180.-theta2)*pi/180.:theta3*pi/180. );
u.t[bottom] = dirichlet(0.);


int main(){
  size (L);
#if AXI
#else
  origin (-L/2.0, -L/2.0);
#endif
  N = 1 << MAXLEVEL;
  init_grid (N);
  double sigma12 = 1.0;
  double sigma23 = 1.0;
  double sigma31 = 1.0;
  f1.sigma 	= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 	= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 	= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 		= 1.0;
  rho2 		= 1.0;
  rho3 		= 1.0;
  mu1 		= OH*1.0;
  mu2		= OH*1.0;
  mu3 		= OH*1.0;
  f1.height = h1;
  f2.height = h2;
  f3.height = h3;

  run();
  // movie generation is broken by the return below. This is a bug of Basilisk.
  //  return 0; 
}

event init (t = 0){
  if(!restore("dumpfile")){
    foreach_vertex()
      phi[] = ((sq(x) + sq(y+L/2.) + sq(z))) - sq(R);
    boundary({phi});
    restriction({phi});
    fractions(phi, f1);
    foreach(){
      if (x<0){
	f3[] = 0.;
	f2[] = 1. - f1[];
      }
      else{
	f2[] = 0.;
	f3[] = 1. - f1[];
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


event end(t=10.){}

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
