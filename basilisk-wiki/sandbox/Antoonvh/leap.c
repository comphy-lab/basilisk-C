/**
   ![Leap into the Deep is an artwork designed by [Tijs Rooijakkers](https://tijsrooijakkers.nl/portfolio/leap-into-the-deep-tu-e/), representing a bath-tub vortex. It is located on the TU/e campus, near the Gemini building (right in picture).](https://tijsrooijakkers.nl/wp-content/uploads/2019/07/Enghel-Tijs-Rooijakkers-BIG-ART-2018-01-kl-1200x550-02.png)

# Leap into the Deep

Here we visualize a parametric curve that looks a bit like it:

![A Bwatch rendering with a rendering of the renovated Gemini building](leap/lid.png)

 */

#include "grid/octree.h"
// Multiple definition of `cross`-product function
#define cross cross_tree
#include "tree/treegen.h"
#undef cross
#include "view.h"
#include "bwatch.h"
#define SECT 200

Branch LiD[SECT];
double t0 = 0, te = pi;

coord leap_into_the_deep (double ti) {
  double fr = ti*(ti - pi)*(ti - 2*pi);
  double phi = 7*ti;
  double fz = 20*sqrt(sin(ti));
  return (coord){fz, fr*cos(phi), fr*sin(phi)}; 
}

int main() {
  L0 = 30;
  Z0 = Y0 = -L0/2;
  init_grid (64);
  for (int i = 0; i < SECT; i++) {
    double ts = t0 + i*(te - t0)/(SECT);
    double tse = t0 + (i + 1)*(te - t0)/(SECT);
    
    LiD[i] = (Branch){.start = leap_into_the_deep (ts),
		    .end = leap_into_the_deep (tse),
		    .R = 0.5, .parent = 1};
  }
  nb = SECT;
  scalar cs[];
  cs.refine = cs.prolongation = fraction_refine;
  face vector fs[];

  
  view (phi = pi/2);

  tree_interface_adapt (LiD, cs, fs, alist = {cs}, crit = (double[]){0.01}, maxlevel = 9, stop = 200);
  squares ("cs", n = {0,0,1}, alpha = 5);
  save ("fig.png");

  watch (O = {10, 40, 4}, poi = {9.5, 0, 0.5}, up = {1, 1e-6, 1e-6},
	 fov = 40, nx = 1024, ny = 1024);
  foreach()
    cs[] = 1 - cs[];
  equiplane (cs, vof = true, mat = {.col = {230, 230, 230}, .R = 0.3});
  //sketch_vof (cs, mat = {.col = {230, 230, 230}, .R = 0.3});
  system ("wget -q \
https://assets.w3.tue.nl/w/fileadmin/content/universiteit/gemini/csm_Visual_Gemini_zuid_noord_697b7b493a.jpg");
  image ("csm_Visual_Gemini_zuid_noord_697b7b493a.jpg", 
  	 n = {.0, 1, 0}, alpha = -14, res = 650);
  system ("wget -q \
https://www.ctm.co.za/media/catalog/product/cache/962f0fa93f6ebc4f041fe71bedef3e3f/a/r/artistry_terrazzo.jpg");
  image ("artistry_terrazzo.jpg", 
  	 n = {.1, 0.001, 0}, alpha = .01, res = 650);
  sphere (R = 50, mat = {.dull = true});
  store (fopen ("nja.ppm", "w"));
  system ("convert nja.ppm -resize 80% lid.png");
}
