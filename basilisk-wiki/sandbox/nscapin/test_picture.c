#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "lambda2.h"

double h_ = 1; 
double k_ = 4; 
double ak = 0.3; 

int main(int argc, char *argv[]) {

  L0 = 2*pi;
  origin (-L0/2., 0, -L0/2.);

  periodic (right);
  periodic (front);
  
  init_grid (1 << 5);
  run();

}

/** 
   Specify the interface shape using a third-order Stokes wave. */

double WaveProfile (double x, double y) {

  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;

}

event init (i = 0) {
  do {
    fraction (f, WaveProfile(x,z)-y);
  }
  //while (adapt_wavelet ({f}, (double[]){1e-3}, 5, 4).nf); // if not adapting anymore, return zero
  while (0); // no adaptation
}

/** 
   **First issue when running in parallel**

   If I remove the box():

   1. If N_proc<7: the picture is visible but the interface is broken;
   2. If N_proc>=7: the picture is not correctly output, i.e. white screen 
   
   ![Test 1](test_picture/3D_p1.png)
*/

event pic_p1 (i = 0) {
  clear();
  //box(); // if I remove the box(), it does not output anything
  view(tx = -0.125, ty = -0.25, fov = 40, 
       samples = 4, bg = {1,1,1}, camera = "iso");
  draw_vof ("f");
  cells (n = {0,0,1}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  save ("3D_p1.png");
}

/** 
   **Second issue when running in parallel**

   If I keep the box() and I specify the width and height in view():

   1. If N_proc<7: the picture is visible but the interface is broken;
   2. If N_proc>=7: the picture is not correctly output, i.e. white screen 

   ![Test 2](test_picture/3D_p2.png)
*/

event pic_p2 (i = 0) {
  clear();
  box();
  view(tx = -0.125, ty = -0.25, fov = 40, width = 1200, height = 1200, 
  //view(tx = -0.125, ty = -0.25, fov = 40, 
       samples = 4, bg = {1,1,1}, camera = "iso");
  draw_vof ("f");
  cells (n = {0,0,1}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  save ("3D_p2.png");
}
