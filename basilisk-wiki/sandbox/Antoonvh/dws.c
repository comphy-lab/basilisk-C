/**
# A droplet water slide

A heavy droplet is placed on a slippery slope and descends into a
water pool. The drop does not mix with water.

![A water wave is created in 2D](dws/f.mp4)

<div class="figure"> <video controls="" preload="metadata"> 
<source
src="https://surfdrive.surf.nl/files/index.php/s/01xfat5oFX20B75/download"
type="video/mp4"> Your browser does not support the video
tag. </video> <p class="caption"> In 3D as well (Video via
surfdrive) </p> </div>
*/
//#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
// Using Nelson Joubert's three-phase formalism
#include "../joubert/three-phase.h"
#include "view.h"

int maxlevel = 9;
// Slope angle, drop size and distance to the pool
double angle = 20.*pi/180., Rd = 1, Dist = 5;

int main () {
  L0 = 25.;
  mu3 = 0.1;   // Drop
  mu2 = 0.01;  // Water
  mu1 = 0.001; // Air
  rho1 = 0.1;  // Air
  rho2 = 1.;   // Water
  rho3 = 1.5;  // Drop
  const face vector av[] = {sin(angle), -cos(angle), 0};
  a = av;
  N = 64;
  run();
}

event init (t = 0) {
  do {
    fraction (f2, tan(angle)*(x - L0/2) - y);
    fraction (f3, sq(Rd) - sq(x - L0/2. + Dist) - sq(y - Y0) - sq(z));
  } while (adapt_wavelet ({f2, f3}, (double[]){0.01, 0.01}, maxlevel).nf);  
  foreach() 
    f1[] = 1. - f2[] - f3[]; 
  DT = 0.01;
}

event adapt (i++) 
  adapt_wavelet ({f1, f2, u}, (double[]){0.01, 0.01, 0.1, 0.1, 0.1}, maxlevel);

event mov (t += 0.1) {
  if (dimension == 2)
    view (fov = 3*Dist, psi = angle, tx = -0.55, ty = -0.05);
  else { 
    view (fov = 4*Dist, psi = angle, tx = -0.55, ty = 0.05,
	  phi = 0.2, phi = 0.1);
    squares ("x + z", n = {0, 1, 0}, min = -2*L0, max = 4*L0);
  }
  box();
  draw_vof ("f3", filled = 1, fc = {205./255., 205./255., 0});
  draw_vof ("f2", filled = 1, fc = {50./255. , 100./255., 205./255.});
  save ("f.mp4");
}

event stop (t = 15);
