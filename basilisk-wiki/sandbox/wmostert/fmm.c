/**
 This is a basic 2D implementation of the Noumir et al (2015) method for fast marching method in GSD. 
Wouter Mostert
University of Oxford
Begun 20 July 2022
*/

/** choose a grid type, Cartesian */
#include "grid/cartesian.h"

/** set up basic running plumbing */
#include "run.h"

/** and useful apparatus for the levelsetting */
#include "fractions.h"

/** now set up the field of alpha (surface levelset) and Mach number */
scalar al[];
scalar M[];
double eps = 1e-6;

/** User initialisation happens here. For this default case, set up a basic Riemann problem, sec 4.2.1 Noumir et al (2015).
   The domain is $[0, 5] \times [0, 5]$.
   The Riemann problem is a shock with a piecewise linear profile, with vertex at $(x,y) = (0.5, 2.5)$. For $x=0.5, y<2.5$ the shock has profile $x = 0.5$ with $M_l=10$ (normal angle $\theta_l=0$ to the horizontal). For $y \geq 2.5$ the shock has a normal angle $\theta_r = 20.97^\circ$ to the horizontal and $M_l=8.5$. The shock surface is therefore at an angle $\chi_r = 20.97^\circ + 90^\circ$. 
   The gradient of the shock profile is $dy/dx = \tan \chi_r$.
 */

int main() {
  origin(0., 0.);
  size(5.);
  init_grid(512);
  run();
  return 0;
}

double deg2rad(double deg){
  return deg*pi/180.;
}

double surface(double x, double y, double chir){
  double x0 = 0.5;
  double y0 = 2.5;
  double m =tan(chir);
  double c = y0 - m*x0;
  if (y > 2.5)
    return (y - c)/m - x;
  else return x0 - x;
}
  
event init (i = 0) {
  printf("hello/n");
  double chir = deg2rad(20.97 + 90.);
  fraction (al, surface(x, y, chir));
/**
  Now cells valued 0 or 1 (i.e. not "on the interface") are treated as far: */
  foreach() {
    if (al[] > eps || al[] < (1. - eps) )
      al[] = HUGE;
    else al[] = 0.;
  }
  //output_ppm(al, n = 256, file="al.png");
}

event end (i = 1) {
  printf("done/n");
}



//![Contour of alpha](fmm/al.png)