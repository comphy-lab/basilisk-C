/**
# A Puck

After the succes of 
[Alexis](http://basilisk.fr/sandbox/aberny/csgBool.c) we try to
construct a shape.
*/
#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"

/**
A puck is defined here as the volume between the perpendicular
intersections of two planes and an infinitely long cylinder.

![Its shape corresponds to that of a fan?](puck/puck.mp4)
 */

int main(){
  init_grid(128);
  double R = 0.25;
  double w = 0.05;
  X0 = Y0 = Z0 = -L0/2;
  double xf = 0;
  double yf = 0;
  double zf = 0;
  scalar fan[], cyl[], plane[];
  fan.prolongation = fraction_refine;
  view(theta = 0.1, phi = 0.3);
  for (double dir = 0; dir < 3*pi; dir += 0.01){
    fraction(cyl,  sq(R) - sq(cos(dir)*(x - xf) - sin(dir)*(y - yf)) - sq(z - zf));
    fraction(plane, -cos(dir)*(y - yf - cos(dir)*w/2.) - sin(dir)*(x - xf - sin(dir)*w/2.));
    foreach()
      fan[] = cyl[] * plane[];
    fraction(plane, cos(dir)*(y - yf + cos(dir)*w/2.) + sin(dir)*(x - xf + sin(dir)*w/2.));
    foreach()
      fan[] *= plane[];
    clear();
    boundary({fan});
    cells(alpha = 0);
    draw_vof("fan");
    save("puck.mp4");
    adapt_wavelet({fan}, (double[]){0.01}, 9, 5);
  }
}
