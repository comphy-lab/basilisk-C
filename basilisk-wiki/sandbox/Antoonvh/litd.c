/**
 ![Leap into the Deep is an artwork designed by [Tijs Rooijakkers](https://tijsrooijakkers.nl/portfolio/leap-into-the-deep-tu-e/), representing a bath-tub vortex. It is located on the TU/e campus, near the Gemini building (right in picture).](https://tijsrooijakkers.nl/wp-content/uploads/2019/07/Enghel-Tijs-Rooijakkers-BIG-ART-2018-01-kl-1200x550-02.png)

# Leap into the deep

 Here we illustrate the tangent, normal and binormal vector of a
 curve, based on an approximation of the Leap into the deep artwork.

 ![One may use a local coordinate system when walking along a paramterized curve](litd/lid.mp4)
 */
#include "grid/octree.h"
#include "view.h"
#include "some_primitives.h"

// Paramterization
coord leap_into_the_deep (double ti) {
  double fr = ti*(ti - pi)*(ti - 2*pi);
  double phi = 7*ti;
  double fz = 20*sqrt(max(sin(ti), 0));
  return (coord){fr*sin(phi), fz, fr*cos(phi)}; 
}

coord knot (double theta) {
  coord C;
  C.x = sin(theta) + 2.*sin(2.*theta);
  C.y = cos(theta) - 2.*cos(2.*theta);
  C.z = -sin(3.*theta) - 4;
  return C;
}

int main() {
  L0 = 10;
  X0 = Z0 = -L0/2.;
  init_grid (16); // Tile the floor
  float fc[3] = {.8, .8, .8};
  double dt = pi/500;
  for (t = pi/100; t <= 99*pi/100; t += dt) {
    double Deltat = dt/10; // Vector rejection is very sensitive 
    view (ty = -1.2, fov = 50, phi = -0.1, theta = 0.2);
    draw_paramterization (0.1, 0., pi, leap_into_the_deep, 500, fc);
    coord tangent = coord_diff (leap_into_the_deep (t ), leap_into_the_deep (t - Deltat));
    coord tangent_next = coord_diff (leap_into_the_deep (t + Deltat), leap_into_the_deep (t));
    coord dtdt = coord_diff (tangent_next, tangent); // Quasi normal
    foreach_dimension() // Central approximation
      tangent.x = (tangent.x + tangent_next.x)/2.;
    normalize (&tangent);
    coord normal = rejection (tangent, dtdt); // normal
    normalize (&normal);
    coord binormal = cross (tangent, dtdt);
    draw_arrow (tangent, leap_into_the_deep (t), rad = .2, fc = {0.8, 0.1, 0.1});
    draw_arrow (normal, leap_into_the_deep (t), rad = .2, fc = {0.8, 0.8, 0.1});
    draw_arrow (binormal, leap_into_the_deep (t), rad = .2, fc = {0.1, 0.8, 0.1});
    cells (n = {0,1,0}, alpha = 0.01);
    draw_string (" Tangent,", lc = {0.8, 0.1, 0.1}, size = 32, lw = 3);
    draw_string ("          Normal, ", lc = {0.8, 0.8, 0.1}, size = 32, lw = 3);
    draw_string ("                  Binormal", lc = {0.1, 0.8, 0.1}, size = 32, lw = 3);
    save ("lid.mp4");
    printf ("%g %% \n", 100*t / (pi));
  }
    /**
    
    ## Bonus
    
    We repeat this excersize for a [trefoil](trefoil4.c), which is periodic.
    
    ![The trefoil goemetry](litd/trefoil.mp4)(loop)
    */
  dt = 2*pi/500.;
  for (t = 0; t < 2*pi; t += dt) {
    view (ty = 0, fov = 15, phi = 0.2, theta = 0);
    double Deltat = dt/10; // Vector rejection is very sensitive 
    //draw_paramterization (0.1, 0., 2*pi, knot, 500, fc);
    coord tangent = coord_diff (knot (t), knot (t - Deltat));
    coord tangent_next = coord_diff (knot (t + Deltat), knot (t));
    coord dtdt = coord_diff (tangent_next, tangent); // Quasi normal
    foreach_dimension() // Central approximation
      tangent.x = (tangent.x + tangent_next.x)/2.;
    normalize (&tangent);
    coord normal = rejection (tangent, dtdt); // normal
    normalize (&normal);
    coord binormal = cross (tangent, dtdt);
    draw_arrow (tangent, knot(t), len = 1, rad = .15, fc = {0.8, 0.1, 0.1});
    draw_arrow (normal, knot(t), len = 1, rad = .15, fc = {0.8, 0.8, 0.1});
    draw_arrow (binormal, knot(t), len = 1, rad = .15, fc = {0.1, 0.8, 0.1});
    draw_string (" Tangent,", lc = {0.8, 0.1, 0.1}, size = 32, lw = 3);
    draw_string ("          Normal, ", lc = {0.8, 0.8, 0.1}, size = 32, lw = 3);
    draw_string ("                  Binormal", lc = {0.1, 0.8, 0.1}, size = 32, lw = 3);
    draw_paramterization (0.1, 0., 2*pi, knot, 200, fc);
    save ("trefoil.mp4");
    printf ("%g %% \n", 100*t / (2*pi));
  }
  
}
