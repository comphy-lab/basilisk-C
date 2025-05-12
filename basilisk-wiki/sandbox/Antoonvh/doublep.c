/**
# Doubly-periodic shear layer instability

The results can be compared against Fig. 7 of [Hokpunna and Manhart, 2010](#hokpunna2010).

![Vorticity isolines](doublep/mov.mp4)

![Reference result](doublep/cropped.png)
*/
//#include "grid/multigrid.h" //Todo
#include "nsf4t.h"
#include "view.h"
scalar * tracers = NULL;

double Re = 1e4;

int main() {
  foreach_dimension()
    periodic (left);
  const scalar muz[] = 1./Re;
  nu = muz;
  N = 256;
  run();
}
double sigma = 30, eps = 0.05;

double u_x (double x, double y) {
  if (y <= 0.5)
    return tanh(sigma*(y - 0.25));
  else
    return tanh(sigma*(0.75 - y));
}

double u_y (double x, double y) {
  return eps*sin(2.*pi*x);
}

event init (t = 0) {
  TOLERANCE = 1e-5;
  foreach_face()
    u.x[] = Gauss6_x(x, y, Delta, u_x);
}

event mov (t += 0.025) {
  scalar omg[], omgv[]; //Cell averages and vertex-point vorticity
  vorticityf (u, omg);
  coord f = {1, -1};
  foreach() {
    omgv[] = 0;
    foreach_dimension()
      omgv[] += f.x*(15*(u.y[] - u.y[-1])
		     - u.y[1] + u.y[-2])/(12.*Delta);
  }
  view (tx = -0.5, ty = -0.5, samples = 4);
  for (double val = -36; val <= 37 ; val += 6.) 
    isoline("omgv", val, lw = 2);
  squares ("omg");
  box();
  save ("mov.mp4");
}

/**
We stop and generate the reference result image from the available
image. */

event stop (t = 1.2) {
  system ("wget https://ars.els-cdn.com/content/image/1-s2.0-S0021999110003141-gr7.jpg && \
  convert 1-s2.0-S0021999110003141-gr7.jpg -crop 300x300+280+600 -resize 180% cropped.png");	  
}

/**
## Reference 

~~~bib
@article{hokpunna2010,
title = "Compact fourth-order finite volume method for numerical solutions of Navier–Stokes equations on staggered grids",
journal = "Journal of Computational Physics",
volume = "229",
number = "20",
pages = "7545 - 7570",
year = "2010",
issn = "0021-9991",
doi = "https://doi.org/10.1016/j.jcp.2010.05.042",
url = "http://www.sciencedirect.com/science/article/pii/S0021999110003141",
author = "Arpiruk Hokpunna and Michael Manhart"
}
~~~
*/
