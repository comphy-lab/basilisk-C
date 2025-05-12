#include "radial.h"
#include "run.h"
#include "view.h"
//#include "display.h"

int main() {
  dtheta=2*pi;
  init_grid (32);
  run();
}

void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}


event init (i = 0) {
  refine (sq(x - 0.5) + sq(y - 0.5) < sq(0.3) && level < 7);
  foreach() {
        fprintf (stderr, "%g %g %g\n", r, theta*180./M_PI, cm[]);
    assert (fabs(cm[] - r*dtheta/L0) < 1e-12);
  }

  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.013, ty = -0.011, tz = -3.829,
      width = 1548, height = 936);
  box ();
  cells ();
  squares(color = "cm");
  save("scale-factor.png");
     
  scalar angle[];
  foreach()
    angle[] = y*dtheta/L0*180./M_PI;
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000}, map = radial,
      fov = 30, near = 0.01, far = 1000,
      tx = 0.013, ty = -0.011, tz = -3.829,
      width = 1548, height = 936);
  box ();
  cells ();
  squares(color = "angle", min=0, max=360);
  save("theta.png");  
}

/**
 ![cm](radial2/scale-factor.png)
 ![angle in degree](radial2/theta.png)
 
How to plot in polar coordinate with gnuplot
  ~~~gnuplot (requires gnuplot>=5.4.0, because of smooth zsort)
### plot circle segments to mimic polar heatmap
reset session

FILE = "log"

da = 5.625        # grid angle step
dr = 0.015625      # grid radius step

set table $Sorted
    plot FILE u 1:2:(-$3):1 smooth zsort lc var       # sort by decreasing radius
set table $Grid
    plot $Sorted u 1:2:3:($1/dr==int($1/dr)) w table  # check if on the grid
unset table

set size ratio -1
set xrange[-1:1]
set yrange[-1:1]
set style fill transparent solid 1.0 border lc "black"
set key noautotitle
set palette defined (0 "#00d0ff", 1 "#40fffc", 2 "#b7ffc4", 3 "#f3ff65", 4 "#ffff00")

plot FILE  u (NaN):(NaN):3 w p palette z, \
     $Grid u (0):(0):(g=($4==1?1:0.5), $1+dr*g):($2-da*g):($2+da*g):3 w circles lc palette z
### end of script
  ~~~
  */