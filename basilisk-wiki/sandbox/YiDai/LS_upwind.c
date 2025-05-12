/**
# Simple Level set method trial using upwinding 
 */

// #include "grid/octree.h"
#include "fractions.h"
#include "view.h"
#include "run.h"

scalar f[], df[];
double U = 0.5;

int main() {
  /**
  We first define a volume fraction field. */
  L0 = 10;
  // init_grid (64);
  N = 1 << 8;
  origin (-0.5*L0,-0.5*L0,-0.5*L0);
  // origin (0,0,0);
  DT = L0/N;
  run();
}

event init(t = 0){
  fraction (f, sq(x) + sq(y) + sq(z) - sq(1));
  boundary ({f});
}

event advec (i++) {
  double dt = DT;
  dt = dtnext (dt);
  foreach()
    df[] =  U * (f[] - f[-1])/Delta;
  foreach()
    f[] = f[] + df[] * dt;
  boundary ({f});
}

event printdata (t = 0; t <= 0.3; t += 0.01) {
  static FILE * fp = fopen ("case0.dat","w");
  for (double y = -L0/2; y<L0/2; y+=L0/N){
    fprintf (fp, "%g %g\n", y, interpolate(f, 0, y));}
  fprintf (fp, "\n \n");
  fflush (fp);
}

event watch(t = 0; t <= 0.3; t += 0.01){
  // view (width = 800, height = 800, tx = -0.5, ty = -0.5, theta = -0.3, phi = 0.5, cache = 10);
  view (fov = 23, width = 800, height = 800, tx = 0, ty = 0, theta = 0, phi = 0);
  box();
  draw_vof ("f");
  save ("out.mp4");
}

/**
you can see the break line of the curve
![interface evolution](LS_upwind/out.mp4)
**/
