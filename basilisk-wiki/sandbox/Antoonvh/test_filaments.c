/**
# Test for `filaments.h`

We want to test the initialization of a `Lamb_Oseen` vorticity field
using a varying number of filament sections for its discretization.

![Vortex tube: The induvidual segments can be notable. The color
 represent the error](test_filaments/omg.mp4)

~~~gnuplot Use about one segment per core radius
set xr [0.09:11]
set logscale xy
set grid
set size square
set xlabel 'segments per core diameter'
set ylabel 'error'
plot 'out'
~~~

~~~gnuplot Average vorticity
reset
set xr [0.005: 0.015]
set grid
set size square
set xlabel 'Average vorticity'
set ylabel 'z'
plot 'prof12' u 2:1 t '12 segments', 'prof22' u 2:1 t '22 segements', 'prof43' u 2:1 w l lw 2 t '43 segements'
~~~

The convergence with the number of segments is limited. This is due to
the finite length of the vortex tube and the limited grid resolution.
*/
#include "grid/multigrid3D.h"
#include "filaments.h"
#include "profile6.h"
#include "view.h"
double xp = 0.1, yp = -0.2;

coord line (double t) {
  coord C;
  C.x = xp;
  C.y = yp;
  C.z = t;
  return C;
}

int main () {
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2;
  init_grid (64);
  vector omega[];
  for (int nsegs = 5; nsegs < 200; nsegs *= sqrt(2)) {
    get_vor_vector (omega, line, -2*L0, 2*L0, nsegs, Lamb_Oseen);
    char fname[99];
    sprintf (fname, "prof%d", nsegs);
    profile ({omega.z}, z, fname);
    double err = 0;
    scalar error[];
    foreach() {
      double val = 0;
      foreach_child()
	val += exp(-sq(x - xp) - sq(y - yp))/(8*pi);
      err += dv()*fabs(val - omega.z[]);
      error[] = (val - omega.z[]);
    }
    printf ("%g %g\n", nsegs/(4*L0), err);
    view (theta = pi/2, psi = pi/2);
    isosurface ("omega.z", exp(-1)/pi, color = "error", min = -0.2, max = 0.2);
    save ("omg.mp4", opt = "-r 3");
  }
}


  

