/**
# Simple test of Basilisk View

Ran in parallel on four MPI processes. */

#include "fractions.h"
#include "view.h"

int main() {

  /**
  We first define a volume fraction field. */
  
  init_grid (16);
  origin (-0.5,-0.5,-0.5);
  scalar f[];
  fraction (f, sq(x) + sq(y) + sq(z) - sq(0.3));

  /**
  Then display it using Basilisk view functions. 
  
  <table>
  <tr>
  <td>![2D Basilisk view](view/out.png)</td>
  <td>![3D Basilisk view](view.3D/out.png)</td>
  </tr>
  <tr>
  <td>2D Basilisk view</td>
  <td>3D Basilisk view</td>
  </tr>
  </table>
  */
  
  view (width = 550, height = 400, tx = -0.2, cache = 10);
  box();
  draw_vof ("f");
  cells();
  squares ("x < 0 && y < 0 ? sin(6*pi*x)*cos(8*pi*y) : nodata", spread = -1,
	   cbar = true, border = true, pos = {0.47, - 0.7},
	   label = "A colorbar!", mid = true, format = "%9.5f", levels = 10);
  squares ("(f[0,1] - f[0,-1])/(2.*Delta)", spread = -1);
#if dimension == 2  
  isoline ("sqrt(x^2 + y^2)", n = 10, spread = -1, lc = {1,0,0});
#endif
  save ("out.png");

  /**
  A few more files just for testing. */
  
  output_facets (f, qerr);
  dump (file = "dump");
}
