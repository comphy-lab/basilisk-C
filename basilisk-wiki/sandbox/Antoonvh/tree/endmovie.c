/**
# Inspecting a dumpfile with `bview`

The result is presented [elsewhere](fractaltree.c).
*/
#include "grid/octree.h"
#include "utils.h"
#include "view.h"
#include "lambda2.h"
#include "../iso3D.h"

scalar cs[];
vector u[];

int main() {
  printf ("Restore...\n");
  restore ("dump3D150");
  boundary (all); 
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});

  printf ("change angle, translate dist., and iso val\n");
  double thst = -0.8, phist = 0.4;  // theta
  double then = 0, phien = pi/2.;   // phi
  double tyst = -0.15, tyen = 0;    // ty
  double pyst = 20, pyen = 0;       // py (slice)
  double valst = -0.1, valen = -10; // isosurf. val
  double fovst = 12, foven = 8;     // fov
  double txst = 0, txen = -0.1;     // tx
  int jtrans = 100;
  view (bg = {65./256,157./256,217./256}, width = 1080, height = 1080);
  double theta, phi, ty, py, val, tx, fov;
  for (int j = 0; j <= jtrans; j++) {
    theta = thst  + (double)j/jtrans*(then  - thst);
    phi   = phist + (double)j/jtrans*(phien - phist);
    ty    = tyst  + (double)j/jtrans*(tyen  - tyst);
    py    = pyst  + (double)j/jtrans*(pyen  - pyst);
    val   = valst + (double)j/jtrans*(valen - valst);
    fov   = fovst + (double)j/jtrans*(foven - fovst);
    tx    = txst  + (double)j/jtrans*(txen  - txst);
    view (theta = theta, phi = phi, tx = tx, ty = ty, fov = fov);
    translate (y = -py)
      squares ("u.x", n = {0,1,0}, alpha = py,
	       min = -1.1, max = 1.1, map = cool_warm);
    draw_vof("cs", fc = {98./256,78./256,44./256});
    isosurface ("l2", val);
    draw_string ("Changing perspective", 1, lw = 3);
    save ("endmov.mp4");
  }
  
  printf ("Raise the slice\n");
  pyst = pyen;
  pyen = 27;
  jtrans = 70;
  for (int j = 0; j <= jtrans; j++) {
    py = pyst + (double)j/jtrans*(pyen - pyst);
    squares ("u.x", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
    draw_vof("cs", fc = {98./256,78./256,44./256});
    draw_string ("The horizontal slice...", 1, lw = 3);
    save ("endmov.mp4");
  }
  
  printf ("Switch to vorticity_y\n");
  scalar omgy[];
  foreach() 
    omgy[] = cs[] ? (u.x[0,0,1] - u.x[0,0,-1]
		     - u.z[1] + u.z[-1])/(2.*Delta) : 0;  
  jtrans = 40;
  for (int j = 0; j <= jtrans; j++) {
    foreach()
      if ((x - X0)/L0 < (double)j/jtrans)
	u.x[] = omgy[];
    squares ("u.x", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
    translate (y = 0.1)
      isoline2 ("cs", val = 0.5, np = {0,1,0}, alpha = py,
		lc = {98./256,78./256,44./256}, lw = 3);
    draw_string ("... shows the vertical vorticity field", 1, lw = 3);
    save ("endmov.mp4");
  }
  
  printf ("Add cells\n");
  jtrans = 20;
  for (int j = 0; j <= jtrans; j++) {
    squares ("omgy", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
    cells (n = {0,1,0}, alpha = py);
    translate (y = 0.1)
      isoline2 ("cs", val = 0.5, np = {0,1,0}, alpha = py,
		lc = {98./256,78./256,44./256}, lw = 3);
    draw_string ("... shows the vertical vorticity field", 1, lw = 3);
    save ("endmov.mp4");
  }
  
  printf ("Translate the slice downwards\n");
  pyst = pyen;
  pyen = Y0; 
  jtrans = 120;
  for (int j = 0; j <= jtrans; j++) {
    py = pyst + (double)j/jtrans*(pyen - pyst);
    squares ("omgy", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
    translate (y = 0.1)
      isoline2 ("cs", val = 0.5, np = {0,1,0}, alpha = py,
		lc = {98./256,78./256,44./256}, lw = 5);
    cells (n = {0,1,0}, alpha = py);
    draw_string ("Vorticity and cells", 1, lw = 3);
    save ("endmov.mp4");
  }

  printf ("Raise the slice again\n");
  pyst = pyen;
  pyen = 20; 
  jtrans = 150;
  for (int j = 0; j <= jtrans; j++) {
    py = pyst + (double)j/jtrans*(pyen - pyst);
    squares ("omgy", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
    translate (y = 0.1)
      isoline2 ("cs", val = 0.5, np = {0,1,0}, alpha = py,
		lc = {98./256,78./256,44./256}, lw = 5);
    cells (n = {0,1,0}, alpha = py);
    draw_string ("Vorticity and cells", 1, lw = 3);
    save ("endmov.mp4");
  }
  
  printf ("Some Last frames\n");
  jtrans = 30;
  squares ("omgy", n = {0,1,0}, alpha = py,
	   min = -1.1, max = 1.1, map = cool_warm);
  translate (y = 0.1)
    isoline2 ("cs", val = 0.5, np = {0,1,0}, alpha = py,
	      lc = {98./256,78./256,44./256}, lw = 3);
  cells (n = {0,1,0}, alpha = py);
  draw_string ("Vorticity and cells", 1, lw = 3);
  for (int j = 0; j <= jtrans; j++)
    save ("endmov.mp4");
}
