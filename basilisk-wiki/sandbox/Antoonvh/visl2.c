/**
# Bwatch usage for visualization of a dump file

The dumpfile that contains a velocity field is generated on a super
computer with MPI. We render the data on a local machine, accelerated
with openMP.

![Volumetric rendering inspired by [Palas Kumar
 Farsoiya](https://gfm.aps.org/meetings/dfd-2020/5f53b6bb199e4c091e67bac0)](visl2/rings.png)
 
 ![The same rendering with voricity-vector-direction color coding
 appears to be a bit over the top](visl2/color_rings.png)

 ![The same rendering with a grid](visl2/rings_cells.png)


 ![Visualization of the $\lambda_2$ field as volume, and later as iso
 surface, color coded with the local velocity vector](visl2/mov.mp4)(width = "600")
*/

#include "grid/octree.h"
#include "utils.h"
#include "bwatch.h"
#include "lambda2.h"

scalar dummy[];
vector uc[];

void vorticity3D (vector u, vector omg) {
  foreach() {
    foreach_dimension() {
      omg.z[] += (u.y[1] - u.y[-1] -
		  u.x[0,1] + u.x[0,-1])/(2*Delta);
    if (isnan(omg.z[])) //fixme
      omg.z[] = 0.;
    }
  }
  boundary ((scalar*){omg});
 
}

void swap_to_cylinder (vector a) {
  foreach() {
    double phi = fabs(x) > 0. ? atan(y/x) : pi/2.;
    coord at = {a.x[], a.y[], a.z[]};
    a.x[] = cos(phi)*at.x + sin(phi)*a.y[];
    a.y[] = cos(phi)*at.y - sin(phi)*a.x[];
   
  }
  boundary ({a.x, a.y});
}

int main() {
  //Obtain the dump file from elsewhere
  system ("test -f download || wget https://surfdrive.surf.nl/files/index.php/s/2UgDG4iEjjZVaiL/download"); 
  restore ("download");
  scalar l2[];
  vector omega[];
  lambda2 (uc, l2);
  vorticity3D (uc, omega);
  swap_to_cylinder (omega);
  foreach()
    l2[] = l2[] > 0 ? 0 : -l2[];
  boundary ({l2});
  
  watch (fov = 50, nx = 2400, ny = 2400); //Highres
  watch (O = {30, 30., 100}, poi = {10,10,0}, fov = L0/4);
  
  volume (l2, cols = true, sc = 0.005, mval = 0.001,
	  min = -.1, max = .1, shading = 1);
  lattice(width = 0.02);
  store (fopen("smkcells.ppm", "w"));
  plain();

  watch (O = {30, 30., 100}, poi = {0.01,0.01,0.01}, fov = 50);
  
  volume (l2, cols = true, sc = 0.005, mval = 0.001,
	  colorv = omega, max = 10, shading = 1);
  store (fopen("smkc.ppm", "w"));
  plain();
  
  volume (l2, cols = true, sc = 0.005, mval = 0.001,
	  min = -.1, max = .1, shading = 1);
  store (fopen("smk.ppm", "w"));
  plain();
  
  system( "convert smkc.ppm -resize 50% color_rings.png"); //4x MSAA
  system( "convert smk.ppm -resize 50% rings.png"); //4x MSAA
  system( "convert smkcells.ppm -resize 50% rings_cells.png"); //4x MSAA
#if 0
  watch (nx = 900, ny = 900);
  FILE * fp = popen ("ppm2mp4 -r 20 mov.mp4", "w");
  for (double a = 0; a <= 2*pi; a += 0.025*pi) {
    watch (O = {100*sin(a), 30., 100*cos(a)});
    volume (l2, sc = 0.01, mval = 0.001, col = {255, 20, 30});
    sphere (R = 200, mat = {dull = true, col = {255, 255, 255}});
    store (fp);
    plain();
  }
  for (double a = 0; a <= 2*pi; a += 0.025*pi) {
    watch (O = {100*sin(a), 30., 100*cos(a)});
    equiplane (l2, 0.005, insideout = true,
	       mat = {v = uc, min = -0.05, max = 0.05});
    sphere (R = 200, mat = {dull = true, col = {255, 255, 255}});
    store (fp);
    plain();
  }
  
  pclose(fp);
#endif
} 
