/**
# Distance field computation from a 3D model

The goal is to build the skeleton of a 3D
[CAD](https://en.wikipedia.org/wiki/Computer-aided_design) model. */

#include "grid/octree.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "../thinning.h"
int n_part = 0;
#include "../../Antoonvh/scatter.h"
#include "../../aberny/output_stl.h"



int main()
{
  system ("test -f Jaekelopterus.stl || "
  "wget \"https://drive.google.com/uc?export=download&id=1x3b3oc-GFDw4VWd4n3i2YFg5bLwlFdWO\" -O Jaekelopterus.stl");


  coord * p = input_stl (fopen ("Jaekelopterus.stl", "r"));
  // taken from https://www.thingiverse.com/thing:4702654
  coord min, max;
  bounding_box (p, &min, &max);  
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  
  init_grid (8);
  size (1.2*maxl);
  origin ((max.x + min.x)/2. - L0/2,
    (max.y + min.y)/2. - L0/2,
    (max.z + min.z)/2. - L0/2);

  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){3.e-4*L0}, 11).nf);

  /**
  We also compute the volume and surface fractions from the distance
  field. We first construct a vertex field interpolated from the
  centered field and then call the appropriate VOF functions. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
       d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  scalar f[];
  face vector s[];
  fractions (phi, f, s);
  
  // clear();
  view (fov = 17.1577, quat = {0.509824,-0.191867,-0.256451,0.798436}, tx = -0.0771716, ty = -0.131169, 
    width = 600, height = 600, bg = {1,1,1});
  /**
  We display an isosurface of the distance function coloured with the
  level of refinement and the surface reconstructed from volume fractions. */

  isosurface ("d", 0, color = "level", min = 5, max = 11);
  save ("Jaekelopterus.png");
  draw_vof ("f", "s", edges = true, lw = 0.5);
  save ("vof.png");

  view (fov = 8.76256, quat = {0.297217,-0.150457,-0.437522,0.835224}, tx =
    -0.0413951, ty = 0.0647216);
  draw_vof ("f", "s", edges = true, lw = 0.5);
  save ("vof2.png");
  isosurface ("d", 0, color = "level", min = 5, max = 11);
  save ("Jaekelopterus2.png");

  scalar c[];
  foreach(){
    if(f[])c[] =1;
    else c[] = 0;
  }
  boundary({c});
  thinning3D(c);

/**
To output the skeleton, we will use the Lagrangian particles of Antoonvh, the
remaining cells will be shown as spheres.
*/

  Cache skeleton = {0}; 

  foreach(){
    if(c[]){
      cache_append (&skeleton, point, 0);
      n_part++;
    }
  }

  coord * loc = malloc (n_part*sizeof(coord));
  int n = 0;
  foreach_cache(skeleton) {
    coord cc = {x, y, z};
    foreach_dimension()
      loc[n].x = cc.x;
    n++;
  }
  view (fov = 17.1577, quat = {0.509824,-0.191867,-0.256451,0.798436}, tx = -0.0771716, ty = -0.131169, 
    width = 600, height = 600, bg = {1,1,1});
  scatter(loc, s = 40);
  save("medial_axis.png");
  view (fov = 8.76256, quat = {0.297217,-0.150457,-0.437522,0.835224}, tx =
    -0.0413951, ty = 0.0647216);
  scatter(loc, s = 40);
  save("medial_axis2.png"); 

  free (loc);
  free (skeleton.p);
  stl_output_binary(f, "f_field.stl");
}

/**
We correctly detect the spikes on the claws of our creature and so it should be
able to eat well, the thinning algorithm seem to behave in these parts.

However, the thinning operation is more noisy on the shell of the
creature, maybe a topological analysis with persistance diagram would give
better results in these regions.

|  Isosurface of the distance function coloured with level of refinement   | VOF   |      Skeleton      | 
|:-------------:|:-------------:|:-------------:|
| ![](distance/Jaekelopterus.png) | ![](distance/vof.png) |![](distance/medial_axis.png)
| ![](distance/Jaekelopterus2.png) | ![](distance/vof2.png) | ![](distance/medial_axis2.png)

## See also

* [Computation of a levelset field from a contour](/src/test/basilisk.c)
*/
