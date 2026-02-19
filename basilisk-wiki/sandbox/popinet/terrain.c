/**
# How to construct a 3D terrain from a KDT database

The goal is to initialize a 3D surface (for example an embedded solid
boundary) using a KDT database, similarly to what is done in two
dimensions in the [tsunami
example](/src/examples/tsunami.c#solver-setup).

We will use a 3D multigrid, the terrain module to read the KDT
database and Basilisk view for visualisation. */

#include "grid/multigrid3D.h"
#include "terrain.h"
#include "view.h"

/**
This is an example of a function, which maybe tuned according to
coordinate systems etc. */

void distance_from_terrain (scalar d, const char * name, double ox, double oy, double width)
{

  /**
  We scale the domain according to the piece of data we want to select
  in the database. */
  
  origin (ox, oy);
  size (width);

  /**
  We call the terrain module to read the data in a (temporary)
  field. Note that since the KDT database is only 2D, all (3D) points
  having the same (x,y) coordinates will have identical values of
  zb. */
  
  scalar zb[];
  terrain (zb, name, NULL);

  /**
  We output an image, on the $z = 0$ plane, just to check that this
  went OK. */
  
  output_ppm (zb, file = "zb.png", spread = -1, n = 256);  

  /**
  ![This indeed looks like India, the Himalayas etc.](terrain/zb.png)

  We rescale vertically to compute the distance of each 3D point: this
  will need to be adapted according to your choices of vertical
  coordinates. */
  
  stats s = statsf (zb);
  size (s.max - s.min); // this needs to be tuned according to vertical/horizontal scaling
  origin (z = s.min);

  /**
  The "distance" (not Euclidean) is simply the difference between the
  height of the topography and the local z coordinate. */
  
  foreach()
    d[] = zb[] - z;  
}

int main()
{

  /**
  We use an uniform grid. Note that you may want to use an "adaptive
  initial refinement" to minimize the startup cost when you are
  pushing the resolution. */
  
  init_grid (64);

  /**
  We call our function with the appropriate choice of coordinates:
  here the same piece of land as in the [tsunami
  example](/src/examples/tsunami.c). */
  
  scalar d[];
  distance_from_terrain (d, "~/terrain/etopo2", 94. - 54./2., 8. - 54/2., 54.);

  /**
  We can rescale again as desired. */
  
  size (1.);
  origin (0, 0, 0);

  /**
  We display an isosurface (at zero) of our distance function, colored with z. */
  
  view (quat = {0.463, 0.107, 0.198, 0.857},
	fov = 30, near = 0.01, far = 1000,
	tx = -0.589, ty = -0.483, tz = -3.172,
	width = 1037, height = 936, bg = {1,1,1});
  box();
  isosurface ("d", 0, color = "z", spread = -1, linear = true);
  save ("iso.png");

  /**
  ![Isosurface of the distance function](terrain/iso.png){width="80%"}

  We can also use this distance function to compute the volume and
  surface fractions required for embedded solid boundaries. */
  
  scalar cs[];
  face vector fs[];
  solid (cs, fs, (d[] + d[-1] + d[0,-1] + d[-1,-1] +
		  d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.);

  /**
  And check that this is indeed OK. */
  
  box();
  draw_vof ("cs", "fs");
  draw_vof ("cs", "fs", edges = true, lw = 0.5);
  save ("solid.png");

  /**
  ![Reconstructed solid surface](terrain/solid.png){width="80%"}
  
  And finally we dump the result in case we want to [visualize with
  bview](/src/bview.c). */
  
  dump();
}

/**
## See also

* [Distance field computation from a 3D model](/src/examples/distance.c)
* [xyz2kdt examples](xyz2kdt.md)
*/
