/**
# Loading a dump file for visualization purposes */

/**
The first step is to load the relevant header files (*.h*) used to
generate the dump file. For example, I load for embedded
boundaries: */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"

/**
Then include [/src/view.h]() along with any other data visualization
header files ([/src/lambda2.h](), ...). */

#include "view.h"
#include "lambda2.h"

/**
The *main()* is where the magic happens. Indeed, since we are not
running a simulation per se, we cannot use the event structures of
Basilisk. */

int main ()
{
  /**
  The size of the domain and the origin must match those used to
  generate the dump file. */

  L0 = 32.;
  size (L0);
  origin (-L0/2., 0., -L0/2.);

  /**
  We now [restore](/src/output.h#restore) the dump file.

  Note that when using the [centered Navier-Stokes
  solver](/src/navier-stokes/centered.h), the pressure is by default
  not dumped. To change this, simply set *p.nodump = false;* in the
  **.c* file used to generate the dump file. */

  char name[80];
  sprintf (name, "mydump");
  bool restart = restore (name);
  assert (restart == true);
	
  /**
  We now perform all the visualization we want. */
  
  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});

  scalar omega[];
  vorticity (u, omega); // Vorticity in xy plane
  boundary ({omega});

  clear();
  view (fov = 20, camera = "front",
	tx = -(0.)/(L0), ty = -(0.)/(L0),
	bg = {1,1,1},
	width = 800, height = 800);

  cells (n = {0., 0., 1.}, alpha = 1.e-12);
  stats somega = statsf (omega);
  squares ("omega",
	   n = {0., 0., 1.}, alpha = 1.e-12, // To avoid roundoff errors.
	   min  = somega.min, max = somega.max,
	   map = cool_warm);
  save ("omega.png")

  /**
  We finaly free the grid the avoid memory leaks. */
  
  free_grid();
  return 0;
}
