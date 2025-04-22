/**
# Computation of volume fractions on a variable-resolution grid

This test case is similar to [fractions.c]() but with a refined band
in the middle. */

#include "fractions.h"

/**
The refinement band is defined by this function. */

int main()
{
  origin (-0.5, -0.5);
  init_grid (16);

  refine (level == 4 && fabs (x) < 0.25);

  /**
  We use a circle of radius 0.3 and initialise the fractions. */

  scalar c[];
  c.refine = c.prolongation = fraction_refine;
  face vector s[];
  solid (c, s, sq(0.3) - sq(x) - sq(y));

  /**
  Output the reconstruced facets and cells. */

  output_facets (c, stdout, s);
  output_facets (c, stderr);
  FILE * fp = fopen ("cells", "w");
  output_cells (fp);
  fclose (fp);

  /**
  Finally, we reconstruct the interface and display the reconstructed
  facets only in the "halo cells". This is a check of the consistency of
  the boundary conditions applied to $\mathbf{n}$ and $\alpha$. */

  vector n[];
  scalar alpha[];
  reconstruction (c, n, alpha);
  boundary ({n, alpha});
  coord p[2];
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l) 
      foreach_child() {
        coord m = {n.x[],n.y[]};
	if (facets (m, alpha[], p) == 2)
	  fprintf (stderr, "halo %g %g\nhalo %g %g\nhalo\n", 
		   x + p[0].x*Delta, y + p[0].y*Delta, 
		   x + p[1].x*Delta, y + p[1].y*Delta);
      }
}

/**
This gives this figure where "exact" uses *c* and *s*, "VOF" uses
only *c* and "halo" is the halo cell reconstuction.

~~~gnuplot Exact and VOF-reconstucted interface
set size ratio -1
set key out
plot 'cells' w l t '', 'out' w l t "exact", '< grep -v halo log' w l t "VOF", \
     '< grep halo log' u 2:3 w p t 'halo'
~~~
*/
