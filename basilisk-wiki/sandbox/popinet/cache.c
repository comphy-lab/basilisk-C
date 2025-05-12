/**
# Traversing a subset of the mesh

The goal is to build an array/list/cache of cells which can then be
used to traverse only a subset of the entire mesh. */

int main()
{

  /**
  We first setup the mesh and shift the origin. */

  init_grid (32);
  origin (-0.5, -0.5);

  /**
  We will use the pre-defined structure *Cache* and associated
  functions. Note that for the moment these are defined in
  [/src/grid/tree.h]() i.e. this will only work on
  bi/quad/octrees. All fields of *a* are initially zero. */
  
  Cache a = {0};

  /**
  We then traverse the entire mesh and add the cells which are inside
  a circle to the cache. *point* is the index of the current cell. The
  last argument is an optional integer associated with the cell. */
  
  foreach()
    if (x*x + y*y < sq(0.25))
      cache_append (&a, point, 0);

  /**
  We "optimize" the size of the cache i.e. shrink the dynamic array of
  indices if necessary. */
  
  cache_shrink (&a);

  /**
  We can then simply traverse the cells contained in cache *a* like
  this: */
  
  foreach_cache (a)
    fprintf (stderr, "%g %g\n", x, y);

  /**
  We output the mesh to make the figure below. */
  
  output_cells (stdout);

  /**
  Finally, we free the cache. */
  
  free (a.p);
}

/**
~~~gnuplot Cells and cached cells
set size ratio -1
unset key
plot 'out' w l, 'log'
~~~
*/
