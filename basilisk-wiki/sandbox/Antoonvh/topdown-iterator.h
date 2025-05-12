/**
We write grid-cell iterators the go from top to bottom, rather than
the default bottom to top. This could be useful in exotic scenarios. 

Normal N-order for trees:

2   4
| \ |
1   3

new mirrored-N order:

1   3
| / |
2   4
 */
#include "foreach_cell_usd.h" 

/**
Now we can create usd versions of `foreach_level(int l)` and
`foreach_level_or_leaf()` and `foreach_usd()`, which concludes the
excersize.
 */

Grid * grid2 = NULL;
#define tree2 ((Tree *) grid2)

@def foreach_level2(l) {
  if (l <= depth()) {
    update_cache_f2();
    CacheLevel _active = tree2->active[l];
    foreach_cache_level (_active,l)
@
@define end_foreach_level2() end_foreach_cache_level(); }}

@define foreach_coarse_level2(l) foreach_level2(l) if (!is_leaf(cell)) {
@define end_foreach_coarse_level2() } end_foreach_level2()

@def foreach_level_or_leaf2(l) {
  for (int _l1 = l; _l1 >= 0; _l1--)
    foreach_level2(_l1)
      if (_l1 == l || is_leaf (cell)) {
@
@define end_foreach_level_or_leaf2() } end_foreach_level2(); }

@define foreach_usd() update_cache_f2(); foreach_cache(tree2->leaves)
@define end_foreach_usd()  end_foreach_cache()


/**
A function that Updates the caches of `grid2`.  
 */
trace
static void update_cache_f2 (void)
{
  grid2 = grid;
  Tree * q = tree2;
  check_periodic (q);
  
  foreach_cache (q->vertices)
    if (level <= depth() && allocated(0))
      cell.flags &= ~vertex;
  
  /* empty caches */
  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;
  
#if FBOUNDARY
  const unsigned short fboundary = 1 << user;
  foreach_cell2() {
#else    
  foreach_cell_all2() {
#endif
    
    if (is_local(cell) && is_active(cell)) {
      // active cells
      //      assert (is_active(cell));
      cache_level_append (&q->active[level], point);
    }
#if !FBOUNDARY
    if (is_boundary(cell)) {
      // boundary conditions
      bool has_neighbors = false;
      foreach_neighbor (BGHOSTS)
	if (allocated(0) && !is_boundary(cell))
	  has_neighbors = true, break;
      if (has_neighbors)
	cache_level_append (&q->boundary[level], point);
      // restriction for masked cells
      if (level > 0 && is_local(aparent(0)))
	cache_level_append (&q->restriction[level], point);
    }
#else
    // boundaries
    if (!is_boundary(cell)) {
      // look in a 5x5 neighborhood for boundary cells
      foreach_neighbor (BGHOSTS)
	if (allocated(0) && is_boundary(cell) && !(cell.flags & fboundary)) {
	  cache_level_append (&q->boundary[level], point);
	  cell.flags |= fboundary;
	}
    }
    // restriction for masked cells
    else if (level > 0 && is_local(aparent(0)))
      cache_level_append (&q->restriction[level], point);
#endif
    if (is_leaf (cell)) {
      if (is_local(cell)) {
	cache_append (&q->leaves, point, 0);
	// faces
	unsigned short flags = 0;
	foreach_dimension()
	  if (is_boundary(neighbor(-1)) || is_prolongation(neighbor(-1)) ||
	      is_leaf(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, flags);
	foreach_dimension()
	  if (is_boundary(neighbor(1)) || is_prolongation(neighbor(1)) ||
	      (!is_local(neighbor(1)) && is_leaf(neighbor(1))))
	    cache_append (&q->faces, neighborp(1), face_x);
	// vertices
	for (int i = 0; i <= 1; i++)
        #if dimension >= 2
	  for (int j = 0; j <= 1; j++)
        #endif
          #if dimension >= 3
	    for (int k = 0; k <= 1; k++)
	  #endif
	      if (!is_vertex(neighbor(i,j,k))) {
		cache_append (&q->vertices, neighborp(i,j,k), 0);
		neighbor(i,j,k).flags |= vertex;
	      }
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
      }
      else if (!is_boundary(cell) || is_local(aparent(0))) { // non-local
	// faces
	unsigned short flags = 0;
	foreach_dimension()
	  if (allocated(-1) &&
	      is_local(neighbor(-1)) && is_prolongation(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append_face (point, flags);
	foreach_dimension()
	  if (allocated(1) && is_local(neighbor(1)) &&
	      is_prolongation(neighbor(1)))
	    cache_append_face (neighborp(1), face_x);
      }
#if FBOUNDARY // fixme: this should always be included
      continue; 
#endif
    }
  }

  /* optimize caches */
  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}
  
  q->dirty = false;

#if FBOUNDARY
  for (int l = depth(); l >= 0; l--)
    foreach_boundary_level (l)
      cell.flags &= ~fboundary;
#endif
  
  // mesh size
  grid->n = q->leaves.n;
  // for MPI the reduction operation over all processes is done by balance()
@if !_MPI
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
@endif
  
}

  
