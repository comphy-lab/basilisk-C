/**
# Tagging connected neighborhoods

The goal is to associate a unique (strictly positive) index ("tag")
*t* to cells which belong to the same "neighborhood". All cells in a
neighborhood have an initial tag value (given by the user) which is
non-zero and are separated from cells in other neighborhoods by cells
which have an initial (and final) tag value of zero.

We first define the restriction function for tag values. The parent
tag is just the minimum of its children's tags. */

static void restriction_tag (Point point, scalar t)
{
  double min = HUGE;
  foreach_child()
    if (t[] < min)
      min = t[];
  t[] = min;
}

/**
We also need a few helper functions. The function below implements a
[binary search](https://en.wikipedia.org/wiki/Binary_search_algorithm)
of a sorted array. It returns the index in the array so that $a[s] \leq
tag < a[s+1]$. */

static long lookup_tag (Array * a, double tag)
{
  long len = a->len/sizeof(double);
  double * p = (double *) a->p;
  if (tag == p[0])
    return 0;
  if (tag == p[len - 1])
    return len - 1;

  long s = 0, e = len - 1;
  while (s < e - 1) {
    long m = (s + e)/2;
    if (p[m] <= tag)
      s = m;
    else
      e = m;
  }
  return s;
}

#if _MPI
static int compar_double (const void * p1, const void * p2)
{
  const double * a = p1, * b = p2;
  return *a > *b;
}
#endif

/**
The function just takes the scalar field *t* which holds the initial
and final tag values. It returns the maximum neighborhood tag value
(which is also the number of neighborhoods). */

trace
int tag (scalar t)
{

  /**
  We first set the restriction and prolongation functions (on trees). */

  t.restriction = restriction_tag;
#if TREE  
  t.refine = t.prolongation = refine_injection;
  t.dirty = true;
#endif

  /**
  As an initial guess, we set all the (leaf) cells which have a
  non-zero initial tag value to the Z- (or Morton-) index. We thus
  have a different "neighborhood" for each cell which has a non-zero
  initial tag, and a single neighborhood (tagged zero) for all the
  cells which have a zero initial tag. */
  
#if _MPI
  scalar index[];
  z_indexing (index, true);
  foreach()
    t[] = (t[] != 0)*(index[] + 1);
#else // !_MPI
  long i = 1;
  foreach_cell()
    if (is_leaf(cell)) {
      t[] = (t[] != 0)*i++;
      continue;
    }
#endif // !_MPI

  /**
  To gather cells which belong to the same neighborhood, we repeat
  multigrid iterations until none of the tag values changes. */
  
  bool changed;
  do {

    /**
    We first do a restriction from the finest to the coarsest level of
    the multigrid hierarchy, using the `restriction_tag()` function
    above. */
    
    restriction ({t});

    /**
    We then go from the coarsest to the finest level and update the
    tag values. */

    changed = false;
    for (int l = 1; l <= grid->maxdepth; l++) {

      /**
      If the parent tag is non-zero, we set the child tag to the value
      of its parent (i.e. to the minimum tag value of its
      siblings). */
      
      foreach_level(l)
	if (coarse(t))
	  t[] = coarse(t);
      boundary_level ({t}, l);

      /**
      For cells which verify the threshold condition (i.e. for which
      `t[] != 0`), we refine this initial guess by taking the minimum
      (non-zero) tag value of its closest neighbors. We also track
      whether this update changes any of the tag values. */
      
      foreach_level (l, reduction(||:changed))
        if (t[]) {
	  double min = t[];
	  foreach_neighbor(1)
	    if (t[] && t[] < min)
	      min = t[];

	  /**
	  On trees, we need to take into account the minimum tag value
	  of neighboring fine cells. */
	  
#if TREE
	  foreach_dimension()
	    for (int i = -1; i <= 2; i += 3)
	      if (is_refined (neighbor((2*i - 1)/3)))
		for (int j = 0; j <= 1; j++)
		  for (int k = 0; k <= 1; k++)
		    if (fine(t,i,j,k) && fine(t,i,j,k) < min)
		      min = fine(t,i,j,k);
#endif // TREE

	  if (t[] != min) {
	    changed = true;
	    t[] = min;
	  }
	}
      boundary_level ({t}, l);
    }
  } while (changed);

  /**
  ## Reducing the range of indices

  Each neighborhood is now tagged with a unique index. The range of
  indices is large however (between one and the total number of
  leaves). The goal of this step is to reduce this range to between
  one and the number of neighborhoods. To do this, we create an ordered
  array of unique indices. */

  Array * a = array_new();
  foreach (serial)
    if (t[] > 0) {
  
      /**
      We first check whether the index is larger than the maximum or
      smaller than the minimum value in the array. *s* is the
      position of the (possibly new) index in the array. A negative
      value means that the index is already in the array.  */
  
      double tag = t[], * ap = (double *) a->p;
      long s = -1;
      if (a->len == 0 || tag > ap[a->len/sizeof(double) - 1])
	s = a->len/sizeof(double);
      else if (tag < ap[0])
	s = 0;
      else {
	
	/**
	We find the range of existing indices [s-1:s] which contains
	the index. We check whether the index is already in the
	array. */
	
	s = lookup_tag (a, tag) + 1;
	if (tag == ap[s - 1] || tag == ap[s])
	  s = -1;
      }
      if (s >= 0) {

	/**
	If the index is new, we add it to the array in the correct
	position (s). */
	
	array_append (a, &tag, sizeof(double)), ap = (double *) a->p;
	for (int i = a->len/sizeof(double) - 1; i > s; i--)
	  ap[i] = ap[i-1];
	ap[s] = tag;
      }
    }

  /**
  ## Parallel reduction

  Each process now has its own local correspondence map between the
  neighborhood index and its rank in the array *a* (i.e. its reduced
  index). In parallel, we now need to build a global correspondence map. 

  We first get the maximum size over all processes of the local map
  and increase the size of the local map to this value. */

#if _MPI
  long lmax = a->len;
  mpi_all_reduce (lmax, MPI_LONG, MPI_MAX);
  a->p = realloc (a->p, lmax);
  lmax /= sizeof(double);

  /**
  We then gather all the local maps into a global map and sort it. All
  local arrays need to be of the same size to be able to use
  MPI_Allgather(), so we first pad the arrays with -1. */
  
  double * q = a->p;
  for (int i = a->len/sizeof(double); i < lmax; i++)
    q[i] = -1;
  double p[lmax*npe()];
  MPI_Allgather (a->p, lmax, MPI_DOUBLE, p, lmax, MPI_DOUBLE, MPI_COMM_WORLD);
  qsort (p, lmax*npe(), sizeof(double), compar_double);

  /**
  This sorted global map will (probably) contain duplicated entries
  (i.e. indices of neighborhoods which span multiple processes). To
  build the new global map (stored in *a*), we eliminate these, as well
  as the negative indices which were used to pad the local arrays. */
  
  array_free (a);
  a = array_new();
  double last = -1;
  for (int i = 0; i < lmax*npe(); i++)
    if (p[i] != last) {
      array_append (a, &p[i], sizeof(double));
      last = p[i];
    }
#endif

  /**
  Once we have the (global) map, we can replace the neighborhood
  indices with their index in the global map (+1). */
  
  foreach()
    if (t[] > 0)
      t[] = lookup_tag (a, t[]) + 1;

  /**
  We return the maximum index value. */
  
  int n = a->len/sizeof(double);
  array_free (a);
  return n;
}

/**
# Removing (small) droplets/bubbles

Using tag(), the function below can identify and remove droplets (or
bubbles) defined by VOF tracer *f* (resp. $1 - f$), smaller than a
given diameter (*minsize*) expressed in number of cells. */

void remove_droplets (scalar f,
		      int minsize = 3,
		      double threshold = 1e-4,
		      bool bubbles = false)
{
  scalar d[];
  foreach()
    d[] = (bubbles ? 1. - f[] : f[]) > threshold;
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++)
    size[i] = 0;
  foreach (serial)
    if (d[] > 0)
      size[((int) d[]) - 1]++;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  minsize = pow (minsize, dimension);
  foreach()
    if (d[] > 0 && size[((int) d[]) - 1] < minsize)
      f[] = bubbles;
}
