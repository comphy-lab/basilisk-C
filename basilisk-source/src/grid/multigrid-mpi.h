#define MULTIGRID_MPI 1

#if dimension == 1

macro2 foreach_slice_x (int start, int end, int l) {
  {
    int ig = 0; NOT_UNUSED(ig);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = start; point.i < end; point.i++)
      {...}
  }
}

#elif dimension == 2

macro2 foreach_slice_x (int start, int end, int l) {
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = start; point.i < end; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*GHOSTS; point.j++)
	{...}
  }
}

macro2 foreach_slice_y (int start, int end, int l) {
  {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = 0; point.i < point.n.x + 2*GHOSTS; point.i++)
      for (point.j = start; point.j < end; point.j++)
	{...}
  }
}

#elif dimension == 3

macro2 foreach_slice_x (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = start; point.i < end; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*GHOSTS; point.j++)
	for (point.k = 0; point.k < point.n.z + 2*GHOSTS; point.k++)
	  {...}
  }
}

macro2 foreach_slice_y (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = 0; point.i < point.n.x + 2*GHOSTS; point.i++)
      for (point.j = start; point.j < end; point.j++)
	for (point.k = 0; point.k < point.n.z + 2*GHOSTS; point.k++)
	  {...}
  }
}

macro2 foreach_slice_z (int start, int end, int l) {
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l; SET_DIMENSIONS();
    for (point.i = 0; point.i < point.n.x + 2*GHOSTS; point.i++)
      for (point.j = 0; point.j < point.n.y + 2*GHOSTS; point.j++)
	for (point.k = start; point.k < end; point.k++)
	  {...}
  }
}

#endif // dimension == 3

typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;

foreach_dimension()
static void * snd_x (int i, int dst, int tag, int level, scalar * list,
		     MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  for (scalar s in list)
    size += s.block;
  size *= pow((1 << level) + 2*GHOSTS, dimension - 1)*GHOSTS*sizeof(real);
  double * buf = (double *) malloc (size), * b = buf;
  foreach_slice_x (i, i + GHOSTS, level)
    for (scalar s in list)
      for (scalar sb = s; sb.i < s.i + s.block; sb.i++, b++)
	memcpy (b, &sb[], sizeof(real));
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}

foreach_dimension()
static void rcv_x (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  for (scalar s in list)
    size += s.block;
  size *= pow((1 << level) + 2*GHOSTS, dimension - 1)*GHOSTS*sizeof(real);
  double * buf = (double *) malloc (size), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
  foreach_slice_x (i, i + GHOSTS, level)
    for (scalar s in list)
      for (scalar sb = s; sb.i < s.i + s.block; sb.i++, b++)
	memcpy (&sb[], b, sizeof(real));
  free (buf);
}

trace
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  prof_start ("mpi_boundary_level");
  
  if (level < 0) level = depth();  
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension() {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);  
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*GHOSTS, nr = 0;
    if ((buf[0] = snd_x (npl - 2*GHOSTS, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left,  1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left,  0, level, list1);
    rcv_x (npl - GHOSTS,   right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    free (buf[0]); free (buf[1]);
  }

  free (list1);

  prof_stop();
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  free (m);
}

static void mpi_dimensions_error (int n)
{
  fprintf (stderr,
	   "%s:%d: error: the number of MPI processes must be equal to ",
	   __FILE__, LINENO);
  if (n > 1)
    fprintf (stderr, "%dx", n);
  fprintf (stderr, "%d^i\n", 1 << dimension);
  exit (1);  
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = qcalloc (1, MpiBoundary);
  int n = 1;
  foreach_dimension()
    n *= Dimensions.x;
  if (npe() % n)
    mpi_dimensions_error (n);
  int j = npe()/n, i = 0;
  while (j > 1) {
    if (j % (1 << dimension))
      mpi_dimensions_error (n);
    j /= 1 << dimension;
    i++;
  }
  foreach_dimension()
    Dimensions.x *= 1 << i;
  MPI_Dims_create (npe(), dimension, &Dimensions.x);
  MPI_Cart_create (MPI_COMM_WORLD, dimension,
		   &Dimensions.x, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), dimension, mpi_coords);

  // make sure other boundary conditions are not applied
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension() {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  }

  // rescale the resolution
  Dimensions_scale = Dimensions.x;
  N /= Dimensions.x;
  int r = 0;
  while (N > 1)
    N /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = Dimensions.x*(1 << r);
  grid->n = 1 << dimension*depth();
  grid->tn = npe()*grid->n;
  
  // setup boundary methods and add to list of boundary conditions
  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}

trace
double z_indexing (scalar index, bool leaves)
{
  long i;
  if (leaves)
    i = pid()*(1 << dimension*depth());
  else
    i = pid()*((1 << dimension*(depth() + 1)) - 1)/((1 << dimension) - 1);
  foreach_cell() {
    if (!leaves || is_leaf(cell))
      index[] = i++;
    if (is_leaf(cell))
      continue;
  }
  boundary ({index});
  return pid() == 0 ? i*npe() - 1 : -1;
}
