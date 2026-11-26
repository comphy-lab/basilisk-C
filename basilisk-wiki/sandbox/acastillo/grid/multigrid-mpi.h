#define MULTIGRID_MPI 1

static struct { void * var, * saved; } _mpi_sum_tab[64];
static int _mpi_sum_sp = 0;

static inline void _mpi_sum_push (void * v, MPI_Datatype datatype,
				  MPI_Op op, int elem)
{
  assert (_mpi_sum_sp < 64);
  _mpi_sum_tab[_mpi_sum_sp].var = v;
  _mpi_sum_tab[_mpi_sum_sp].saved = mpi_sum_reduce_init0 (v, datatype, op, elem);
  _mpi_sum_sp++;
}

static inline void _mpi_sum_pop (void * v, MPI_Datatype datatype,
				 MPI_Op op, int elem)
{
  for (int i = _mpi_sum_sp - 1; i >= 0; i--)
    if (_mpi_sum_tab[i].var == v) {
      void * saved = _mpi_sum_tab[i].saved;
      for (int j = i; j < _mpi_sum_sp - 1; j++)
	_mpi_sum_tab[j] = _mpi_sum_tab[j + 1];
      _mpi_sum_sp--;
      mpi_sum_reduce_array (saved, v, datatype, op, elem);
      return;
    }
  assert (0);
}

@undef mpi_sum_reduce_init
@define mpi_sum_reduce_init(s,v,datatype,op,elem) _mpi_sum_push(v,datatype,op,elem)
@undef mpi_sum_reduce_array
@define mpi_sum_reduce_array(s,v,datatype,op,elem) _mpi_sum_pop(v,datatype,op,elem)

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
		     MPI_Comm comm, MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  for (scalar s in list)
    size += s.block;
  Point point = { .level = level };
  SET_DIMENSIONS();
#if dimension == 2
  size *= (size_t)(point.n.y + 2*GHOSTS);
#elif dimension == 3
  size *= (size_t)(point.n.y + 2*GHOSTS)*(point.n.z + 2*GHOSTS);
#endif
  size *= GHOSTS*sizeof(real);
  real * buf = (real *) malloc (size), * b = buf;
  foreach_slice_x (i, i + GHOSTS, level)
    for (scalar s in list)
      for (scalar sb = s; sb.i < s.i + s.block; sb.i++, b++)
	memcpy (b, &sb[], sizeof(real));
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, comm, req);
  return buf;
}

foreach_dimension()
static void rcv_x (int i, int src, int tag, int level, scalar * list,
		   MPI_Comm comm)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  for (scalar s in list)
    size += s.block;
  Point point = { .level = level };
  SET_DIMENSIONS();
#if dimension == 2
  size *= (size_t)(point.n.y + 2*GHOSTS);
#elif dimension == 3
  size *= (size_t)(point.n.y + 2*GHOSTS)*(point.n.z + 2*GHOSTS);
#endif
  size *= GHOSTS*sizeof(real);
  real * buf = (real *) malloc (size), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, comm, &s);
  foreach_slice_x (i, i + GHOSTS, level)
    for (scalar s in list)
      for (scalar sb = s; sb.i < s.i + s.block; sb.i++, b++)
	memcpy (&sb[], b, sizeof(real));
  free (buf);
}

#if _GPU && _CUDA
# include "cuda/multigrid-mpi.h"
#endif // _GPU && _CUDA

trace
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{
  scalar * list1 = NULL;
#if _GPU && _CUDA
  scalar * list_gpu = NULL;
#endif
  for (scalar s in list)
    if (!is_constant(s) && s.block > 0) {
#if _GPU && _CUDA
      if (s.gpu.stored <= 0)
	list_gpu = list_add (list_gpu, s);
      else
#endif
	list1 = list_add (list1, s);
    }
#if _GPU && _CUDA
  if (list_gpu) {
    mpi_boundary_level_gpu (list_gpu, level);
    free (list_gpu);
  }
#endif
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
    Point point = { .level = level };
    SET_DIMENSIONS();
    int npl = point.n.x + 2*GHOSTS, nr = 0;
    if ((buf[0] = snd_x (npl - 2*GHOSTS, right, 0, level, list1,
			 mpi->cartcomm, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (GHOSTS, left, 1, level, list1,
			 mpi->cartcomm, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list1, mpi->cartcomm);
    rcv_x (npl - GHOSTS, right, 1, level, list1, mpi->cartcomm);
    MPI_Waitall (nr, reqs, MPI_STATUSES_IGNORE);
    free (buf[0]); free (buf[1]);
  }

  free (list1);

  prof_stop();
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
#if _GPU && _CUDA
  if (gpu_mpi_boundary == m)
    gpu_mpi_boundary = NULL;
#endif
  MPI_Comm_free (&m->cartcomm);
  free (m);
}

/**
## Process-grid factorization

Given a user aspect ratio `UserDims` $= (u_x, u_y, u_z)$ and
$k_{\text{cells}} = 2^{\text{depth}}$ cells per root, we partition the
domain across a process grid $p = (p_x, p_y, p_z)$ ($\prod p_d = \text{npe}$)
by refining the global root grid by a power-of-two factor $S = 2^s$ so
that each rank owns an integer number of roots
`LocalDims`$_d = (u_d S) / p_d$.

Divisibility $p_d \mid (u_d S)$ holds iff $r_d = p_d / \gcd(u_d, p_d)$
is a power of two and $S \ge \max_d r_d$, leaving a local multigrid
depth of $\log_2(k_{\text{cells}}/S)$. Among all valid factorizations,
we minimize $S$ first (maximizing local multigrid coarsening depth)
and then the per-rank halo surface area. */

static int mpi_gcd (int a, int b)
{
  while (b != 0) {
    int t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static bool eval_mpi_factorization (const int * p, int k_cells,
				    int * out_S, double * out_cost)
{
  const int * ud = (const int *) &UserDims;
  const int * per = (const int *) &Period;
  int S = 1;
  for (int d = 0; d < dimension; d++) {
    int r = p[d]/mpi_gcd (ud[d], p[d]);
    if ((r & (r - 1)) != 0)
      return false;
    if (r > S)
      S = r;
  }
  if (k_cells >= S && (k_cells % S) != 0)
    return false;
  int rem = k_cells < S ? 1 : (k_cells/S);
  if ((rem & (rem - 1)) != 0)
    return false;

  double ld[dimension];
  for (int d = 0; d < dimension; d++)
    ld[d] = (double)(ud[d]*S)/p[d];

  double cost = 0.;
  for (int d = 0; d < dimension; d++)
    if (p[d] > 1) {
      int nfaces = (per[d] || p[d] > 2) ? 2 : 1;
      double area = 1.;
      for (int j = 0; j < dimension; j++)
	if (j != d)
	  area *= ld[j];
      cost += nfaces*area;
    }

  *out_S = S;
  *out_cost = cost;
  return true;
}

static int choose_mpi_dims (int k_cells)
{
  int np = npe(), best_p[3] = {0, 0, 0}, best_S = 1 << 30;
  double best_cost = 1e300;

  if (UserMpiDims.x > 0) {
    int prod = 1;
    foreach_dimension()
      prod *= UserMpiDims.x;
    if (prod != np ||
	!eval_mpi_factorization ((const int *) &UserMpiDims, k_cells,
				 &best_S, &best_cost)) {
      if (pid() == 0)
	fprintf (stderr,
		 "%s:%d: error: requested mpi_dims() is incompatible"
		 " with npe=%d and N=%d\n", __FILE__, LINENO, np, N);
      exit (1);
    }
    MpiDims = UserMpiDims;
    return best_S;
  }

  for (int px = 1; px <= np; px++) {
    if (np % px) continue;
    int rem_x = np/px, py_max = dimension >= 2 ? rem_x : 1;
    for (int py = 1; py <= py_max; py++) {
      if (rem_x % py) continue;
      int pz = dimension >= 3 ? (rem_x/py) : 1;
      if (px*py*pz != np) continue;
      int cand_p[3] = {px, py, pz}, cand_S = 0;
      double cand_cost = 0.;
      if (eval_mpi_factorization (cand_p, k_cells, &cand_S, &cand_cost) &&
	  (cand_S < best_S ||
	   (cand_S == best_S && cand_cost < best_cost - 1e-12))) {
	best_S = cand_S;
	best_cost = cand_cost;
	memcpy (best_p, cand_p, sizeof(best_p));
      }
    }
  }
  if (best_p[0] == 0) {
    if (pid() == 0)
      fprintf (stderr,
	       "%s:%d: error: cannot factor npe=%d for current dimensions()"
	       " at N=%d\n", __FILE__, LINENO, np, N);
    exit (1);
  }
  for (int d = 0; d < dimension; d++)
    (&MpiDims.x)[d] = best_p[d];
  return best_S;
}

Boundary * mpi_boundary_new()
{
#if _GPU && _CUDA
  {
    MPI_Comm local_comm;
    MPI_Comm_split_type (MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
			 MPI_INFO_NULL, &local_comm);
    int local_rank = 0;
    MPI_Comm_rank (local_comm, &local_rank);
    MPI_Comm_free (&local_comm);
    cuda_set_device (local_rank);
  }
#endif
  MpiBoundary * m = qcalloc (1, MpiBoundary);
#if _GPU && _CUDA
  gpu_mpi_boundary = m;
#endif

  // choose process grid MpiDims and root-refinement factor S
  int k_cells = 1 << grid->depth;
  int S = choose_mpi_dims (k_cells);
  foreach_dimension() {
    Dimensions.x = UserDims.x*S;
    LocalDims.x = Dimensions.x/MpiDims.x;
  }

  // create Cartesian communicator and rank coordinates
  MPI_Cart_create (MPI_COMM_WORLD, dimension,
		   &MpiDims.x, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), dimension, mpi_coords);
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension()
    MpiOffset.x = mpi_coords[dir.x]*LocalDims.x;
#if _GPU && _CUDA
  int ld[3] = {1, 1, 1}, off[3] = {0, 0, 0};
  foreach_dimension() {
    ld[dir.x] = LocalDims.x;
    off[dir.x] = MpiOffset.x;
  }
  cuda_set_mpi_domain (ld[0], ld[1], ld[2], off[0], off[1], off[2]);
#endif

  // make sure other boundary conditions are not applied on internal faces
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
  int rem = k_cells < S ? 1 : (k_cells/S), r = 0;
  while (rem > 1)
    rem /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = Dimensions.x*(1 << r);
  grid->n = 1L << (dimension*depth());
  foreach_dimension()
    grid->n *= LocalDims.x;
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
#if _GPU
  for (scalar s in baseblock)
    if (s.gpu.stored < 0) {
      gpu_cpu_sync_scalar (s.i, s.block, multigrid->d,
			   multigrid->shift[depth() + 1], GPU_READ);
      s.gpu.stored = 0;
    }
  index.gpu.stored = 1;
#endif
  long cells_per_root = leaves ? (1L << (dimension*depth())) :
    ((1L << (dimension*(depth() + 1))) - 1)/((1L << dimension) - 1);
  long total_roots = 1;
  foreach_dimension()
    total_roots *= Dimensions.x;

  for (int ox = 0; ox < LocalDims.x; ox++)
#if dimension >= 2
    for (int oy = 0; oy < LocalDims.y; oy++)
#endif
#if dimension >= 3
      for (int oz = 0; oz < LocalDims.z; oz++)
#endif
      {
#if dimension == 1
	Point root = {GHOSTS + ox, 0};
#elif dimension == 2
	Point root = {GHOSTS + ox, GHOSTS + oy, 0};
#else
	Point root = {GHOSTS + ox, GHOSTS + oy, GHOSTS + oz, 0};
#endif
	long idx = global_root_index (root)*cells_per_root;
	foreach_cell_root (root) {
	  if (!leaves || is_leaf(cell))
	    index[] = idx++;
	  if (is_leaf(cell))
	    continue;
	}
      }
  boundary ({index});
  return pid() == 0 ? total_roots*cells_per_root - 1 : -1;
}
