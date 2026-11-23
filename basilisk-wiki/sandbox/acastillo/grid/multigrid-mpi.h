#define MULTIGRID_MPI 1

static struct { void * var, * saved; } _mpi_sum_tab[64];
static int _mpi_sum_sp = 0;

static inline void _mpi_sum_push (void * v, MPI_Datatype datatype, MPI_Op op, int elem)
{
  assert (_mpi_sum_sp < 64);
  _mpi_sum_tab[_mpi_sum_sp].var = v;
  _mpi_sum_tab[_mpi_sum_sp].saved = mpi_sum_reduce_init0 (v, datatype, op, elem);
  _mpi_sum_sp++;
}

static inline void _mpi_sum_pop (void * v, MPI_Datatype datatype, MPI_Op op, int elem)
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

static MpiBoundary * current_mpi_boundary = NULL;

foreach_dimension()
static void * snd_x (int i, int dst, int tag, int level, scalar * list,
		     MPI_Request * req)
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
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, current_mpi_boundary->cartcomm, req);
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
  MPI_Recv (buf, size, MPI_BYTE, src, tag, current_mpi_boundary->cartcomm, &s);
  foreach_slice_x (i, i + GHOSTS, level)
    for (scalar s in list)
      for (scalar sb = s; sb.i < s.i + s.block; sb.i++, b++)
	memcpy (&sb[], b, sizeof(real));
  free (buf);
}

#if _GPU && _CUDA
static int use_cuda_aware_mpi (void)
{
  static int cached = -1;
  if (cached < 0)
    cached = cuda_query_mpi_support ();
  return cached;
}

trace
void mpi_boundary_level_gpu (scalar * list, int level)
{
  if (!current_mpi_boundary || !list)
    return;
  int lvl = level < 0 ? depth() : level;
  unsigned int bit = (1u << lvl);
  int nmax = 0;
  for (scalar s0 in list)
    if (!is_constant(s0) && s0.i >= 0) {
      scalar s = s0.block < 0 ? (scalar){s0.i + s0.block} : s0;
      if (s.block > 0)
	nmax += s.block;
    }
  if (!nmax)
    return;

  int * fields = (int *) malloc (nmax*sizeof(int)), nfields = 0;
  for (scalar s0 in list)
    if (!is_constant(s0) && s0.i >= 0) {
      scalar s = s0.block < 0 ? (scalar){s0.i + s0.block} : s0;
      if (s.block <= 0)
	continue;
      if (!(s.stencil.bc & ((1 << 0) | (1 << 1) | (1 << 6))))
	s.gpu.mpi_clean = 0;
      if (s.gpu.mpi_clean & bit)
	continue;
      if (s.gpu.stored > 0)
	gpu_cpu_sync_scalar (s.i, s.block, multigrid->d,
			     multigrid->shift[depth() + 1], GPU_WRITE);
      s.gpu.stored = -1;
      s.gpu.mpi_clean |= bit;
      for (int b = 0; b < s.block; b++)
	fields[nfields++] = s.i + b;
    }
  if (!nfields) {
    free (fields);
    return;
  }

  prof_start ("mpi_boundary_level_gpu");

  int nx = (1 << lvl)*LocalDims.x;
  int ny = (1 << lvl)*LocalDims.y;
#if dimension == 3
  int nz = (1 << lvl)*LocalDims.z;
#endif
  int nx_tot = nx + 2*GHOSTS;
  int ny_tot = ny + 2*GHOSTS;
#if dimension == 3
  int nz_tot = nz + 2*GHOSTS;
#else
  int nz_tot = 1;
#endif
  size_t fsize = multigrid->shift[depth() + 1];
  size_t shift = multigrid->shift[lvl];
  int stage_host = !use_cuda_aware_mpi ();

  for (int axis = 0; axis < dimension; axis++) {
    int left = MPI_PROC_NULL, right = MPI_PROC_NULL;
    MPI_Cart_shift (current_mpi_boundary->cartcomm, axis, 1, &left, &right);
    if (left == MPI_PROC_NULL && right == MPI_PROC_NULL)
      continue;

    int ci = nx_tot, cj = ny_tot, ck = nz_tot;
    int sr_i = 0, sr_j = 0, sr_k = 0;
    int sl_i = 0, sl_j = 0, sl_k = 0;
    int rl_i = 0, rl_j = 0, rl_k = 0;
    int rr_i = 0, rr_j = 0, rr_k = 0;
    if (axis == 0) {
      ci = GHOSTS;
      sr_i = nx; sl_i = GHOSTS; rl_i = 0; rr_i = nx + GHOSTS;
    }
    else if (axis == 1) {
      cj = GHOSTS;
      sr_j = ny; sl_j = GHOSTS; rl_j = 0; rr_j = ny + GHOSTS;
    }
#if dimension == 3
    else {
      ck = GHOSTS;
      sr_k = nz; sl_k = GHOSTS; rl_k = 0; rr_k = nz + GHOSTS;
    }
#endif

    if (left == pid() && right == pid()) {
      size_t nb0 = cuda_halo_pack (fields, nfields, fsize, shift, ny_tot, nz_tot,
				   sr_i, sr_j, sr_k, ci, cj, ck, 0);
      size_t nb1 = cuda_halo_pack (fields, nfields, fsize, shift, ny_tot, nz_tot,
				   sl_i, sl_j, sl_k, ci, cj, ck, 1);
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rl_i, rl_j, rl_k, ci, cj, ck, 0, -1, nb0);
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rr_i, rr_j, rr_k, ci, cj, ck, 1, -1, nb1);
      cuda_halo_sync ();
      continue;
    }

    size_t nbytes = (size_t)nfields*ci*cj*ck*sizeof(real);
    MPI_Request reqs[4];
    int nreq = 0;

    if (left != MPI_PROC_NULL) {
      void * p_rl = cuda_halo_rcv_ptr (0, stage_host, nbytes);
      MPI_Irecv (p_rl, nbytes, MPI_BYTE, left, 0,
		 current_mpi_boundary->cartcomm, &reqs[nreq++]);
    }
    if (right != MPI_PROC_NULL) {
      void * p_rr = cuda_halo_rcv_ptr (1, stage_host, nbytes);
      MPI_Irecv (p_rr, nbytes, MPI_BYTE, right, 1,
		 current_mpi_boundary->cartcomm, &reqs[nreq++]);
    }

    if (right != MPI_PROC_NULL)
      cuda_halo_pack (fields, nfields, fsize, shift, ny_tot, nz_tot,
		      sr_i, sr_j, sr_k, ci, cj, ck, 0);
    if (left != MPI_PROC_NULL)
      cuda_halo_pack (fields, nfields, fsize, shift, ny_tot, nz_tot,
		      sl_i, sl_j, sl_k, ci, cj, ck, 1);

    if (right != MPI_PROC_NULL) {
      void * p_sr = cuda_halo_snd_ptr (0, stage_host, nbytes);
      MPI_Isend (p_sr, nbytes, MPI_BYTE, right, 0,
		 current_mpi_boundary->cartcomm, &reqs[nreq++]);
    }
    if (left != MPI_PROC_NULL) {
      void * p_sl = cuda_halo_snd_ptr (1, stage_host, nbytes);
      MPI_Isend (p_sl, nbytes, MPI_BYTE, left, 1,
		 current_mpi_boundary->cartcomm, &reqs[nreq++]);
    }

    MPI_Waitall (nreq, reqs, MPI_STATUSES_IGNORE);

    if (left != MPI_PROC_NULL)
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rl_i, rl_j, rl_k, ci, cj, ck, 0, stage_host, nbytes);
    if (right != MPI_PROC_NULL)
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rr_i, rr_j, rr_k, ci, cj, ck, 1, stage_host, nbytes);
    cuda_halo_sync ();
  }

  free (fields);
  prof_stop();
}
#endif // _GPU && _CUDA

trace
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{
  scalar * list_cpu = NULL;
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
	list_cpu = list_add (list_cpu, s);
    }

#if _GPU && _CUDA
  if (list_gpu) {
    mpi_boundary_level_gpu (list_gpu, level);
    free (list_gpu);
  }
#endif

  if (!list_cpu)
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
    if ((buf[0] = snd_x (npl - 2*GHOSTS, right, 0, level, list_cpu, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (GHOSTS, left, 1, level, list_cpu, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list_cpu);
    rcv_x (npl - GHOSTS, right, 1, level, list_cpu);
    MPI_Waitall (nr, reqs, MPI_STATUSES_IGNORE);
    free (buf[0]); free (buf[1]);
  }

  free (list_cpu);

  prof_stop();
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  if (current_mpi_boundary == m)
    current_mpi_boundary = NULL;
  MPI_Comm_free (&m->cartcomm);
  free (m);
}

static int mpi_gcd (int a, int b)
{
  while (b != 0) {
    int t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static bool eval_mpi_factorization (const int p[3], const int ud[3],
				    const int per[3], int k_cells,
				    int * out_S, double * out_cost)
{
  int S = 1;
  for (int d = 0; d < dimension; d++) {
    int g = mpi_gcd (ud[d], p[d]);
    int r = p[d]/g;
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

  double ld[3] = {1., 1., 1.};
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
  current_mpi_boundary = m;

  int ud[3] = { UserDims.x,
#if dimension > 1
		UserDims.y,
#else
		1,
#endif
#if dimension > 2
		UserDims.z
#else
		1
#endif
  };
  int per[3] = { Period.x,
#if dimension > 1
		 Period.y,
#else
		 0,
#endif
#if dimension > 2
		 Period.z
#else
		 0
#endif
  };

  int k_cells = 1 << grid->depth;
  int np = npe();
  int best_p[3] = {0, 0, 0}, best_S = 1 << 30;
  double best_cost = 1e300;

  int user_p[3] = { UserMpiDims.x,
#if dimension > 1
		    UserMpiDims.y,
#else
		    1,
#endif
#if dimension > 2
		    UserMpiDims.z
#else
		    1
#endif
  };

  if (user_p[0] > 0) {
    int prod = 1;
    for (int d = 0; d < dimension; d++) prod *= user_p[d];
    if (prod != np ||
	!eval_mpi_factorization (user_p, ud, per, k_cells, &best_S, &best_cost)) {
      if (pid() == 0)
	fprintf (stderr,
		 "%s:%d: error: requested mpi_dims(%d,%d,%d) is incompatible"
		 " with npe=%d, dimensions(%d,%d,%d) and N=%d\n",
		 __FILE__, LINENO, user_p[0], user_p[1], user_p[2],
		 np, ud[0], ud[1], ud[2], N);
      exit (1);
    }
    memcpy (best_p, user_p, sizeof(best_p));
  }
  else {
    for (int px = 1; px <= np; px++) {
      if (np % px) continue;
      int rem_x = np/px;
      int py_max = dimension >= 2 ? rem_x : 1;
      for (int py = 1; py <= py_max; py++) {
	if (rem_x % py) continue;
	int pz = dimension >= 3 ? (rem_x/py) : 1;
	if (px*py*pz != np) continue;
	int cand_p[3] = {px, py, pz}, cand_S = 0;
	double cand_cost = 0.;
	if (eval_mpi_factorization (cand_p, ud, per, k_cells, &cand_S, &cand_cost)) {
	  if (cand_S < best_S ||
	      (cand_S == best_S && cand_cost < best_cost - 1e-12)) {
	    best_S = cand_S;
	    best_cost = cand_cost;
	    memcpy (best_p, cand_p, sizeof(best_p));
	  }
	}
      }
    }
    if (best_p[0] == 0) {
      if (pid() == 0)
	fprintf (stderr,
		 "%s:%d: error: cannot factor npe=%d across dimensions(%d,%d,%d)"
		 " at N=%d\n",
		 __FILE__, LINENO, np, ud[0], ud[1], ud[2], N);
      exit (1);
    }
  }

  MpiDims.x = best_p[0];
  Dimensions.x = ud[0]*best_S;
  LocalDims.x = Dimensions.x/MpiDims.x;
#if dimension > 1
  MpiDims.y = best_p[1];
  Dimensions.y = ud[1]*best_S;
  LocalDims.y = Dimensions.y/MpiDims.y;
#endif
#if dimension > 2
  MpiDims.z = best_p[2];
  Dimensions.z = ud[2]*best_S;
  LocalDims.z = Dimensions.z/MpiDims.z;
#endif

  MPI_Cart_create (MPI_COMM_WORLD, dimension,
		   &MpiDims.x, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), dimension, mpi_coords);

  mpi_off0 = mpi_coords[0]*LocalDims.x;
#if dimension > 1
  mpi_off1 = mpi_coords[1]*LocalDims.y;
#endif
#if dimension > 2
  mpi_off2 = mpi_coords[2]*LocalDims.z;
#endif

  for (int i = 0; i < 6; i++)
    mpi_neighbor[i] = -1;
  struct { int x, y, z; } dir = {0,1,2};
  foreach_dimension() {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    mpi_neighbor[left]  = (l == MPI_PROC_NULL ? -1 : l);
    mpi_neighbor[right] = (r == MPI_PROC_NULL ? -1 : r);
  }

  Dimensions_scale = Dimensions.x;
  int rem = k_cells < best_S ? 1 : (k_cells/best_S), r = 0;
  while (rem > 1)
    rem /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = Dimensions.x*(1 << r);

  grid->n = 1L << (dimension*depth());
  foreach_dimension()
    grid->n *= LocalDims.x;
  grid->tn = npe()*grid->n;

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

  for_each_root (LocalDims) {
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
