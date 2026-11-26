/**
# GPU-resident MPI halo exchange for CUDA multigrids

Included from `grid/multigrid-mpi.h` when `_GPU && _CUDA` is enabled.
Exchanges ghost layers directly between GPU buffers (via CUDA-aware MPI
when available, or pinned-host staging otherwise) using the pack/unpack
kernels in `grid/cuda/cuda-mpi.c`. */

static MpiBoundary * gpu_mpi_boundary = NULL;

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
  if (!gpu_mpi_boundary || !list)
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

  int * fields = qmalloc (nmax, int), nfields = 0;
  for (scalar s0 in list)
    if (!is_constant(s0) && s0.i >= 0) {
      scalar s = s0.block < 0 ? (scalar){s0.i + s0.block} : s0;
      if (s.block <= 0)
	continue;
      if (!(s.stencil.bc & (s_centered | s_face | s_phys_gpu)))
	s.gpu.mpi_clean = 0;
      if (s.gpu.mpi_clean & bit)
	continue;
      if (s.gpu.stored > 0)
	gpu_cpu_sync_scalar (s.i, s.block, grid_data(),
			     field_size(), GPU_WRITE);
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
  size_t fsize = field_size();
  size_t shift = shift_level(lvl);
  int stage_host = !use_cuda_aware_mpi ();

  for (int axis = 0; axis < dimension; axis++) {
    int left = MPI_PROC_NULL, right = MPI_PROC_NULL;
    MPI_Cart_shift (gpu_mpi_boundary->cartcomm, axis, 1, &left, &right);
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
      size_t nb0 = cuda_halo_pack (fields, nfields, fsize, shift,
				   ny_tot, nz_tot, sr_i, sr_j, sr_k,
				   ci, cj, ck, 0);
      size_t nb1 = cuda_halo_pack (fields, nfields, fsize, shift,
				   ny_tot, nz_tot, sl_i, sl_j, sl_k,
				   ci, cj, ck, 1);
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rl_i, rl_j, rl_k, ci, cj, ck, 0, -1, nb0);
      cuda_halo_unpack (fields, nfields, fsize, shift, ny_tot, nz_tot,
			rr_i, rr_j, rr_k, ci, cj, ck, 1, -1, nb1);
      cuda_halo_sync ();
      continue;
    }

    size_t nbytes = (size_t) nfields*ci*cj*ck*sizeof(real);
    MPI_Request reqs[4];
    int nreq = 0;

    if (left != MPI_PROC_NULL) {
      void * p_rl = cuda_halo_rcv_ptr (0, stage_host, nbytes);
      MPI_Irecv (p_rl, nbytes, MPI_BYTE, left, 0,
		 gpu_mpi_boundary->cartcomm, &reqs[nreq++]);
    }
    if (right != MPI_PROC_NULL) {
      void * p_rr = cuda_halo_rcv_ptr (1, stage_host, nbytes);
      MPI_Irecv (p_rr, nbytes, MPI_BYTE, right, 1,
		 gpu_mpi_boundary->cartcomm, &reqs[nreq++]);
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
		 gpu_mpi_boundary->cartcomm, &reqs[nreq++]);
    }
    if (left != MPI_PROC_NULL) {
      void * p_sl = cuda_halo_snd_ptr (1, stage_host, nbytes);
      MPI_Isend (p_sl, nbytes, MPI_BYTE, left, 1,
		 gpu_mpi_boundary->cartcomm, &reqs[nreq++]);
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
