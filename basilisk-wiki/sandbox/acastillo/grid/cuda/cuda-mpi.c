/**
# Multi-GPU MPI support for the CUDA backend

This file is included by `cuda.c` and provides:

1. Per-rank GPU device binding (`cuda_set_device`) and local MPI
   subdomain metadata (`cuda_set_mpi_domain`) so `libbuda.a` /
   `libbuda64.a` remain free of any compile-time `mpi.h` dependency.
2. Runtime detection of CUDA-aware MPI (`cuda_query_mpi_support`).
3. GPU-resident halo packing and unpacking kernels (`cuda_halo_pack`,
   `cuda_halo_unpack`) with reusable device and pinned-host staging
   buffers (`cuda_halo_snd_ptr`, `cuda_halo_rcv_ptr`, `cuda_halo_sync`).
*/

void cuda_set_device (int local_rank)
{
  cuda_local_rank = local_rank;
}

void cuda_set_mpi_domain (int lx, int ly, int lz, int ox, int oy, int oz)
{
  cuda_local_dims.x = lx; cuda_local_dims.y = ly; cuda_local_dims.z = lz;
  cuda_mpi_off[0] = ox; cuda_mpi_off[1] = oy; cuda_mpi_off[2] = oz;
}

int cuda_query_mpi_support (void)
{
  const char * env = getenv ("BASILISK_CUDA_MPI");
  if (env)
    return atoi (env) != 0;
  typedef int (* query_fn_t) (void);
  query_fn_t fn = (query_fn_t) dlsym (RTLD_DEFAULT, "MPIX_Query_cuda_support");
  if (fn)
    return fn () == 1;
  return 0;
}

/**
## GPU-resident halo pack & unpack kernels */

static CUdeviceptr halo_d_snd[2] = {0, 0}, halo_d_rcv[2] = {0, 0};
static CUdeviceptr halo_d_fields = 0;
static void * halo_h_snd[2] = {NULL, NULL}, * halo_h_rcv[2] = {NULL, NULL};
static size_t halo_snd_cap[2] = {0, 0}, halo_rcv_cap[2] = {0, 0};
static size_t halo_hsnd_cap[2] = {0, 0}, halo_hrcv_cap[2] = {0, 0};
static size_t halo_fields_cap = 0;
static CUfunction k_halo_pack = 0, k_halo_unpack = 0;

static void cuda_halo_free (void)
{
  for (int s = 0; s < 2; s++) {
    if (halo_d_snd[s]) {
      CUDA_CHECK (cuMemFree (halo_d_snd[s]));
      halo_d_snd[s] = 0; halo_snd_cap[s] = 0;
    }
    if (halo_d_rcv[s]) {
      CUDA_CHECK (cuMemFree (halo_d_rcv[s]));
      halo_d_rcv[s] = 0; halo_rcv_cap[s] = 0;
    }
    if (halo_h_snd[s]) {
      CUDA_CHECK (cuMemFreeHost (halo_h_snd[s]));
      halo_h_snd[s] = NULL; halo_hsnd_cap[s] = 0;
    }
    if (halo_h_rcv[s]) {
      CUDA_CHECK (cuMemFreeHost (halo_h_rcv[s]));
      halo_h_rcv[s] = NULL; halo_hrcv_cap[s] = 0;
    }
  }
  if (halo_d_fields) {
    CUDA_CHECK (cuMemFree (halo_d_fields));
    halo_d_fields = 0; halo_fields_cap = 0;
  }
}

static void ensure_halo_kernels (void)
{
  if (k_halo_pack && k_halo_unpack)
    return;
  static const char halo_src[] =
    "typedef " CUDA_REAL_STR " rtype;\n"
    "extern \"C\" __global__\n"
    "void halo_pack (const rtype * __restrict__ data,\n"
    "                rtype * __restrict__ buf,\n"
    "                const int * __restrict__ fields, int nfields,\n"
    "                size_t fsize, size_t shift, int ny_tot, int nz_tot,\n"
    "                int i0, int j0, int k0, int ci, int cj, int ck) {\n"
    "  size_t idx = (size_t) blockIdx.x*blockDim.x + threadIdx.x;\n"
    "  size_t plane = (size_t) ci*cj*ck;\n"
    "  size_t total = plane*(size_t) nfields;\n"
    "  if (idx >= total) return;\n"
    "  int f = (int)(idx / plane);\n"
    "  size_t rem = idx % plane;\n"
    "  int dk = (int)(rem % (size_t) ck); rem /= (size_t) ck;\n"
    "  int dj = (int)(rem % (size_t) cj);\n"
    "  int di = (int)(rem / (size_t) cj);\n"
    "  size_t off = (size_t) fields[f]*fsize + shift +\n"
    "    (size_t) nz_tot*((size_t)(j0 + dj) +\n"
    "    (size_t)(i0 + di)*(size_t) ny_tot) + (size_t)(k0 + dk);\n"
    "  buf[idx] = data[off];\n"
    "}\n"
    "extern \"C\" __global__\n"
    "void halo_unpack (rtype * __restrict__ data,\n"
    "                  const rtype * __restrict__ buf,\n"
    "                  const int * __restrict__ fields, int nfields,\n"
    "                  size_t fsize, size_t shift, int ny_tot, int nz_tot,\n"
    "                  int i0, int j0, int k0, int ci, int cj, int ck) {\n"
    "  size_t idx = (size_t) blockIdx.x*blockDim.x + threadIdx.x;\n"
    "  size_t plane = (size_t) ci*cj*ck;\n"
    "  size_t total = plane*(size_t) nfields;\n"
    "  if (idx >= total) return;\n"
    "  int f = (int)(idx / plane);\n"
    "  size_t rem = idx % plane;\n"
    "  int dk = (int)(rem % (size_t) ck); rem /= (size_t) ck;\n"
    "  int dj = (int)(rem % (size_t) cj);\n"
    "  int di = (int)(rem / (size_t) cj);\n"
    "  size_t off = (size_t) fields[f]*fsize + shift +\n"
    "    (size_t) nz_tot*((size_t)(j0 + dj) +\n"
    "    (size_t)(i0 + di)*(size_t) ny_tot) + (size_t)(k0 + dk);\n"
    "  data[off] = buf[idx];\n"
    "}\n";
  size_t ptx_size;
  char arch[] = "--gpu-architecture=compute_????";
  architecture (arch);
  char * ptx = compile_ptx (halo_src, arch, "halo_pack",
                            __FILE__, __LINE__, &ptx_size);
  CUmodule mod;
  CUDA_CHECK (cuModuleLoadData (&mod, ptx));
  free (ptx);
  CUDA_CHECK (cuModuleGetFunction (&k_halo_pack, mod, "halo_pack"));
  CUDA_CHECK (cuModuleGetFunction (&k_halo_unpack, mod, "halo_unpack"));
}

static void upload_halo_fields (const int * fields, int nfields)
{
  size_t need = (size_t) nfields*sizeof(int);
  if (need > halo_fields_cap) {
    if (halo_d_fields)
      CUDA_CHECK (cuMemFree (halo_d_fields));
    CUDA_CHECK (cuMemAlloc (&halo_d_fields, need));
    halo_fields_cap = need;
  }
  CUDA_CHECK (cuMemcpyHtoDAsync (halo_d_fields, fields, need, stream));
}

size_t cuda_halo_pack (const int * fields, int nfields, size_t fsize,
                       size_t shift, int ny_tot, int nz_tot,
                       int i0, int j0, int k0, int ci, int cj, int ck,
                       int slot)
{
  ensure_halo_kernels ();
  size_t nelem = (size_t) nfields*ci*cj*ck;
  size_t nbytes = nelem*sizeof(real);
  if (nbytes > halo_snd_cap[slot]) {
    if (halo_d_snd[slot])
      CUDA_CHECK (cuMemFree (halo_d_snd[slot]));
    CUDA_CHECK (cuMemAlloc (&halo_d_snd[slot], nbytes));
    halo_snd_cap[slot] = nbytes;
  }
  upload_halo_fields (fields, nfields);
  const size_t threads = 256;
  size_t blocks = (nelem + threads - 1)/threads;
  void * args[] = { &ssbo, &halo_d_snd[slot], &halo_d_fields, &nfields,
                    &fsize, &shift, &ny_tot, &nz_tot,
                    &i0, &j0, &k0, &ci, &cj, &ck };
  CUDA_CHECK (cuLaunchKernel (k_halo_pack, blocks, 1, 1, threads, 1, 1, 0,
                              stream, args, NULL));
  return nbytes;
}

void * cuda_halo_snd_ptr (int slot, int to_host, size_t nbytes)
{
  if (to_host) {
    if (nbytes > halo_hsnd_cap[slot]) {
      if (halo_h_snd[slot])
        CUDA_CHECK (cuMemFreeHost (halo_h_snd[slot]));
      CUDA_CHECK (cuMemAllocHost (&halo_h_snd[slot], nbytes));
      halo_hsnd_cap[slot] = nbytes;
    }
    CUDA_CHECK (cuMemcpyDtoHAsync (halo_h_snd[slot], halo_d_snd[slot],
                                   nbytes, stream));
    CUDA_CHECK (cuStreamSynchronize (stream));
    return halo_h_snd[slot];
  }
  CUDA_CHECK (cuStreamSynchronize (stream));
  return (void *)(uintptr_t) halo_d_snd[slot];
}

void * cuda_halo_rcv_ptr (int slot, int from_host, size_t nbytes)
{
  if (nbytes > halo_rcv_cap[slot]) {
    if (halo_d_rcv[slot])
      CUDA_CHECK (cuMemFree (halo_d_rcv[slot]));
    CUDA_CHECK (cuMemAlloc (&halo_d_rcv[slot], nbytes));
    halo_rcv_cap[slot] = nbytes;
  }
  if (from_host) {
    if (nbytes > halo_hrcv_cap[slot]) {
      if (halo_h_rcv[slot])
        CUDA_CHECK (cuMemFreeHost (halo_h_rcv[slot]));
      CUDA_CHECK (cuMemAllocHost (&halo_h_rcv[slot], nbytes));
      halo_hrcv_cap[slot] = nbytes;
    }
    return halo_h_rcv[slot];
  }
  return (void *)(uintptr_t) halo_d_rcv[slot];
}

/* `from_host`: <0 = local device-to-device periodic copy (from halo_d_snd),
   0 = CUDA-aware MPI (already in halo_d_rcv), >0 = pinned-host staging. */
void cuda_halo_unpack (const int * fields, int nfields, size_t fsize,
                       size_t shift, int ny_tot, int nz_tot,
                       int i0, int j0, int k0, int ci, int cj, int ck,
                       int slot, int from_host, size_t nbytes)
{
  ensure_halo_kernels ();
  CUdeviceptr src_buf;
  if (from_host < 0)
    src_buf = halo_d_snd[slot];
  else {
    if (from_host > 0)
      CUDA_CHECK (cuMemcpyHtoDAsync (halo_d_rcv[slot], halo_h_rcv[slot],
                                     nbytes, stream));
    src_buf = halo_d_rcv[slot];
  }
  upload_halo_fields (fields, nfields);
  size_t nelem = (size_t) nfields*ci*cj*ck;
  const size_t threads = 256;
  size_t blocks = (nelem + threads - 1)/threads;
  void * args[] = { &ssbo, &src_buf, &halo_d_fields, &nfields,
                    &fsize, &shift, &ny_tot, &nz_tot,
                    &i0, &j0, &k0, &ci, &cj, &ck };
  CUDA_CHECK (cuLaunchKernel (k_halo_unpack, blocks, 1, 1, threads, 1, 1, 0,
                              stream, args, NULL));
}

void cuda_halo_sync (void)
{
  if (ctx)
    CUDA_CHECK (cuStreamSynchronize (stream));
}
