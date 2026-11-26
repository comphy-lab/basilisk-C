#ifndef DUMP_RESTORE_H
#define DUMP_RESTORE_H

#include "output.h"

#if MULTIGRID && defined(_LOCALDIMS)

/**
## Canonical multi-root *dump()* and *restore()* for MPI / GPU multigrids

Writes and reads the standard Basilisk `dump` binary format (byte-for-byte
identical to [src/output.h](/src/output.h#dump-basilisk-snapshots)) using C
`stdio` buffered one `level = 0` root subtree at a time:

* **Canonical root order**: `(ox, oy, oz)` lexicographic with
  `header.n = Dimensions` and `header.npe = nx*ny*nz`, producing bit-identical
  files regardless of `npe()` or `MpiDims`.
* **Zero temporary GPU scalars**: Evaluates uniform-multigrid subtree sizes in
  closed form without allocating `scalar size[]` or `scalar index[]` (avoiding
  GPU SSBO reallocation and `SINGLE_PRECISION` 24-bit mantissa overflow above
  $2^{24} \approx 16.7\times 10^6$ cells).
* **Cross-rank and cross-backend restarts**: Checkpoints written in `FP64` or
  `FP32` on any GPU/CPU rank count restore onto any valid GPU/CPU rank count,
  stock CPU `multigrid.h` / `multigrid3D.h`, or MPI `acastillo/grid/octree.h`.

> **Note**: When combined with `maxruntime.h`, include this header *before*
> `#include "maxruntime.h"` so `maxruntime.h`'s checkpoint call also resolves
> to `dump_multigrid()`.
*/

static inline long multigrid_subtree_cells (int h)
{
  return ((1L << (dimension*(h + 1))) - 1)/((1L << dimension) - 1);
}

static long dump_multigrid_header (FILE * fp, const struct DumpHeader * header,
				   scalar * dlist)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  const char * size_name = "size";
  unsigned slen = 4;
  if (fwrite (&slen, sizeof(unsigned), 1, fp) < 1 ||
      fwrite (size_name, sizeof(char), slen, fp) < slen) {
    perror ("dump(): error while writing s.name");
    exit (1);
  }
  long nbytes = sizeof(struct DumpHeader) + sizeof(unsigned) + slen;
  for (scalar s in dlist) {
    unsigned len = strlen (s.name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1 ||
	fwrite (s.name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
    nbytes += sizeof(unsigned) + len;
  }
  double o[4] = {X0, Y0, Z0, L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
  return nbytes + 4*sizeof(double);
}

static inline char * pack_multigrid_root (Point root, scalar * dlist, char * buf)
{
  char * p = buf;
  foreach_cell_root (root) {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    memcpy (p, &flags, sizeof(unsigned));
    p += sizeof(unsigned);
    double sz = (double) multigrid_subtree_cells (depth() - point.level);
    memcpy (p, &sz, sizeof(double));
    p += sizeof(double);
    for (scalar s in dlist) {
      double val = s[];
      memcpy (p, &val, sizeof(double));
      p += sizeof(double);
    }
    if (is_leaf(cell))
      continue;
  }
  return p;
}

trace
void dump_multigrid (const char * file = "dump",
		     scalar * list = all,
		     FILE * fp = NULL,
		     bool unbuffered = false,
		     bool zero = true)
{
#if _MPI
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit (1);
  }
#endif

  scalar * dlist = dump_list (list, zero);
#if _GPU
  for (scalar s in dlist)
    if (s.gpu.stored < 0) {
      gpu_cpu_sync_scalar (s.i, s.block, grid_data(), field_size(), GPU_READ);
      s.gpu.stored = 0;
    }
#endif

  int nroots_tot = 1;
  struct DumpHeader header = {
    t, (long)(1 + list_len (dlist)), iter, depth(), 1, dump_version
  };
  foreach_dimension() {
    header.n.x = Dimensions.x;
    nroots_tot *= Dimensions.x;
  }
  header.npe = nroots_tot;

  long sizeofheader = sizeof(struct DumpHeader) + 4*sizeof(double) +
		      sizeof(unsigned) + 4;
  for (scalar s in dlist)
    sizeofheader += sizeof(unsigned) + strlen (s.name);

  char * name = NULL;
  FILE * fh = fp;
  if (!fh) {
    name = (char *) malloc (strlen (file) + 2);
    strcpy (name, file);
    if (!unbuffered)
      strcat (name, "~");
    if (pid() == 0) {
      if ((fh = fopen (name, "wb")) == NULL) {
	perror (name);
	exit (1);
      }
      long written = dump_multigrid_header (fh, &header, dlist);
      assert (written == sizeofheader);
      fflush (fh);
    }
#if _MPI
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() != 0 && (fh = fopen (name, "r+b")) == NULL) {
      perror (name);
      exit (1);
    }
#endif
  }
  else {
    long written = dump_multigrid_header (fh, &header, dlist);
    assert (written == sizeofheader);
  }

  long cells_per_root = multigrid_subtree_cells (depth());
  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  size_t root_bytes = (size_t)cells_per_root*(size_t)cell_size;
  char * root_buf = (char *) malloc (root_bytes);
  if (!root_buf) {
    perror ("dump(): malloc");
    exit (1);
  }

  long pos = (pid() == 0 && fp == NULL) ? sizeofheader : -1;
  for_each_root (LocalDims) {
    char * endp = pack_multigrid_root (root, dlist, root_buf);
    assert ((size_t)(endp - root_buf) == root_bytes);
    if (fp == NULL) {
      long offset = sizeofheader + global_root_index (root)*(long)root_bytes;
      if (pos != offset) {
	if (fseek (fh, offset, SEEK_SET) < 0) {
	  perror ("dump(): error while seeking");
	  exit (1);
	}
	pos = offset;
      }
    }
    if (fwrite (root_buf, 1, root_bytes, fh) < root_bytes) {
      perror ("dump(): error while writing scalars");
      exit (1);
    }
    if (pos >= 0)
      pos += (long)root_bytes;
  }

  free (root_buf);
  free (dlist);

  if (!fp) {
    fflush (fh);
    fclose (fh);
#if _MPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    if (!unbuffered && pid() == 0)
      rename (name, file);
    free (name);
  }
  else
    fflush (fh);
}

static inline const char * unpack_multigrid_root (Point root, scalar * slist,
						  const char * buf)
{
  const char * p = buf;
  foreach_cell_root (root) {
    p += sizeof(unsigned) + sizeof(double); // skip flags + subtree size
    for (scalar s in slist) {
      double val;
      memcpy (&val, p, sizeof(double));
      p += sizeof(double);
      if (s.i != INT_MAX)
	s[] = isfinite(val) ? val : nodata;
    }
    if (is_leaf(cell))
      continue;
  }
  return p;
}

trace
bool restore_multigrid (const char * file = "dump",
			scalar * list = NULL,
			FILE * fp = NULL)
{
  FILE * fh = fp;
  if (!fh && (fh = fopen (file, "rb")) == NULL)
    return false;
  assert (fh);

  struct DumpHeader header = {0};
  if (fread (&header, sizeof(header), 1, fh) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }

  int fnx = header.n.x >= 1. ? (int)(header.n.x + 0.5) : 1;
  int fny = (dimension >= 2 && header.n.y >= 1.) ? (int)(header.n.y + 0.5) : 1;
  int fnz = (dimension >= 3 && header.n.z >= 1.) ? (int)(header.n.z + 0.5) : 1;
  int target_n = (1 << header.depth)*fnx;

  bool same_grid = (grid != NULL && N == target_n &&
		    Dimensions.x*(1 << depth()) == target_n
#if dimension >= 2
		    && Dimensions.y*(1 << depth()) == (1 << header.depth)*fny
#endif
#if dimension >= 3
		    && Dimensions.z*(1 << depth()) == (1 << header.depth)*fnz
#endif
		    );
  if (!same_grid) {
    dimensions (fnx, fny, fnz);
    if (grid)
      free_grid();
    init_grid (target_n);
  }

  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all, true);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
	       "restore(): error: the list lengths don't match: "
	       "%ld (file) != %d (code)\n",
	       header.len - 1, list_len (slist));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
	       "restore(): error: file version mismatch: "
	       "%d (file) != %d (code)\n",
	       header.version, dump_version);
      exit (1);
    }
    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fh) < 1) {
	fprintf (ferr, "restore(): error: expecting len\n");
	exit (1);
      }
      char sname[len + 1];
      if (fread (sname, sizeof(char), len, fh) < len) {
	fprintf (ferr, "restore(): error: expecting s.name\n");
	exit (1);
      }
      sname[len] = '\0';
      if (i > 0) { // skip subtree size
	bool found = false;
	for (scalar s in slist)
	  if (!strcmp (s.name, sname)) {
	    input = list_append (input, s);
	    found = true; break;
	  }
	if (!found) {
	  if (restore_all) {
	    scalar s = new scalar;
	    free (s.name);
	    s.name = strdup (sname);
	    input = list_append (input, s);
	  }
	  else
	    input = list_append (input, (scalar){INT_MAX});
	}
      }
    }
    free (slist);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fh) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }

  long data_start = ftell (fh);
  long fcells_per_root = multigrid_subtree_cells (header.depth);
  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  size_t froot_bytes = (size_t)fcells_per_root*(size_t)cell_size;
  char * root_buf = (char *) malloc (froot_bytes);
  if (!root_buf) {
    perror ("restore(): malloc");
    exit (1);
  }

  int shift_lvl = header.depth - depth();
  assert (shift_lvl >= 0);

  if (shift_lvl == 0) {
    long pos = data_start;
    for_each_root (LocalDims) {
      long r_glob = global_root_index (root);
      long offset = data_start + r_glob*(long)froot_bytes;
      if (pos != offset) {
	if (fseek (fh, offset, SEEK_SET) < 0) {
	  perror ("restore(): error while seeking");
	  exit (1);
	}
	pos = offset;
      }
      if (fread (root_buf, 1, froot_bytes, fh) < froot_bytes) {
	fprintf (ferr, "restore(): error: expecting scalars\n");
	exit (1);
      }
      pos += (long)froot_bytes;
      const char * endp = unpack_multigrid_root (root, slist, root_buf);
      assert ((size_t)(endp - root_buf) == froot_bytes);
    }
  }
  else {
    // Coarser root factor in file than current MPI split (e.g. 1x1x1 -> 1x2x2)
    for (int fox = 0; fox < fnx; fox++)
      for (int foy = 0; foy < fny; foy++)
	for (int foz = 0; foz < fnz; foz++) {
	  long r_file = ((long)fox*fny + foy)*fnz + foz;
	  long offset = data_start + r_file*(long)froot_bytes;
	  if (fseek (fh, offset, SEEK_SET) < 0 ||
	      fread (root_buf, 1, froot_bytes, fh) < froot_bytes) {
	    fprintf (ferr, "restore(): error: expecting scalars\n");
	    exit (1);
	  }
	  const char * p = root_buf;
	  struct { int l, i, j, k, stg; } stk[STACKSIZE];
	  int sp = 0;
	  stk[0].l = 0; stk[0].i = GHOSTS; stk[0].j = GHOSTS;
	  stk[0].k = GHOSTS; stk[0].stg = 0;
	  while (sp >= 0) {
	    int fl = stk[sp].l, fi = stk[sp].i, fj = stk[sp].j,
		fk = stk[sp].k, stg = stk[sp].stg;
	    sp--;
	    if (stg == 0) {
	      p += sizeof(unsigned) + sizeof(double);
	      if (fl >= shift_lvl) {
		int loc_l = fl - shift_lvl;
		int li = fox*(1 << fl) + (fi - GHOSTS) - MpiOffset.x*(1 << loc_l);
#if dimension >= 2
		int lj = foy*(1 << fl) + (fj - GHOSTS) - MpiOffset.y*(1 << loc_l);
#endif
#if dimension >= 3
		int lk = foz*(1 << fl) + (fk - GHOSTS) - MpiOffset.z*(1 << loc_l);
#endif
		bool inside = (li >= 0 && li < LocalDims.x*(1 << loc_l))
#if dimension >= 2
			   && (lj >= 0 && lj < LocalDims.y*(1 << loc_l))
#endif
#if dimension >= 3
			   && (lk >= 0 && lk < LocalDims.z*(1 << loc_l))
#endif
		;
		Point point = {0};
		point.level = loc_l;
		point.i = li + GHOSTS;
#if dimension >= 2
		point.j = lj + GHOSTS;
#endif
#if dimension >= 3
		point.k = lk + GHOSTS;
#endif
		for (scalar s in slist) {
		  double val;
		  memcpy (&val, p, sizeof(double));
		  p += sizeof(double);
		  if (inside && s.i != INT_MAX)
		    s[] = isfinite(val) ? val : nodata;
		}
	      }
	      else
		p += (header.len - 1)*sizeof(double);
	      if (fl < header.depth) {
		sp++;
		stk[sp].l = fl; stk[sp].i = fi; stk[sp].j = fj;
		stk[sp].k = fk; stk[sp].stg = 1;
		sp++;
		stk[sp].l = fl + 1; stk[sp].i = 2*fi - GHOSTS;
		stk[sp].j = 2*fj - GHOSTS; stk[sp].k = 2*fk - GHOSTS;
		stk[sp].stg = 0;
	      }
	    }
	    else {
#if dimension == 1
	      sp++;
	      stk[sp].l = fl + 1; stk[sp].i = 2*fi - GHOSTS + 1;
	      stk[sp].j = GHOSTS; stk[sp].k = GHOSTS; stk[sp].stg = 0;
#elif dimension == 2
	      if (stg < 3) {
		sp++;
		stk[sp].l = fl; stk[sp].i = fi; stk[sp].j = fj;
		stk[sp].k = fk; stk[sp].stg = stg + 1;
	      }
	      int cx = (stg >= 2) ? 1 : 0, cy = (stg == 1 || stg == 3) ? 1 : 0;
	      sp++;
	      stk[sp].l = fl + 1; stk[sp].i = 2*fi - GHOSTS + cx;
	      stk[sp].j = 2*fj - GHOSTS + cy; stk[sp].k = GHOSTS;
	      stk[sp].stg = 0;
#else
	      if (stg < 7) {
		sp++;
		stk[sp].l = fl; stk[sp].i = fi; stk[sp].j = fj;
		stk[sp].k = fk; stk[sp].stg = stg + 1;
	      }
	      int cx = (stg >> 2) & 1, cy = (stg >> 1) & 1, cz = stg & 1;
	      sp++;
	      stk[sp].l = fl + 1; stk[sp].i = 2*fi - GHOSTS + cx;
	      stk[sp].j = 2*fj - GHOSTS + cy; stk[sp].k = 2*fk - GHOSTS + cz;
	      stk[sp].stg = 0;
#endif
	    }
	  }
	  assert ((size_t)(p - root_buf) == froot_bytes);
	}
  }

  free (root_buf);

#if _GPU
  for (scalar s in slist)
    if (s.i != INT_MAX)
      s.gpu.stored = 1; // stored on CPU
#endif
  for (scalar s in all) {
    set_dirty_stencil (s);
#if _GPU && _MPI
    s.gpu.mpi_clean = 0;
#endif
  }

  scalar * listm = is_constant(cm) ? NULL : (scalar *){fm};
  scalar * other = NULL;
  for (scalar s in all)
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  free (other);

  free (slist);
  if (file && !fp)
    fclose (fh);

  // the events are advanced to catch up with the time
  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  return true;
}

#define dump    dump_multigrid
#define restore restore_multigrid

#endif // MULTIGRID && _LOCALDIMS
#endif // DUMP_RESTORE_H
