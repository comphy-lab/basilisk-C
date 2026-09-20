/**
# Sampling a plane or slab stack at a coarser level than the local grid

Level-based twins of [spectra_sample.h](spectra_sample.h)'s
`sample_scalar_plane()` and `sample_scalar_stack_sum()`. That header samples
through `locate()`, which always descends to the finest cell under the query
point -- correct when the lattice is finer than or equal to the local grid,
but not when the lattice is coarser (the diagnostic level) while part of the
domain is refined further (the interface). There a query point lands on a
child leaf, not the coarse-cell average, exactly the bias
`profiles_slab_restrict.h` fixes for slab profiles.

Both functions here apply the same fix: like `profile_scalar_slab()`,
restrict the sampled fields first, then bin at exactly `slablevel` with
`foreach_level`. `m1 = m2 = 2^slablevel` -- the lattice is the level, not a
caller-chosen resolution, since `foreach_level` visits cells, not points.

`slablevel` has no default, for the same reason `profiles_slab_restrict.h`
gives none: `foreach_level(l)` visited where the tree is coarser than `l`
would come back empty, indistinguishable from a hole.

Shares `snap_to_cell()` with [spectra_sample.h](spectra_sample.h) via
[spectra_utils.h](spectra_utils.h), so a test can include both headers
without a duplicate-symbol error.
*/

#if dimension != 3
  #error sample_scalar_plane_restrict() is 3D only
#endif

#include "spectra_utils.h"

/**
## Sampling a plane at a coarser level than the local grid

Restrict the sampled fields first, then bin at exactly `slablevel` with
`foreach_level`, filtered to the row of cells whose $z$ contains `h`. The
caller allocates `plane` with `m1*m2*len` doubles, laid out as
`sample_scalar_plane()`'s output. Returns the number of unfilled points, 0
for a complete plane. */

int sample_scalar_plane_restrict (scalar * list, double * plane, double h,
                                  int slablevel)
{
  int m1 = 1 << slablevel, m2 = m1, len = list_len (list), n = m1*m2*len;
  for (int i = 0; i < n; i++)
    plane[i] = nodata;

  restriction (list);

  double del = L0/m1;
  foreach_level (slablevel, reduction(min:plane[:n])) {
    if (fabs (z - h) > Delta/2.) continue;
    int i = (x - X0)/del, j = (y - Y0)/del;
    if (i < 0) i = 0; if (i >= m1) i = m1 - 1;
    if (j < 0) j = 0; if (j >= m2) j = m2 - 1;
    double * alias = plane;
    int k = 0;
    for (scalar s in list)
      alias[(i*m2 + j)*len + k++] = s[];
  }

  int holes = 0;
  for (int i = 0; i < n; i++)
    if (plane[i] == nodata)
      holes++;
  return holes;
}

/**
## Foliating a stack of planes at a coarser level than the local grid

One `restriction()` and one `foreach_level(slablevel)` pass over the whole
box, every cell binned by its own `z` into `iz` (alongside `x,y` into `i,j`)
rather than one `locate()`-based plane per slab
(`sample_scalar_stack_sum()`'s approach). `hmin`/`hmax` must align to the
`slablevel` lattice in `z` (`nz` slabs spanning a whole multiple of
`L0/2^slablevel`), so each `iz` is exactly one row of cells -- the same
alignment `profile_scalar_slab()` requires of its own `hmin`/`hmax`/`n`.

`means` is laid out as `means[iz*len + k]` for field `k` at slab `iz`, the
same convention `sample_scalar_stack_sum()` uses. `plane` accumulates
`sum_iz (q(x,y,z_iz) - means[iz,k])`, laid out as
`sample_scalar_plane_restrict()`'s output. `m1 = m2 = 2^slablevel`. Returns
the count of `(i,j,iz)` cells with no contributing tree cell, summed over
slabs -- like `sample_scalar_stack_sum()`, this can't distinguish an empty
point from one whose contributions summed to exactly zero. */

int sample_scalar_stack_sum_restrict (scalar * list, double * plane,
                                      const double * means, double * zout,
                                      double hmin, double hmax, int nz,
                                      int slablevel)
{
  int m1 = 1 << slablevel, m2 = m1, len = list_len (list), n = m1*m2*len;
  for (int i = 0; i < n; i++)
    plane[i] = 0.;
  for (int iz = 0; iz < nz; iz++)
    zout[iz] = snap_to_cell (nz > 1 ? hmin + (hmax - hmin)*(iz + 0.5)/nz : hmin,
                             m1);

  restriction (list);

  double del = L0/m1;
  double dz = nz > 1 ? (hmax - hmin)/nz : L0;
  int * count = malloc (nz*sizeof(int));
  for (int iz = 0; iz < nz; iz++) count[iz] = 0;

  foreach_level (slablevel, reduction(+:plane[:n]) reduction(+:count[:nz])) {
    int iz = nz > 1 ? (int)floor ((z - hmin)/dz) : 0;
    if (iz < 0 || iz >= nz) continue;
    int i = (x - X0)/del, j = (y - Y0)/del;
    if (i < 0) i = 0; if (i >= m1) i = m1 - 1;
    if (j < 0) j = 0; if (j >= m2) j = m2 - 1;
    int idx = i*m2 + j;
    count[iz]++;
    double * alias = plane;
    int k = 0;
    for (scalar s in list)
      alias[idx*len + k] += s[] - means[iz*len + k], k++;
  }

  int holes = 0;
  for (int iz = 0; iz < nz; iz++)
    if (count[iz] < m1*m2)
      holes += m1*m2 - count[iz];
  free (count);
  return holes;
}
