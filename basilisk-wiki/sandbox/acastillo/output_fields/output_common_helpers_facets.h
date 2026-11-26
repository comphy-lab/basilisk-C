#ifndef OUTPUT_COMMON_HELPERS_FACETS_H
#define OUTPUT_COMMON_HELPERS_FACETS_H

#include "geometry.h"
#include "fractions.h"
#include "acastillo/output_fields/output_common_helpers.h"

#if dimension == 1
coord mycs(Point point, scalar c)
{
  coord n = {1.};
  return n;
}
#elif dimension == 2
#include "myc2d.h"
#define mfacets int m = facets(n, alpha, v);
#else // dimension == 3
#include "myc.h"
#define mfacets int m = facets(n, alpha, v, 1.);
#endif

/** Macro to simplify facets calculation */ 
#define shortcut_facets               \
  coord n;                            \
  n = mycs(point, c);                 \
  double alpha = plane_alpha(c[], n); \
  coord v[12];                        \
  mfacets;

/**
## Other functions

### Synchronize VOF and facet scalar fields for host-side extraction

Ensures boundary/MPI halos of *c* are current (needed by the $3^d$ stencil in
`mycs`) and, on GPU backends, downloads *c* and *slist* to host memory before
`foreach (serial, noauto)` facet extraction. */

static inline void sync_facets_fields (scalar c, scalar * slist) {
  boundary ({c});
#if _GPU
  if (!is_constant (c) && c.gpu.stored < 0)
    gpu_cpu_sync_scalar (c.i, c.block, grid_data(), field_size(), GPU_READ);
  for (scalar _s in slist)
    if (!is_constant (_s) && _s.gpu.stored < 0)
      gpu_cpu_sync_scalar (_s.i, _s.block, grid_data(), field_size(), GPU_READ);
#endif
}

/** ### Count the number of vertices and facets */

void count_vertices_and_facets(scalar c, long *nverts, long *nfacets) {
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets
      for (int i = 0; i < m; i++)(*nverts)++;
      if (m > 0)
        (*nfacets)++;
    }
  }
}

/** ### Count the number of vertices and facets within a box

Cell-center-in-box variant of `count_vertices_and_facets()`. */
void count_vertices_and_facets_box(scalar c, coord box[2], long *nverts, long *nfacets) {
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1 && IN_BOX_CC(box))
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6 && IN_BOX_CC(box))
    #endif
    {
      shortcut_facets
      for (int i = 0; i < m; i++)(*nverts)++;
      if (m > 0)
        (*nfacets)++;
    }
  }
}

#endif // OUTPUT_COMMON_HELPERS_FACETS_H
