#ifndef OUTPUT_COMMON_HELPERS_FACETS_H
#define OUTPUT_COMMON_HELPERS_FACETS_H

#include "geometry.h"
#include "fractions.h"

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

### Count the number of vertices and facets */

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

#endif // OUTPUT_COMMON_HELPERS_FACETS_H
