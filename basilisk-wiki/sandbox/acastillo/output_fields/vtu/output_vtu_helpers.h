#ifndef OUTPUT_VTU_HELPERS_H
#define OUTPUT_VTU_HELPERS_H

/**
# Helper functions for output_vtu.h

This file acts as a bundle header for all VTU-related helper functions.
It defines the core macros (`shortcut_facets`, `shortcut_slice`) and geometry 
functions (`count_vertices_and_facets`, `populate_vertex_and_offset_arrays`)
that are required for VOF reconstruction. 

It also includes helper modules for writing XML and binary data.
*/

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

/** Macro to simplify dealing with slices */ 
#define shortcut_slice(n, _alpha)                                \
  double alpha = (_alpha - n.x * x - n.y * y - n.z * z) / Delta; \
  if (fabs(alpha) > 0.87)                                        \
    continue;

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

/** ### Populate vertex coordinates and facet offsets */
void populate_vertex_and_offset_arrays(scalar c, long nverts, long nfacets, double *pt_array_x, double *pt_array_y, double *pt_array_z, long *offsets) {
  long iverts = 0, ifacet = 0, offset = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets // we cycle if cell is not at the interface

      // Calculate and store vertex coordinates
      coord _p = {x, y, z};
      for (int i = 0; i < m; i++){
        pt_array_x[iverts] = _p.x + v[i].x * Delta;
        pt_array_y[iverts] = _p.y + v[i].y * Delta;
        #if dimension == 3
          pt_array_z[iverts] = _p.z + v[i].z * Delta;
        #endif
        iverts++;
      }
      // Store facet offset if there are vertices in the facet
      if (m > 0){
        offset += m;
        offsets[ifacet] = offset;
        ifacet++;
      }
    }
  }
}

/** ### Populate arrays with values at the facets */
void populate_facet_arrays(scalar c, scalar s, long nverts, long nfacets, double *val_array_s) {
  long ifacet = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets // we cycle if cell is not at the interface

      if (m > 0){
        val_array_s[ifacet] = s[];
        ifacet++;
      }
    }
  }
}


/**
## XML Formatting Helpers 
Includes functions for writing VTU/PVTU headers and XML metadata (light data).
*/
#include "output_vtu_helpers_xml.h"

/**
## Binary Data Writing Helpers
Includes functions for iterating over the grid and writing heavy binary data.
*/
#include "output_vtu_helpers_data.h"

#endif
