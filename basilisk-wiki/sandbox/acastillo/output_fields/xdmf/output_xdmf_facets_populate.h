#ifndef OUTPUT_XDMF_FACETS_POPULATE_H
#define OUTPUT_XDMF_FACETS_POPULATE_H

#ifdef HAVE_HDF5

/** ### Populate points_dset for facets in xdmf (interleaved) */
trace
void populate_points_dset_facets_xdmf(scalar c, double **points_dset, long num_points, long *offset_points, hsize_t *count, hsize_t *offset) {
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_points;
  count[1] = 3;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_points[i - 1];
    }
  }

  // Allocate memory for points_dset
  *points_dset = (double *)malloc(count[0] * count[1] * sizeof(double));

  long iverts = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets // we cycle if cell is not at the interface
      coord _p = {x, y, z};
      for (int i = 0; i < m; i++){
        long ii = iverts * 3;
        (*points_dset)[ii + 0] = _p.x + v[i].x * Delta;
        (*points_dset)[ii + 1] = _p.y + v[i].y * Delta;
        #if dimension == 2
          (*points_dset)[ii + 2] = 0.;
        #else
          (*points_dset)[ii + 2] = _p.z + v[i].z * Delta;
        #endif
        iverts++;
      }
    }
  }
}

/** ### Populate topo_dset for facets in xdmf (Mixed Topology) */
trace
void populate_topo_dset_facets_xdmf(scalar c, long **topo_dset, long topo_size, long *offset_topo, long *offset_points, hsize_t *count, hsize_t *offset) {
  
  // Hyperslab parameters for the 1D topology array
  count[0] = topo_size;
  count[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_topo[i - 1];
    }
  }

  // Allocate memory for topo_dset
  *topo_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  // Determine XDMF type
  #if dimension == 2
    long type = 2; // Polyline
  #else
    long type = 3; // Polygon
  #endif

  // Offset points determines the starting node ID for this process
  long node_offset = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      node_offset += offset_points[i - 1];
    }
  }

  long idata = 0;
  long iverts = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets 
      if (m > 0) {
        (*topo_dset)[idata++] = type;
        (*topo_dset)[idata++] = m;
        for (int i = 0; i < m; i++){
          (*topo_dset)[idata++] = node_offset + iverts;
          iverts++;
        }
      }
    }
  }
}

/** ### Populate scalar_dset for facets */
trace
void populate_scalar_dset_facets_xdmf(scalar c, scalar s, double *scalar_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset) {
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_cells;
  count[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_cells[i - 1];
    }
  }

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
        scalar_dset[ifacet] = s[];
        ifacet++;
      }
    }
  }
}

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_FACETS_POPULATE_H
