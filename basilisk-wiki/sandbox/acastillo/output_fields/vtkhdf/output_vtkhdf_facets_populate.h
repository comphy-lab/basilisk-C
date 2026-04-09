#ifndef OUTPUT_VTKHDF_FACETS_POPULATE_H
#define OUTPUT_VTKHDF_FACETS_POPULATE_H

#ifdef HAVE_HDF5

/** ### Populate points_dset for facets in vtkhdf (interleaved) */
trace
void populate_points_dset_facets_vtkhdf(scalar c, double **points_dset, long num_points, long *offset_points, hsize_t *count, hsize_t *offset) {
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

/** ### Populate types_dset for facets */
void populate_types_dset_facets_vtkhdf(char **types_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset) {
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

  // Allocate memory for types_dset
  *types_dset = (char *)malloc(count[0] * count[1] * sizeof(char));

  #if dimension == 2
    char type = 4; // VTK_POLY_LINE
  #else
    char type = 7; // VTK_POLYGON
  #endif

  // All cells are of the same polygon type
  for (long i = 0; i < num_cells; ++i){
    (*types_dset)[i] = type;
  }
}

/** ### Populate offsets_dset for facets (cumulative points per cell) */
void populate_offsets_dset_facets_vtkhdf(scalar c, long **offsets_dset, long num_cells_loc, long *offset_offset, hsize_t *count, hsize_t *offset) {
  // For offsets, the array has length (num_cells + npe()) because each core writes an extra element?
  // Wait, in VTKHDF offset arrays have length N+1 per partition!
  // BUT the total length is N + npe().
  count[0] = num_cells_loc + 1; // write N+1 values per process
  count[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_offset[i - 1]; // Offset in the global arrays 
    }
  }

  *offsets_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  long current_offset = 0;
  long ifacet = 0;
  
  // The first offset is ALWAYS 0 for this partition in VTKHDF because point IDs are zero-indexed per partition!
  (*offsets_dset)[0] = 0;
  ifacet++;

  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets 
      if (m > 0){
        current_offset += m;
        (*offsets_dset)[ifacet] = current_offset;
        ifacet++;
      }
    }
  }
}

/** ### Populate topo_dset (Connectivity) for facets */
trace
void populate_topo_dset_facets_vtkhdf(scalar c, long **topo_dset, long num_points_loc, long *offset_points, hsize_t *count, hsize_t *offset) {
  // Connectivity length is exactly the number of points (since each point is unique)
  count[0] = num_points_loc;
  count[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_points[i - 1]; // using points array offset because size is same
    }
  }

  *topo_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

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
        for (int i = 0; i < m; i++){
          // Local vertex index
          (*topo_dset)[iverts] = iverts;
          iverts++;
        }
      }
    }
  }
}

/** ### Populate scalar_dset for facets */
trace
void populate_scalar_dset_facets_vtkhdf(scalar c, scalar s, double *scalar_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset) {
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
      shortcut_facets 
      if (m > 0){
        scalar_dset[ifacet] = s[];
        ifacet++;
      }
    }
  }
}

#endif // HAVE_HDF5
#endif // OUTPUT_VTKHDF_FACETS_POPULATE_H
