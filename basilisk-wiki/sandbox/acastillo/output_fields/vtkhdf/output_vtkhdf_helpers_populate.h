#ifndef OUTPUT_VTKHDF_HELPERS_POPULATE_H
#define OUTPUT_VTKHDF_HELPERS_POPULATE_H

#ifdef HAVE_HDF5
/** ### Populate points_dset based on markers and dimensions */ 
void populate_points_dset(double **points_dset, long num_points, long *offset_points, hsize_t *count, hsize_t *offset) {
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

  // Iterate over each vertex
  long num_points_iter = 0;
  foreach_vertex(serial, noauto){
    // Calculate starting index
    long ii = num_points_iter * 3;

    // Store coordinates
    (*points_dset)[ii + 0] = x;
    (*points_dset)[ii + 1] = y;
    #if dimension == 2
      (*points_dset)[ii + 2] = 0.;
    #else 
      (*points_dset)[ii + 2] = z;
    #endif
    num_points_iter++;
  }
}

void populate_points_dset_box(vertex scalar mask, vertex scalar marker, double **points_dset, long num_points, long *offset_points, hsize_t *count, hsize_t *offset) {
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

  // Iterate over each vertex
  foreach_vertex(serial, noauto){
    if (mask[] >= 0.5){
      long ii = marker[] * 3;

      // Store coordinates
      (*points_dset)[ii + 0] = x;
      (*points_dset)[ii + 1] = y;
      #if dimension == 2
        (*points_dset)[ii + 2] = 0.;
      #else 
        (*points_dset)[ii + 2] = z;
      #endif
    }
  }
}

void populate_points_dset_slice(double **points_dset, long num_points, long *offset_points, hsize_t *count,
                                hsize_t *offset, coord n = {0, 0, 1}, double _alpha = 0)
{
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

  // Iterate over each vertex
  num_points = 0;
  foreach_vertex(serial, noauto){
    shortcut_slice(n, _alpha);

    // Calculate starting index
    long ii = num_points * 3;

    // Store coordinates
    (*points_dset)[ii + 0] = x;
    (*points_dset)[ii + 1] = y;
    (*points_dset)[ii + 2] = z;
    num_points++;
  }
}

/** ### Populate types_dset  */ 
void populate_types_dset(char **types_dset, char type, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset) {
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

  for (long i = 0; i < num_cells; ++i){
    (*types_dset)[i] = type;      
  }
}

/** ### Populate scalar_dset using the the scalar s */ 
void populate_scalar_dset(scalar s, double *scalar_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask) {
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

  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      scalar_dset[ii] = s[];
      ii++;
    }
  }
}

void populate_scalar_dset_slice(scalar s, double *scalar_dset, long num_cells, long *offset_cells, hsize_t *count,
                                hsize_t *offset, scalar per_mask, coord n = {0, 0, 1}, double _alpha = 0)
{
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

  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      if (n.x == 1)
        scalar_dset[ii] = 0.5 * (val(s) + val(s, 1, 0, 0));
      else if (n.y == 1)
        scalar_dset[ii] = 0.5 * (val(s) + val(s, 0, 1, 0));
      else
        scalar_dset[ii] = 0.5 * (val(s) + val(s, 0, 0, 1));
      ii++;
    }
  }
}

/** ### Populate vector_dset using the vector v */ 
void populate_vector_dset(vector v, double *vector_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask) {
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_cells;
  count[1] = 3;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_cells[i - 1];
    }
  }

  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate starting index
      long ij = ii * 3;

      // Store each component
      vector_dset[ij + 0] = v.x[];
      vector_dset[ij + 1] = v.y[];
      #if dimension == 2
        vector_dset[ij + 2] = 0.;
      #else
        vector_dset[ij + 2] = v.z[];
      #endif
      ii++;
    }
  }
}

#if dimension == 3
void populate_vector_dset_slice(vector v, double *vector_dset, long num_cells, long *offset_cells, hsize_t *count,
                                hsize_t *offset, scalar per_mask, coord n = {0, 0, 1}, double _alpha = 0){
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_cells;
  count[1] = 3;
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_cells[i - 1];
    }
  }

  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      long ij = ii * 3;
      if (n.x == 1){
        vector_dset[ij + 0] = 0.5 * (val(v.x) + val(v.x, 1, 0, 0));
        vector_dset[ij + 1] = 0.5 * (val(v.y) + val(v.y, 1, 0, 0));
        vector_dset[ij + 2] = 0.5 * (val(v.z) + val(v.z, 1, 0, 0));
      }
      else if (n.y == 1){
        vector_dset[ij + 0] = 0.5 * (val(v.x) + val(v.x, 0, 1, 0));
        vector_dset[ij + 1] = 0.5 * (val(v.y) + val(v.y, 0, 1, 0));
        vector_dset[ij + 2] = 0.5 * (val(v.z) + val(v.z, 0, 1, 0));
      }
      else{
        vector_dset[ij + 0] = 0.5 * (val(v.x) + val(v.x, 0, 0, 1));
        vector_dset[ij + 1] = 0.5 * (val(v.y) + val(v.y, 0, 0, 1));
        vector_dset[ij + 2] = 0.5 * (val(v.z) + val(v.z, 0, 0, 1));
      }
      ii++;
    }
  }
}
#endif

/** ### Populate offsets_dset  */ 
void populate_offsets_dset(long **offsets_dset, char noffset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset) {
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

  // Allocate memory for topo_dset
  *offsets_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  for (long i = 0; i < num_cells; ++i){   
    (*offsets_dset)[i] = (long)i * (long)noffset;      
  }
}

/** ### Populate topo_dset based on markers and dimensions */ 
void populate_topo_dset(long **topo_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask, vertex scalar marker) {
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_cells;
  count[1] = pow(2, dimension);
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_cells[i - 1];
    }
  }

  // Allocate memory for topo_dset
  *topo_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  // Iterate over each cell
  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate starting index for topo_dset
      long ij = ii * count[1];

      // Assign marker values to topo_dset
      (*topo_dset)[ij + 0] = (long)marker[];
      (*topo_dset)[ij + 1] = (long)marker[1, 0];
      (*topo_dset)[ij + 2] = (long)marker[1, 1];
      (*topo_dset)[ij + 3] = (long)marker[0, 1];

      #if dimension == 3
        // Additional assignments for 3D
        (*topo_dset)[ij + 4] = (long)marker[0, 0, 1];
        (*topo_dset)[ij + 5] = (long)marker[1, 0, 1];
        (*topo_dset)[ij + 6] = (long)marker[1, 1, 1];
        (*topo_dset)[ij + 7] = (long)marker[0, 1, 1];
      #endif
      ii++;
    }
  }
  count[0] = num_cells*(long)pow(2, dimension);
  count[1] = 1;
}

void populate_topo_dset_slice(long **topo_dset, long num_cells, long *offset_cells, hsize_t *count,
                              hsize_t *offset, scalar per_mask, vertex scalar marker, coord n = {0, 0, 1}, double _alpha = 0)
{
  // Each process defines dataset in memory and writes to an hyperslab
  count[0] = num_cells;
  count[1] = pow(2, dimension - 1);
  offset[0] = 0;
  offset[1] = 0;
  if (pid() != 0){
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_cells[i - 1];
    }
  }

  // Allocate memory for topo_dset
  *topo_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  // Iterate over each cell
  long ii = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate index
      long ij = ii * count[1];
      if (n.x == 1){
        (*topo_dset)[ij + 0] = (long)marker[1, 0, 0];
        (*topo_dset)[ij + 1] = (long)marker[1, 1, 0];
        (*topo_dset)[ij + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ij + 3] = (long)marker[1, 0, 1];
      }
      else if (n.y == 1){
        (*topo_dset)[ij + 0] = (long)marker[0, 1, 0];
        (*topo_dset)[ij + 1] = (long)marker[1, 1, 0];
        (*topo_dset)[ij + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ij + 3] = (long)marker[0, 1, 1];
      }
      else{
        (*topo_dset)[ij + 0] = (long)marker[0, 0, 1];
        (*topo_dset)[ij + 1] = (long)marker[1, 0, 1];
        (*topo_dset)[ij + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ij + 3] = (long)marker[0, 1, 1];
      }
      ii++;
    }
  }
  count[0] = ii*(long)pow(2, dimension - 1);
  count[1] = 1;
}
#endif // HAVE_HDF5

#endif // OUTPUT_VTKHDF_HELPERS_POPULATE_H
