#ifndef OUTPUT_XDMF_HELPERS_POPULATE_H
#define OUTPUT_XDMF_HELPERS_POPULATE_H

#ifdef HAVE_HDF5

/** ### Populate topo_dset based on markers and dimensions */ 
void populate_topo_dset(long **topo_dset, int num_cells, int *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask, vertex scalar marker) {
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
  int num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){    
      // Calculate starting index for topo_dset
      int ii = num_cells_iter * count[1];

      // Assign marker values to topo_dset
      (*topo_dset)[ii + 0] = (long)marker[];
      (*topo_dset)[ii + 1] = (long)marker[1, 0];
      (*topo_dset)[ii + 2] = (long)marker[1, 1];
      (*topo_dset)[ii + 3] = (long)marker[0, 1];

      #if dimension == 3
        // Additional assignments for 3D
        (*topo_dset)[ii + 4] = (long)marker[0, 0, 1];
        (*topo_dset)[ii + 5] = (long)marker[1, 0, 1];
        (*topo_dset)[ii + 6] = (long)marker[1, 1, 1];
        (*topo_dset)[ii + 7] = (long)marker[0, 1, 1];
      #endif
    }
    num_cells_iter++;
  }
}

void populate_topo_dset_slice(long **topo_dset, int num_cells, int *offset_cells, hsize_t *count,
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
  num_cells = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate index
      int ii = num_cells * count[1];
      if (n.x == 1){
        (*topo_dset)[ii + 0] = (long)marker[1, 0, 0];
        (*topo_dset)[ii + 1] = (long)marker[1, 1, 0];
        (*topo_dset)[ii + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ii + 3] = (long)marker[1, 0, 1];
      }
      else if (n.y == 1){
        (*topo_dset)[ii + 0] = (long)marker[0, 1, 0];
        (*topo_dset)[ii + 1] = (long)marker[1, 1, 0];
        (*topo_dset)[ii + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ii + 3] = (long)marker[0, 1, 1];
      }
      else{
        (*topo_dset)[ii + 0] = (long)marker[0, 0, 1];
        (*topo_dset)[ii + 1] = (long)marker[1, 0, 1];
        (*topo_dset)[ii + 2] = (long)marker[1, 1, 1];
        (*topo_dset)[ii + 3] = (long)marker[0, 1, 1];
      }
      num_cells++;
    }
  }
}

/** ### Populate points_dset based on markers and dimensions */ 
void populate_points_dset(double **points_dset, int num_points, int *offset_points, hsize_t *count, hsize_t *offset) {
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
  int num_points_iter = 0;
  foreach_vertex(serial, noauto){
    // Calculate starting index
    int ii = num_points_iter * 3;

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

void populate_points_dset_slice(double **points_dset, int num_points, int *offset_points, hsize_t *count,
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
    int ii = num_points * 3;

    // Store coordinates
    (*points_dset)[ii + 0] = x;
    (*points_dset)[ii + 1] = y;
    (*points_dset)[ii + 2] = z;
    num_points++;
  }
}

/** ### Populate scalar_dset using the the scalar s */ 
void populate_scalar_dset(scalar s, double *scalar_dset, int num_cells, int *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask) {
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

  int num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Store values
      scalar_dset[num_cells_iter] = s[];
      num_cells_iter++;
    }
  }
}

void populate_scalar_dset_slice(scalar s, double *scalar_dset, int num_cells, int *offset_cells, hsize_t *count,
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

  int num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      if (n.x == 1)
        scalar_dset[num_cells_iter] = 0.5 * (val(s) + val(s, 1, 0, 0));
      else if (n.y == 1)
        scalar_dset[num_cells_iter] = 0.5 * (val(s) + val(s, 0, 1, 0));
      else
        scalar_dset[num_cells_iter] = 0.5 * (val(s) + val(s, 0, 0, 1));
      num_cells_iter++;
    }
  }
}

/** ### Populate vector_dset using the vector v */ 
void populate_vector_dset(vector v, double *vector_dset, int num_cells, int *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask) {
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

  int num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate starting index
      int ii = num_cells_iter * 3;

      // Store each component
      vector_dset[ii + 0] = v.x[];
      vector_dset[ii + 1] = v.y[];
      #if dimension == 2
        vector_dset[ii + 2] = 0.;
      #else
        vector_dset[ii + 2] = v.z[];
      #endif
      num_cells_iter++;
    }
  }
}

#if dimension == 3
void populate_vector_dset_slice(vector v, double *vector_dset, int num_cells, int *offset_cells, hsize_t *count,
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

  int num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      int ii = num_cells_iter * 3;
      if (n.x == 1){
        vector_dset[ii + 0] = 0.5 * (val(v.x) + val(v.x, 1, 0, 0));
        vector_dset[ii + 1] = 0.5 * (val(v.y) + val(v.y, 1, 0, 0));
        vector_dset[ii + 2] = 0.5 * (val(v.z) + val(v.z, 1, 0, 0));
      }
      else if (n.y == 1){
        vector_dset[ii + 0] = 0.5 * (val(v.x) + val(v.x, 0, 1, 0));
        vector_dset[ii + 1] = 0.5 * (val(v.y) + val(v.y, 0, 1, 0));
        vector_dset[ii + 2] = 0.5 * (val(v.z) + val(v.z, 0, 1, 0));
      }
      else{
        vector_dset[ii + 0] = 0.5 * (val(v.x) + val(v.x, 0, 0, 1));
        vector_dset[ii + 1] = 0.5 * (val(v.y) + val(v.y, 0, 0, 1));
        vector_dset[ii + 2] = 0.5 * (val(v.z) + val(v.z, 0, 0, 1));
      }
      num_cells_iter++;
    }
  }
}
#endif

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_POPULATE_H
