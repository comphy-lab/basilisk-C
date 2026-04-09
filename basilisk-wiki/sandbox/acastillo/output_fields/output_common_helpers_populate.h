#ifndef OUTPUT_COMMON_HELPERS_POPULATE_H
#define OUTPUT_COMMON_HELPERS_POPULATE_H

#ifdef HAVE_HDF5

/** ### Populate points_dset based on markers and dimensions */
trace
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

trace
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

/** ### Populate scalar_dset using the the scalar s */
trace
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

  long num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Store values
      scalar_dset[num_cells_iter] = s[];
      num_cells_iter++;
    }
  }
}

trace
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

  long num_cells_iter = 0;
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
trace
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

  long num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate starting index
      long ii = num_cells_iter * 3;

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
trace
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

  long num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      long ii = num_cells_iter * 3;
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

#endif // OUTPUT_COMMON_HELPERS_POPULATE_H
