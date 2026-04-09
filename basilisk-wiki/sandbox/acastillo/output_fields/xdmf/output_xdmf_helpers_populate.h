#ifndef OUTPUT_XDMF_HELPERS_POPULATE_H
#define OUTPUT_XDMF_HELPERS_POPULATE_H

#ifdef HAVE_HDF5

/** ### Populate topo_dset for xdmf (2D topology array [num_cells, 4/8]) */
trace
void populate_topo_dset_xdmf(long **topo_dset, long num_cells, long *offset_cells, hsize_t *count, hsize_t *offset, scalar per_mask, vertex scalar marker) {
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
  long num_cells_iter = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      // Calculate starting index for topo_dset
      long ii = num_cells_iter * count[1];

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
      num_cells_iter++;
    }
  }
}

trace
void populate_topo_dset_slice_xdmf(long **topo_dset, long num_cells, long *offset_cells, hsize_t *count,
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
      long ii = num_cells * count[1];
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

/** ### Populate points_dset_box for xdmf (marker stores global indices) */
void populate_points_dset_box_xdmf(vertex scalar mask, vertex scalar marker, double **points_dset, long num_points, long *offset_points, hsize_t *count, hsize_t *offset) {
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
      long ii = (marker[] - (long)offset[0]) * 3;

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

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_POPULATE_H
