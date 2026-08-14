#ifndef OUTPUT_COMMON_HELPERS_H
#define OUTPUT_COMMON_HELPERS_H

#define shortcut_slice(n, _alpha)                                \
  double alpha = (_alpha - n.x * x - n.y * y - n.z * z) / Delta; \
  if (fabs(alpha) > 0.87)                                        \
    continue;

/** ### Count points and cells in each subdomain and total */
void count_points_and_cells(long *num_points_glob, long *num_cells_glob, long *num_points, long *num_cells, scalar per_mask) {
  foreach_vertex(serial, noauto){
    (*num_points)++;
  }

  foreach (serial, noauto){
    if (per_mask[]){
      (*num_cells)++;
    }
  }

#if _MPI
  long loc_points = *num_points, loc_cells = *num_cells;
  MPI_Allreduce(&loc_points, num_points_glob, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_cells,  num_cells_glob,  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  *num_points_glob = *num_points;
  *num_cells_glob = *num_cells;
#endif
}

void count_points_and_cells_box(long *num_points_glob, long *num_cells_glob, long *num_points, long *num_cells, scalar cell_mask, vertex scalar vertex_needed) {
  foreach_vertex(serial, noauto){
    if (vertex_needed[] > 0.5){
      (*num_points)++;
    }
  }

  foreach (serial, noauto){
    if (cell_mask[] > 0.5){
      (*num_cells)++;
    }
  }

#if _MPI
  long loc_points = *num_points, loc_cells = *num_cells;
  MPI_Allreduce(&loc_points, num_points_glob, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_cells,  num_cells_glob,  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  *num_points_glob = *num_points;
  *num_cells_glob = *num_cells;
#endif
}

void count_points_and_cells_slice(long *num_points_glob, long *num_cells_glob, long *num_points, long *num_cells, scalar per_mask, coord n = {0, 0, 1}, double _alpha = 0) {
  foreach_vertex(serial, noauto){
    shortcut_slice(n, _alpha);
    (*num_points)++;
  }

  foreach (serial, noauto){
    if (per_mask[]){
      (*num_cells)++;
    }
  }

#if _MPI
  long loc_points = *num_points, loc_cells = *num_cells;
  MPI_Allreduce(&loc_points, num_points_glob, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loc_cells,  num_cells_glob,  1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  *num_points_glob = *num_points;
  *num_cells_glob = *num_cells;
#endif
}

#ifdef HAVE_HDF5

/** ### Calculate offsets for points and cells in each subdomain */
void calculate_offsets(long *offset_points, long *offset_cells, long num_points, long num_cells, hsize_t *offset) {
  // Arrays to store the number of points and cells in each subdomain
  long list_points[npe()];
  long list_cells[npe()];

  // Initialize the arrays to zero
  for (int i = 0; i < npe(); ++i){
    list_points[i] = 0;
    list_cells[i] = 0;
  }

  // Set the number of points and cells for the current subdomain
  list_points[pid()] = num_points;
  list_cells[pid()] = num_cells;

#if _MPI
  // Perform an all-reduce operation to gather the number of points and cells from all subdomains
  MPI_Allreduce(list_points, offset_points, npe(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(list_cells, offset_cells, npe(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  // Without MPI, just copy the local values
  for (int i = 0; i < npe(); ++i){
    offset_points[i] = list_points[i];
    offset_cells[i] = list_cells[i];
  }
#endif

  // Calculate the offset for the points in the current subdomain
  offset[0] = 0;
  if (pid() != 0){
    // Sum the offsets of the previous subdomains to get the starting offset for the current subdomain
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_points[i - 1];
    }
  }
}

void calculate_offsets2(long *offset_points, long num_points, hsize_t *offset) {
  // Arrays to store the number of points in each subdomain
  long list_points[npe()];

  // Initialize the arrays to zero
  for (int i = 0; i < npe(); ++i){
    list_points[i] = 0;
  }

  // Set the number of points for the current subdomain
  list_points[pid()] = num_points;

#if _MPI
  // Perform an all-reduce operation to gather the number of points from all subdomains
  MPI_Allreduce(list_points, offset_points, npe(), MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  // Without MPI, just copy the local values
  for (int i = 0; i < npe(); ++i){
    offset_points[i] = list_points[i];
  }
#endif

  // Calculate the offset for the points in the current subdomain
  offset[0] = 0;
  if (pid() != 0){
    // Sum the offsets of the previous subdomains to get the starting offset for the current subdomain
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_points[i - 1];
    }
  }
}

/** ### Initialize marker to rebuild the topology */
void initialize_marker(vertex scalar marker, hsize_t *offset, hsize_t accumulate = 1) {
  long num_points_iter = 0;
  foreach_vertex(serial, noauto){
    marker[] = num_points_iter + offset[0]*accumulate;
    num_points_iter++;
  }
  // The 'dirty' field attribute has been removed
  // marker.dirty = true;
}

void initialize_marker_box(vertex scalar marker, vertex scalar vertex_needed, hsize_t *offset, hsize_t accumulate = 1) {
  long num_points_iter = 0;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    if (vertex_needed[] > 0.5) {
      marker[] = num_points_iter + offset[0]*accumulate;
      num_points_iter++;
    }
  }
}

void initialize_marker_slice(vertex scalar marker, hsize_t *offset, coord n = {0, 0, 1}, double _alpha = 0, hsize_t accumulate = 1) {
  long num_points_iter = 0;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    shortcut_slice(n, _alpha);
    marker[] = num_points_iter + offset[0]*accumulate;
    num_points_iter++;
  }
}

#endif // HAVE_HDF5

#include "acastillo/output_fields/output_common_helpers_data.h"
#include "acastillo/output_fields/output_common_helpers_facets.h"
#include "acastillo/output_fields/output_common_helpers_populate.h"

#endif // OUTPUT_COMMON_HELPERS_H
