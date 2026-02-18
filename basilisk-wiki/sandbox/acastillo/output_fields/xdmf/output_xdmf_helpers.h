/** 
# Helper functions for output_xdmf.h 
*/ 

#ifndef OUTPUT_XDMF_HELPERS_H
#define OUTPUT_XDMF_HELPERS_H

#define shortcut_slice(n, _alpha)                                \
  double alpha = (_alpha - n.x * x - n.y * y - n.z * z) / Delta; \
  if (fabs(alpha) > 0.87)                                        \
    continue;

/** ### Count points and cells in each subdomain and total */ 
void count_points_and_cells(int *num_points_glob, int *num_cells_glob, int *num_points, int *num_cells, scalar per_mask) {
  foreach_vertex(serial, noauto){
    (*num_points)++;
  }

  foreach (serial, noauto){
    if (per_mask[]){
      (*num_cells)++;
    }
  }

#if _MPI
  MPI_Allreduce(num_points, num_points_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_cells, num_cells_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  *num_points_glob = *num_points;
  *num_cells_glob = *num_cells;
#endif
}

void count_points_and_cells_slice(int *num_points_glob, int *num_cells_glob, int *num_points, int *num_cells, scalar per_mask, coord n = {0, 0, 1}, double _alpha = 0) {
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
  MPI_Allreduce(num_points, num_points_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_cells, num_cells_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  *num_points_glob = *num_points;
  *num_cells_glob = *num_cells;
#endif
}

#ifdef HAVE_HDF5

/** ### Calculate offsets for points and cells in each subdomain */ 
void calculate_offsets(int *offset_points, int *offset_cells, int num_points, int num_cells, hsize_t *offset) {
  // Arrays to store the number of points and cells in each subdomain
  int list_points[npe()];
  int list_cells[npe()];

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
  MPI_Allreduce(list_points, offset_points, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(list_cells, offset_cells, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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

/** ### Initialize marker to rebuild the topology */ 
void initialize_marker(vertex scalar marker, hsize_t *offset) {
  int num_points = 0;
  foreach_vertex(serial, noauto){    
    marker[] = num_points + offset[0];
    num_points++;
  }
  marker.dirty = true;
}

void initialize_marker_slice(vertex scalar marker, hsize_t *offset, coord n = {0, 0, 1}, double _alpha = 0) {
  int num_points = 0;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    shortcut_slice(n, _alpha);
    marker[] = num_points + offset[0];
    num_points++;
  }
}

#endif // HAVE_HDF5

#include "output_xdmf_helpers_xml.h"
#include "output_xdmf_helpers_data.h"
#include "output_xdmf_helpers_populate.h"

#endif // OUTPUT_XDMF_HELPERS_H