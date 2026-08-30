#ifndef OUTPUT_VTKHDF_FACETS_BOX_H
#define OUTPUT_VTKHDF_FACETS_BOX_H

#include "acastillo/output_fields/output_common_helpers_facets.h"
#include "output_vtkhdf_facets_populate.h"

/**
## *output_facets_vtkhdf_box()*: Exports VOF facets, restricted to a box.

Same as `output_facets_vtkhdf()` but only includes interface cells whose
**cell center** lies inside `box[2] = {min, max}`. This is the same
cell-center-in-box inclusion test used by `output_vtkhdf_box()` and
`output_facets_xmf_box()`: a facet whose cell center is inside the box is
included whole, even if part of the reconstructed polygon pokes outside; a
facet whose cell center is outside is dropped whole, even if part of the
polygon pokes in.

The arguments and their default values are:

*c*
: vof scalar field.

*s*
: scalar field to store, for instance, curvature.

*name*
: Output file name generally uses the `.vtkhdf` extension.

*box*
: array of two coordinates defining the bounding box [min, max]. Optional;
  defaults to the whole domain (`[X0,Y0,Z0]` to `[X0+L0,Y0+L0,Z0+L0]`), in
  which case this behaves exactly like `output_facets_vtkhdf()`.

*compression_level*
: Level of compression to use when writing data to the HDF5 file (default=9).

*/

#if dimension == 2
trace void output_facets_vtkhdf_box(scalar c, scalar s, char *name = "domain.vtkhdf", coord box[2] = {{X0, Y0}, {X0 + L0, Y0 + L0}}, int compression_level = 9){
#else
trace void output_facets_vtkhdf_box(scalar c, scalar s, char *name = "domain.vtkhdf", coord box[2] = {{X0, Y0, Z0}, {X0 + L0, Y0 + L0, Z0 + L0}}, int compression_level = 9){
#endif
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hid_t subgroup_id; // HDF5 subgroup ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab
  hsize_t dims[1] = {2};

  // Obtain the number of vertices and facets within the box
  long num_points_loc = 0, num_cells_loc = 0;
  long num_points = 0, num_cells = 0;
  count_vertices_and_facets_box(c, box, &num_points_loc, &num_cells_loc);

  // Connectivity size is identical to num_points_loc because every vertex is unique
  long num_ids_loc = num_points_loc;
  long num_ids = 0;

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()], offset_ids[npe()], offset_offset[npe()];
  calculate_offsets2(offset_offset, num_cells_loc+1,  offset);
  calculate_offsets2(offset_ids,    num_ids_loc,      offset);
  calculate_offsets2(offset_cells,  num_cells_loc,    offset);
  calculate_offsets2(offset_points, num_points_loc,   offset);

  // Global counts
#if _MPI
  MPI_Allreduce(&num_points_loc, &num_points, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&num_cells_loc, &num_cells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&num_ids_loc, &num_ids, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  num_points = num_points_loc;
  num_cells = num_cells_loc;
  num_ids = num_ids_loc;
#endif

  // Centralized chunk size calculation
  hsize_t chunk_size = compute_chunk_size(num_cells);

  // Create a new HDF5 file using helper
  file_id = create_hdf5_file(name);
  if (file_id < 0) return;

  // Create group
  group_id = H5Gcreate(file_id, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create "Version", "Type" (and other) attributes
  dims[0] = 2;
  int version_data[2] = {2, 1};
  create_attribute(group_id, "Version", version_data, dims);
  create_attribute_type(group_id, "Type", "UnstructuredGrid", 16);
  create_attribute_type(group_id, "Description", "Simulation perfomed using Basilisk (/)", 57);

  // Write "NumberOfConnectivityIds", "NumberOfPoints", "NumberOfCells"
  dims[0] = npe();
  write_simple_dataset(group_id, "NumberOfConnectivityIds", offset_ids, dims);
  write_simple_dataset(group_id, "NumberOfPoints", offset_points, dims);
  write_simple_dataset(group_id, "NumberOfCells", offset_cells, dims);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset_facets_vtkhdf_box(c, box, &points_dset, num_points_loc, offset_points, count, offset);
  write_dataset(group_id, count, offset, "Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, HDF5_CHUNKED, compute_chunk_size(num_points), compression_level);
  free(points_dset);

  // Populate and write the types dataset
  char * types_dset;
  populate_types_dset_facets_vtkhdf(&types_dset, num_cells_loc, offset_cells, count, offset);
  write_dataset(group_id, count, offset, "Types", num_cells, num_cells_loc, 1, types_dset, H5T_STD_U8LE, HDF5_CHUNKED, chunk_size, compression_level);
  free(types_dset);

  // Populate and write the connectivity dataset
  long *topo_dset;
  populate_topo_dset_facets_vtkhdf_box(c, box, &topo_dset, num_ids_loc, offset_ids, count, offset);
  write_dataset(group_id, count, offset, "Connectivity", num_ids, num_ids_loc, 1, topo_dset, H5T_NATIVE_LONG, HDF5_CHUNKED, compute_chunk_size(num_ids), compression_level);
  free(topo_dset);

  // Populate and write the offsets dataset
  long *offsets_dset;
  populate_offsets_dset_facets_vtkhdf_box(c, box, &offsets_dset, num_cells_loc, offset_offset, count, offset);
  write_dataset(group_id, count, offset, "Offsets", num_cells+npe(), num_cells_loc+1, 1, offsets_dset, H5T_NATIVE_LONG, HDF5_CHUNKED, chunk_size, compression_level);
  free(offsets_dset);

  // Create subgroup "CellData"
  subgroup_id = H5Gcreate(group_id, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Allocate memory and write scalar datasets
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  populate_scalar_dset_facets_vtkhdf_box(c, box, s, scalar_dset, num_cells_loc, offset_cells, count, offset);
  write_dataset(subgroup_id, count, offset, s.name, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, HDF5_CHUNKED, chunk_size, compression_level);
  free(scalar_dset);

  H5Gclose(subgroup_id);

  // Create subgroup "FieldData"
  subgroup_id = H5Gcreate(group_id, "FieldData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(subgroup_id);

  // Create subgroup "PointData"
  subgroup_id = H5Gcreate(group_id, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(subgroup_id);
  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
#else
  // HDF5 not available
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_facets_vtkhdf_box() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif // OUTPUT_VTKHDF_FACETS_BOX_H
