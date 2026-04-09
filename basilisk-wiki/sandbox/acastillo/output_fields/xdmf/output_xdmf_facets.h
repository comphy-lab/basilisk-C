#ifndef OUTPUT_XDMF_FACETS_H
#define OUTPUT_XDMF_FACETS_H

#include "acastillo/output_fields/output_common_helpers_facets.h"
#include "output_xdmf_facets_populate.h"

/** 
## output_facets_xmf(): Exports the VOF interface.

This function writes a surface from field data using the built-in VOF PLIC
surface. Heavy data is stored using the Hierarchical Data Format HDF5 in parallel.
Because the cell sizes vary, it relies on XDMF's `Mixed` topology.

The arguments and their default values are:

*c*
: vof scalar field.

*s*
: scalar field to store, for instance, curvature.

*subname*
: subname to be used for the output file.

*mode*
: HDF5 writing mode.

*compression_level*
: HDF5 compression level.

*/

trace
void output_facets_xmf(scalar c, scalar s, char *subname, int mode = HDF5_CHUNKED, int compression_level = 6){
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  // Construct the HDF5 file name
  char name[260];  // Buffer for file name construction
  snprintf(name, sizeof(name), "%s.h5", subname); 

  // Obtain the number of vertices and facets
  long num_points_loc = 0, num_cells_loc = 0;
  long num_points = 0, num_cells = 0;
  count_vertices_and_facets(c, &num_points_loc, &num_cells_loc);

  // In Mixed topology, each cell requires `2 + m` integers: `[type, m, id1, ..., id_m]`.
  // The total 1D size is exactly `2 * num_cells_loc + num_points_loc`.
  long topo_size_loc = 2 * num_cells_loc + num_points_loc;
  long topo_size = 0;

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()], offset_topo[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  // Re-use offset logic for topology size
  calculate_offsets2(offset_topo, topo_size_loc, offset);

  // Determine global counts
#if _MPI
  MPI_Allreduce(&num_points_loc, &num_points, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&num_cells_loc, &num_cells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&topo_size_loc, &topo_size, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  num_points = num_points_loc;
  num_cells = num_cells_loc;
  topo_size = topo_size_loc;
#endif

  // Write the light data
  if (pid() == 0) {
    // Custom logic to write XDMF header, mixed topology, and attributes for facets
    char xmf_name[111];
    sprintf(xmf_name, "%s.xmf", subname);
    FILE *fp = fopen(xmf_name, "w");

    write_xdmf_header(fp, name);

    fputs("\t<Domain>\n", fp);
    fputs("\t\t<Grid Name=\"Unstructured Grid\" GridType=\"Uniform\">\n", fp);
    fprintf(fp, "\t\t\t<Time Type=\"Single\" Value=\"%g\" />\n", t);

    // Write Mixed Topology
    fprintf(fp, "\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%ld\">\n", num_cells);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" Dimensions=\"%ld\" DataType=\"Int\" Precision=\"8\" >\n", topo_size);
    fputs("\t\t\t\t\t&HeavyData;/Topology\n", fp);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Topology>\n", fp);

    // Write Geometry
    fputs("\t\t\t<Geometry GeometryType=\"XYZ\">\n", fp);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Dimensions=\"%ld 3\" Precision=\"8\" >\n", num_points);
    fputs("\t\t\t\t\t&HeavyData;/Geometry/Points\n", fp);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Geometry>\n", fp);

    // Write Attribute
    fprintf(fp, "\t\t\t<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
    fprintf(fp, "\t\t\t\t<DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", num_cells);
    fprintf(fp, "\t\t\t\t\t&HeavyData;/Cells/%s\n", s.name);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Attribute>\n", fp);

    write_xdmf_footer(fp);
    fflush(fp);
    fclose(fp);
  }

  // Write the heavy data
  file_id = create_hdf5_file(name);
  if (file_id < 0) return; // Exit if file creation failed

  // Centralized chunk size calculation
  hsize_t chunk_size = compute_chunk_size(topo_size);

  // Populate and write the topology dataset
  long *topo_dset;
  populate_topo_dset_facets_xdmf(c, &topo_dset, topo_size_loc, offset_topo, offset_points, count, offset);
  write_dataset(file_id, count, offset, "/Topology", topo_size, topo_size_loc, 1, topo_dset, H5T_NATIVE_LONG, mode, chunk_size, compression_level);
  free(topo_dset);

  // Create group for mesh geometry data
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset_facets_xdmf(c, &points_dset, num_points_loc, offset_points, count, offset);
  write_dataset(group_id, count, offset, "/Geometry/Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, mode, compute_chunk_size(num_points), compression_level);
  free(points_dset);
  H5Gclose(group_id);

  // Create group for cell data
  group_id = H5Gcreate(file_id, "Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Allocate memory and write scalar datasets
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  char substamp[1024];
  sprintf(substamp, "/Cells/%s", s.name);
  populate_scalar_dset_facets_xdmf(c, s, scalar_dset, num_cells_loc, offset_cells, count, offset);
  write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, mode, compute_chunk_size(num_cells), compression_level);
  free(scalar_dset);

  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
#else
  // HDF5 not available - print warning and return
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_facets_xmf() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif // OUTPUT_XDMF_FACETS_H
