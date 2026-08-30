#ifndef OUTPUT_XDMF_FACETS_BOX_H
#define OUTPUT_XDMF_FACETS_BOX_H

#include "acastillo/output_fields/output_common_helpers_facets.h"
#include "output_xdmf_facets_populate.h"

/**
## output_facets_xmf_box(): Exports the VOF interface, restricted to a box.

Same as `output_facets_xmf_list()` (multiple scalar fields on the same VOF
interface mesh, written as separate HDF5 `/Cells/<name>` datasets / XDMF
Attributes) but only includes interface cells whose **cell center** lies
inside `box[2] = {min, max}`. This is the same semantics `output_xmf_box()`
uses for its regular-field export: a per-cell inclusion test, not a
geometric clip of the reconstructed facet polygon — a facet whose cell
center is inside the box is included whole, even if part of the polygon
pokes outside; a facet whose cell center is outside is dropped whole, even
if part of the polygon pokes in.

The arguments and their default values are:

*c*
: vof scalar field.

*slist*
: list of scalar fields to store on the facets (e.g. curvature, pressure, etc.).

*subname*
: subname to be used for the output file.

*box*
: array of two coordinates defining the bounding box [min, max]. Optional;
  defaults to the whole domain (`[X0,Y0,Z0]` to `[X0+L0,Y0+L0,Z0+L0]`), in
  which case this behaves exactly like `output_facets_xmf_list()`.

*mode*
: HDF5 writing mode.

*compression_level*
: HDF5 compression level.

*/

#if dimension == 2
trace
void output_facets_xmf_box(scalar c, scalar * slist, char *subname, coord box[2] = {{X0, Y0}, {X0 + L0, Y0 + L0}}, int mode = HDF5_CHUNKED, int compression_level = 6){
#else
trace
void output_facets_xmf_box(scalar c, scalar * slist, char *subname, coord box[2] = {{X0, Y0, Z0}, {X0 + L0, Y0 + L0, Z0 + L0}}, int mode = HDF5_CHUNKED, int compression_level = 6){
#endif
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  // Construct the HDF5 file name
  char name[260];
  snprintf(name, sizeof(name), "%s.h5", subname);

  // Obtain the number of vertices and facets within the box
  long num_points_loc = 0, num_cells_loc = 0;
  long num_points = 0, num_cells = 0;
  count_vertices_and_facets_box(c, box, &num_points_loc, &num_cells_loc);

  long topo_size_loc = 2 * num_cells_loc + num_points_loc;
  long topo_size = 0;

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()], offset_topo[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);
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
    char xmf_name[111];
    sprintf(xmf_name, "%s.xmf", subname);
    FILE *fp = fopen(xmf_name, "w");

    write_xdmf_header(fp, name);

    fputs("\t<Domain>\n", fp);
    fputs("\t\t<Grid Name=\"Unstructured Grid\" GridType=\"Uniform\">\n", fp);
    fprintf(fp, "\t\t\t<Time Type=\"Single\" Value=\"%g\" />\n", t);

    fprintf(fp, "\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%ld\">\n", num_cells);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" Dimensions=\"%ld\" DataType=\"Int\" Precision=\"8\" >\n", topo_size);
    fputs("\t\t\t\t\t&HeavyData;/Topology\n", fp);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Topology>\n", fp);

    fputs("\t\t\t<Geometry GeometryType=\"XYZ\">\n", fp);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Dimensions=\"%ld 3\" Precision=\"8\" >\n", num_points);
    fputs("\t\t\t\t\t&HeavyData;/Geometry/Points\n", fp);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Geometry>\n", fp);

    for (scalar s in slist) {
      fprintf(fp, "\t\t\t<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
      fprintf(fp, "\t\t\t\t<DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", num_cells);
      fprintf(fp, "\t\t\t\t\t&HeavyData;/Cells/%s\n", s.name);
      fputs("\t\t\t\t</DataItem>\n", fp);
      fputs("\t\t\t</Attribute>\n", fp);
    }

    write_xdmf_footer(fp);
    fflush(fp);
    fclose(fp);
  }

  // Write the heavy data
  file_id = create_hdf5_file(name);
  if (file_id < 0) return;

  hsize_t chunk_size = compute_chunk_size(topo_size);

  // Topology
  long *topo_dset;
  populate_topo_dset_facets_xdmf_box(c, box, &topo_dset, topo_size_loc, offset_topo, offset_points, count, offset);
  write_dataset(file_id, count, offset, "/Topology", topo_size, topo_size_loc, 1, topo_dset, H5T_NATIVE_LONG, mode, chunk_size, compression_level);
  free(topo_dset);

  // Geometry
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  double *points_dset;
  populate_points_dset_facets_xdmf_box(c, box, &points_dset, num_points_loc, offset_points, count, offset);
  write_dataset(group_id, count, offset, "/Geometry/Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, mode, compute_chunk_size(num_points), compression_level);
  free(points_dset);
  H5Gclose(group_id);

  // Cell data, one dataset per field in slist
  group_id = H5Gcreate(file_id, "Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  hsize_t chunk_size_cells = compute_chunk_size(num_cells);
  for (scalar s in slist) {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", s.name);
    populate_scalar_dset_facets_xdmf_box(c, box, s, scalar_dset, num_cells_loc, offset_cells, count, offset);
    write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, mode, chunk_size_cells, compression_level);
  }
  free(scalar_dset);
  H5Gclose(group_id);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
#else
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_facets_xmf_box() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif // OUTPUT_XDMF_FACETS_BOX_H
