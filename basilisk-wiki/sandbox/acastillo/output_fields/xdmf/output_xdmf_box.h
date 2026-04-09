#ifndef OUTPUT_XDMF_BOX_H
#define OUTPUT_XDMF_BOX_H

/** 
## *output_xmf_box()*: Exports 2D (or 3D) fields within a specified box region.

This function writes light data (XML) and heavy data (HDF5) compatible with XDMF.
The output is constrained to a bounding box defined by `box[2]`.

The arguments and their default values are:

 * **slist** : Pointer to an array of scalar data.
 * **vlist** : Pointer to an array of vector data.
 * **subname** : String used to construct the HDF5 file name.
 * **box** : Array of two coordinates defining the bounding box [min, max].
 * **mode** : Writing mode (HDF5_CONTIGUOUS or HDF5_CHUNKED).
 * **compression_level** : Compression level for GZIP or rate for ZFP.
  
 */

trace void output_xmf_box(scalar *slist, vector *vlist, char *subname, coord box[2], int mode = HDF5_CHUNKED, int compression_level = 6){
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  // Construct the HDF5 file name
  char name[260];  // Buffer for file name construction
  snprintf(name, sizeof(name), "%s.h5", subname); 

  /** Define a scalar field for cell selection with consistent boundaries */ 
  scalar cell_mask[];
  foreach () {
    cell_mask[] = 0.;  // Initialize to 0
#if dimension == 2
    if (x >= box[0].x && x < box[1].x &&
        y >= box[0].y && y < box[1].y)   
#elif dimension == 3
    if (x >= box[0].x && x < box[1].x &&
        y >= box[0].y && y < box[1].y &&
        z >= box[0].z && z < box[1].z)
#endif
    {  
      cell_mask[] = 1.;
    }
  }

  vertex scalar vertex_needed[];
  foreach_vertex(){
    vertex_needed[] = 0;
  }

  foreach (serial, noauto){
    if (cell_mask[] > 0.5){
      vertex_needed[0] = 1;
      vertex_needed[1] = 1;
      vertex_needed[1,1] = 1;
      vertex_needed[0,1] = 1;
#if dimension == 3
      vertex_needed[0,0,1] = 1;
      vertex_needed[1,0,1] = 1;
      vertex_needed[1,1,1] = 1;
      vertex_needed[0,1,1] = 1;
#endif      
    }
  }

  // Obtain the number of points and cells and get a marker to reconstruct the topology
  long num_points = 0, num_cells = 0;
  long num_points_loc = 0, num_cells_loc = 0;
  count_points_and_cells_box(&num_points, &num_cells, &num_points_loc, &num_cells_loc, cell_mask, vertex_needed);

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  // Initialize marker for topology reconstruction
  vertex scalar marker[];
  initialize_marker_box(marker, vertex_needed, offset, 1);

  // Write the light data (XML)
  if (pid() == 0) {
    write_xdmf_light_data(slist, vlist, name, subname, num_cells, num_points);
  }

  // Write the heavy data (HDF5)
  file_id = create_hdf5_file(name);
  if (file_id < 0) return; // Exit if file creation failed

  // Centralized chunk size calculation
  hsize_t chunk_size = compute_chunk_size(num_cells);

  // Populate and write the topology dataset
  long *topo_dset;
  populate_topo_dset_xdmf(&topo_dset, num_cells_loc, offset_cells, count, offset, cell_mask, marker);
  write_dataset(file_id, count, offset, "/Topology", num_cells, num_cells_loc, pow(2, dimension), topo_dset, H5T_NATIVE_LONG, mode, chunk_size, compression_level);
  free(topo_dset);

  // Create group for mesh geometry data
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset_box_xdmf(vertex_needed, marker, &points_dset, num_points_loc, offset_points, count, offset);
  write_dataset(group_id, count, offset, "/Geometry/Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, mode, chunk_size, compression_level);
  free(points_dset);
  H5Gclose(group_id);

  // Create group for cell data
  group_id = H5Gcreate(file_id, "Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Allocate memory and write scalar datasets
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  for (scalar s in slist) {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", s.name);
    populate_scalar_dset(s, scalar_dset, num_cells_loc, offset_cells, count, offset, cell_mask);
    write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, mode, chunk_size, compression_level);
  }
  free(scalar_dset);

  // Allocate memory and write vector datasets
  double *vector_dset = (double *)malloc(num_cells_loc * 3 * sizeof(double));
  for (vector v in vlist) {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", v.x.name);
    populate_vector_dset(v, vector_dset, num_cells_loc, offset_cells, count, offset, cell_mask);
    write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 3, vector_dset, H5T_NATIVE_DOUBLE, mode, chunk_size, compression_level);
  }
  free(vector_dset);
  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
#else
  // HDF5 not available
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_xmf_box() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif // OUTPUT_XDMF_BOX_H
