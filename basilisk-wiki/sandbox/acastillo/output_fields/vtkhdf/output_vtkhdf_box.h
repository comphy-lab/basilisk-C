#ifndef OUTPUT_VTKHDF_BOX_H
#define OUTPUT_VTKHDF_BOX_H

/** 
## *output_vtkhdf_box()*: Exports 2D (or 3D) fields within a specified box region.

This function writes one [VTKHDF
file](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html) which
can be read using Paraview. The VTKHDF file format is a file format relying on
HDF5. It is meant to provide good I/O performance as well as robust and flexible
parallel I/O capabilities. The file stores scalar and vector fields defined at
the center points are stored at the cell center of an unstructured grid with
type `VTK_QUAD` in 2D, or `VTK_HEXAHEDRON` in 3D. 

The arguments and their default values are:

*list*
: pointer to the list of scalar fields to be exported.

*vlist*
: pointer to the list of vector fields to be exported.

*name*
: Output file name generally uses the `.vtkhdf` extension.

*box*
: array of two coordinates defining the bounding box [min, max].

*compression_level* 
: Level of compression to use when writing data to the HDF5 file (default=9).
 
### Example Usage

```c
scalar * list = {a,b};
vector * vlist = {c,d};
coord box[2] = {{-0.25, -0.25, -0.025}, {0.25, 0.25, 0.025}};
output_vtkhdf_box(list, vlist, "domain.vtkhdf", box);
```

*/

trace void output_vtkhdf_box(scalar *list, vector *vlist, char *name = "domain.vtkhdf", coord box[2], int compression_level = 9){
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hid_t subgroup_id; // HDF5 subgroup ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab
  hsize_t dims[1] = {2};
  
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

  // VTK cell types: VTK_QUAD (in 2D) or VTK_HEXAHEDRON (in 3D)
  int type, noffset;
  #if dimension == 2
    type = 9;
    noffset = 4;
  #elif dimension == 3
    type = 12;
    noffset = 8;
  #endif

  // Obtain the number of points and cells and get a marker to reconstruct the topology
  long num_points = 0, num_points_loc = 0;
  long num_cells = 0,  num_cells_loc = 0;
  
  count_points_and_cells_box(&num_points, &num_cells, &num_points_loc, &num_cells_loc, cell_mask, vertex_needed);
  
  long num_ids = num_cells*noffset;
  long num_ids_loc = num_cells_loc*noffset;

  // Centralized chunk size calculation
  hsize_t chunk_size = compute_chunk_size(num_cells);

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()], offset_ids[npe()], offset_offset[npe()];
  calculate_offsets2(offset_offset, num_cells_loc+1,  offset);
  calculate_offsets2(offset_ids,    num_ids_loc,      offset);
  calculate_offsets2(offset_cells,  num_cells_loc,    offset);
  calculate_offsets2(offset_points, num_points_loc,   offset);

  // Initialize marker for topology reconstruction
  vertex scalar marker[];
  initialize_marker_box(marker, vertex_needed, offset, 0);
  
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
  populate_points_dset_box_vtkhdf(vertex_needed, marker, &points_dset, num_points_loc, offset_points, count, offset);
  write_dataset(group_id, count, offset, "Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, HDF5_CHUNKED, chunk_size, compression_level);
  free(points_dset);

  // Populate and write the types dataset
  char * types_dset;
  populate_types_dset(&types_dset, type, num_cells_loc, offset_cells, count, offset);
  write_dataset(group_id, count, offset, "Types", num_cells, num_cells_loc, 1, types_dset, H5T_STD_U8LE, HDF5_CHUNKED, chunk_size, compression_level);
  free(types_dset);  

  // Populate and write the connectivity dataset
  long *topo_dset;
  populate_topo_dset_vtkhdf(&topo_dset, num_cells_loc, offset_ids, count, offset, cell_mask, marker);
  write_dataset(group_id, count, offset, "Connectivity", num_ids, num_ids_loc, 1, topo_dset, H5T_NATIVE_LONG, HDF5_CHUNKED, chunk_size, compression_level);
  free(topo_dset);  

  // Populate and write the offsets dataset
  long *offsets_dset;
  populate_offsets_dset(&offsets_dset, noffset, num_cells_loc+1, offset_offset, count, offset);
  write_dataset(group_id, count, offset, "Offsets", num_cells+npe(), num_cells_loc+1, 1, offsets_dset, H5T_NATIVE_LONG, HDF5_CHUNKED, chunk_size, compression_level);
  free(offsets_dset);  
  
  // Create subgroup "CellData"
  subgroup_id = H5Gcreate(group_id, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Allocate memory and write scalar datasets
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  for (scalar s in list) {
    populate_scalar_dset(s, scalar_dset, num_cells_loc, offset_cells, count, offset, cell_mask);
    write_dataset(subgroup_id, count, offset, s.name, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, HDF5_CHUNKED, chunk_size, compression_level);
  }
  free(scalar_dset);

  // Allocate memory and write vector datasets
  double *vector_dset = (double *)malloc(num_cells_loc * 3 * sizeof(double));
  for (vector v in vlist) {    
    populate_vector_dset(v, vector_dset, num_cells_loc, offset_cells, count, offset, cell_mask);
    write_dataset(subgroup_id, count, offset, v.x.name, num_cells, num_cells_loc, 3, vector_dset, H5T_NATIVE_DOUBLE, HDF5_CHUNKED, chunk_size, compression_level);
  }
  free(vector_dset);
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
  // HDF5 not available - print warning and return
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_vtkhdf_box() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif
