/** 
# Helper functions for output_vtkhdf.h 
*/ 

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

  MPI_Allreduce(num_points, num_points_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_cells,  num_cells_glob,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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

  MPI_Allreduce(num_points, num_points_glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(num_cells,  num_cells_glob,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

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

  // Perform an all-reduce operation to gather the number of points and cells from all subdomains
  MPI_Allreduce(list_points, offset_points, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(list_cells, offset_cells, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Calculate the offset for the points in the current subdomain
  offset[0] = 0;
  if (pid() != 0){
    // Sum the offsets of the previous subdomains to get the starting offset for the current subdomain
    for (int i = 1; i <= pid(); ++i){
      offset[0] += offset_points[i - 1];
    }
  }
}

void calculate_offsets2(int *offset_points, int num_points, hsize_t *offset) {
  // Arrays to store the number of points and cells in each subdomain
  int list_points[npe()];
  
  // Initialize the arrays to zero
  for (int i = 0; i < npe(); ++i){
    list_points[i] = 0;
  }

  // Set the number of points and cells for the current subdomain
  list_points[pid()] = num_points;
  
  // Perform an all-reduce operation to gather the number of points and cells from all subdomains
  MPI_Allreduce(list_points, offset_points, npe(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
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
void initialize_marker(vertex scalar marker, hsize_t *offset, hsize_t accumulate = 1 ) {
  int num_points = 0;
  foreach_vertex(serial, noauto){
    #if !TREE
      #if dimension == 2
      int _k = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
      #else
      int _k = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
      #endif
    #else // TREE
      int _k = num_points;
    #endif    
    marker[] = _k + offset[0]*accumulate;
    num_points++;
  }
  marker.dirty = true;
}

void initialize_marker_slice(vertex scalar marker, hsize_t *offset, coord n = {0, 0, 1}, double _alpha = 0, hsize_t accumulate = 1) {
  int num_points = 0;
  foreach_vertex(serial, noauto){
    marker[] = 0.;
    shortcut_slice(n, _alpha);
    marker[] = num_points + offset[0]*accumulate;
    num_points++;
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
  foreach_vertex(serial, noauto){
    #if !TREE
      #if dimension == 2
        int _k = (point.i - 2) * ((1 << point.level) + 1) + (point.j - 2);
      #else
        int _k = (point.i - 2) * sq((1 << point.level) + 1) + (point.j - 2) * ((1 << point.level) + 1) + (point.k - 2);
      #endif
    #endif

    // Calculate starting index
    int ii = _k * 3;

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

/** ### Populate types_dset  */ 
void populate_types_dset(char **types_dset, char type, int num_cells, int *offset_cells, hsize_t *count, hsize_t *offset) {
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

  // Allocate memory for types_dset
  *types_dset = (char *)malloc(count[0] * count[1] * sizeof(char));

  for (int i = 0; i < num_cells; ++i){
    (*types_dset)[i] = type;      
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

  foreach (serial, noauto){
    if (per_mask[]){
      #if !TREE
        #if dimension == 2
          int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
        #else
          int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
        #endif
      #endif

      // Store values
      scalar_dset[_k] = s[];
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

  num_cells = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      if (n.x == 1)
        scalar_dset[num_cells] = 0.5 * (val(s) + val(s, 1, 0, 0));
      else if (n.y == 1)
        scalar_dset[num_cells] = 0.5 * (val(s) + val(s, 0, 1, 0));
      else
        scalar_dset[num_cells] = 0.5 * (val(s) + val(s, 0, 0, 1));
      num_cells++;
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

  foreach (serial, noauto){
    if (per_mask[]){
      #if !TREE
        #if dimension == 2
          int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
        #else
          int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
        #endif
      #endif

      // Calculate starting index
      int ii = _k * 3;

      // Store each component
      vector_dset[ii + 0] = v.x[];
      vector_dset[ii + 1] = v.y[];
      #if dimension == 2
        vector_dset[ii + 2] = 0.;
      #else
        vector_dset[ii + 2] = v.z[];
      #endif
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

  num_cells = 0;
  foreach (serial, noauto){
    if (per_mask[]){
      int ii = num_cells * 3;
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
      num_cells++;
    }
  }
}
#endif

/** ### Populate offsets_dset  */ 
void populate_offsets_dset(long **offsets_dset, char noffset, int num_cells, int *offset_cells, hsize_t *count, hsize_t *offset) {
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

  // Allocate memory for topo_dset
  *offsets_dset = (long *)malloc(count[0] * count[1] * sizeof(long));

  for (int i = 0; i < num_cells; ++i){   
    (*offsets_dset)[i] = (long)i * (long)noffset;      
  }

}

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
  foreach (serial, noauto){
    if (per_mask[]){
      // _k exist by default on quad/octrees, but not on multigrid
      #if !TREE
        #if dimension == 2
          // Calculate index for 2D
          int _k = (point.i - 2) * ((1 << point.level)) + (point.j - 2);
        #else
          // Calculate index for 3D
          int _k = (point.i - 2) * sq((1 << point.level)) + (point.j - 2) * ((1 << point.level)) + (point.k - 2);
        #endif
      #endif

      // Calculate starting index for topo_dset
      int ii = _k * count[1];

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
  }
  count[0] = num_cells*pow(2, dimension);
  count[1] = 1;
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
  count[0] = num_cells*pow(2, dimension - 1);
  count[1] = 1;
}

/** ## Write Dataset 

*/

herr_t create_attribute_type(hid_t group_id, const char *attrname_type, const char *attrvalue_type, size_t str_size) {
  hid_t space_id, strtype, attr_id;
  herr_t status;

  // Create a scalar dataspace
  space_id = H5Screate(H5S_SCALAR);
  if (space_id < 0) {
    fprintf(stderr, "Failed to create scalar dataspace\n");
    return -1;
  }

  // Copy the string datatype and set its properties
  strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    fprintf(stderr, "Failed to copy string datatype\n");
    H5Sclose(space_id);
    return -1;
  }

  status = H5Tset_size(strtype, str_size);
  if (status < 0) {
    fprintf(stderr, "Failed to set string size\n");
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  status = H5Tset_strpad(strtype, H5T_STR_NULLTERM);
  if (status < 0) {
    fprintf(stderr, "Failed to set string padding\n");
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  status = H5Tset_cset(strtype, H5T_CSET_ASCII);
  if (status < 0) {
    fprintf(stderr, "Failed to set character set\n");
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  // Create the attribute
  attr_id = H5Acreate2(group_id, attrname_type, strtype, space_id, H5P_DEFAULT, H5P_DEFAULT);
  if (attr_id < 0) {
    fprintf(stderr, "Failed to create attribute\n");
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  // Write the attribute value
  status = H5Awrite(attr_id, strtype, attrvalue_type);
  if (status < 0) {
    fprintf(stderr, "Failed to write attribute value\n");
    H5Aclose(attr_id);
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  // Close the attribute
  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close attribute\n");
    H5Tclose(strtype);
    H5Sclose(space_id);
    return -1;
  }

  // Close the datatype
  status = H5Tclose(strtype);
  if (status < 0) {
    fprintf(stderr, "Failed to close string datatype\n");
    H5Sclose(space_id);
    return -1;
  }

  // Close the dataspace
  status = H5Sclose(space_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close scalar dataspace\n");
    return -1;
  }

  return 0;
}


// Function to create and write an attribute to an HDF5 group
herr_t create_attribute(hid_t group_id, const char *attrname_version, const int *version_data, const hsize_t *dims) {
  hid_t space_id, attr_id;
  herr_t status;

  // Create a simple dataspace
  space_id = H5Screate_simple(1, dims, NULL);
  if (space_id < 0) {
    fprintf(stderr, "Failed to create simple dataspace\n");
    return -1;
  }

  // Create the attribute
  attr_id = H5Acreate2(group_id, attrname_version, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT);
  if (attr_id < 0) {
    fprintf(stderr, "Failed to create attribute\n");
    H5Sclose(space_id);
    return -1;
  }

  // Write the attribute value
  status = H5Awrite(attr_id, H5T_NATIVE_INT, version_data);
  if (status < 0) {
    fprintf(stderr, "Failed to write attribute value\n");
    H5Aclose(attr_id);
    H5Sclose(space_id);
    return -1;
  }

  // Close the attribute
  status = H5Aclose(attr_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close attribute\n");
    H5Sclose(space_id);
    return -1;
  }

  // Close the dataspace
  status = H5Sclose(space_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close simple dataspace\n");
    return -1;
  }

  return 0;
}

// Function to create and write a simple dataset to an HDF5 group
herr_t write_simple_dataset(hid_t group_id, const char *dataset_name, const int *data, const hsize_t *dims) {
  hid_t space_id, dataset_id;
  herr_t status;

  // Create a simple dataspace
  space_id = H5Screate_simple(1, dims, NULL);
  if (space_id < 0) {
    fprintf(stderr, "Failed to create simple dataspace\n");
    return -1;
  }

  // Create the dataset
  dataset_id = H5Dcreate(group_id, dataset_name, H5T_NATIVE_LONG, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Failed to create dataset\n");
    H5Sclose(space_id);
    return -1;
  }

  // Write the dataset value
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    fprintf(stderr, "Failed to write dataset value\n");
    H5Dclose(dataset_id);
    H5Sclose(space_id);
    return -1;
  }

  // Close the dataset
  status = H5Dclose(dataset_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close dataset\n");
    H5Sclose(space_id);
    return -1;
  }

  // Close the dataspace
  status = H5Sclose(space_id);
  if (status < 0) {
    fprintf(stderr, "Failed to close simple dataspace\n");
    return -1;
  }

  return 0;
}

/** ### create_chunked_dataset(): Creates a chunked dataset in an HDF5 file
  
 The arguments and their default values are:

 *   **file_id**: HDF5 file identifier
 *   **count**: Size of dataset to create
 *   **offset**: Starting position for dataset creation
 *   **dataset_name**: Name of dataset to create
 *   **num_cells**: Total number of cells in dataset
 *   **num_cells_loc**: Number of cells to create in this call
 *   **num_dims**: Number of dimensions in dataset
 *   **topo_dset**: Pointer to data to write to dataset
 *   **datatype**: Data type of data to write to dataset
 *   **chunk_size**: Size of chunks in which dataset will be stored
 *   **compression_level**: Compression level (`default=6`)
 */
void create_chunked_dataset(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                            int num_cells, int num_cells_loc, int num_dims, const void *topo_dset, 
                            hid_t datatype, int chunk_size = num_cells_loc, int compression_level = 9)
{
  hid_t dataspace_id, dataset_id, memspace_id, plist_id, acc_tpl1;
  hsize_t dims2[2];
  hsize_t chunk_dims[2];
  herr_t status;

  // Define dimensions
  dims2[0] = num_cells;
  dims2[1] = num_dims;

  // Create the dataspace
  dataspace_id = H5Screate_simple(2, dims2, NULL);
  if (dataspace_id < 0) {
    fprintf(stderr, "Error creating dataspace\n");
    return;
  }

  // Create the dataset creation property list and set the chunking properties
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  if (plist_id < 0) {
    fprintf(stderr, "Error creating dataset creation property list\n");
    H5Sclose(dataspace_id);
    return;
  }
  chunk_dims[0] = chunk_size;
  chunk_dims[1] = dims2[1];
  status = H5Pset_chunk(plist_id, 2, chunk_dims);
  if (status < 0) {
    fprintf(stderr, "Error setting chunking properties\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Set the compression properties
  status = H5Pset_deflate(plist_id, compression_level);
  if (status < 0) {
    fprintf(stderr, "Error setting compression properties\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Create the dataset with chunking and compression properties
  dataset_id = H5Dcreate2(file_id, dataset_name, datatype, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Error creating dataset\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }
  H5Sclose(dataspace_id);

  // Define memory space for the dataset
  count[0] = num_cells_loc;
  count[1] = dims2[1];
  memspace_id = H5Screate_simple(2, count, NULL);
  if (memspace_id < 0) {
    fprintf(stderr, "Error creating memory space\n");
    H5Dclose(dataset_id);
    H5Pclose(plist_id);
    return;
  }

  // Select hyperslab in the dataset
  dataspace_id = H5Dget_space(dataset_id);
  if (dataspace_id < 0) {
    fprintf(stderr, "Error getting dataspace\n");
    H5Dclose(dataset_id);
    H5Sclose(memspace_id);
    H5Pclose(plist_id);
    return;
  }
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  if (status < 0) {
    fprintf(stderr, "Error selecting hyperslab\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Create property list for collective dataset write
  acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
  if (acc_tpl1 < 0) {
    fprintf(stderr, "Error creating property list for collective dataset write\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(plist_id);
    return;
  }
  status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
  if (status < 0) {
    fprintf(stderr, "Error setting collective dataset write property\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(plist_id);
    H5Pclose(acc_tpl1);
    return;
  }

  // Write data to the dataset
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, acc_tpl1, topo_dset);
  if (status < 0) {
    fprintf(stderr, "Error writing data to dataset\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(plist_id);
    H5Pclose(acc_tpl1);
    return;
  }

  // Close all HDF5 objects to release resources
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(plist_id);
  H5Pclose(acc_tpl1);
}






































