#ifndef OUTPUT_XDMF_HELPERS_DATA_H
#define OUTPUT_XDMF_HELPERS_DATA_H

/** 
## Functions to write heavy data 
*/

#ifdef HAVE_HDF5


/** ## Write Dataset
 
### create_contiguous_dataset(): Create a contiguous dataset in an HDF5 file
  
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
 */
void create_contiguous_dataset(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                               int num_cells, int num_cells_loc, int num_dims, const void *topo_dset,
                               hid_t datatype)
{
  hid_t dataspace_id, dataset_id, memspace_id, acc_tpl1;
  hsize_t dims2[2];
  herr_t status;

  // Define dimensions
  dims2[0] = num_cells;
  dims2[1] = num_dims;

  // Create the dataspace
  dataspace_id = H5Screate_simple(2, dims2, NULL);
  if(dataspace_id < 0) {
    fprintf(stderr, "Failed to create dataspace\n");
    return;
  }

  // Create the dataset with chunking properties
  dataset_id = H5Dcreate2(file_id, dataset_name, datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(dataspace_id);
  if(dataset_id < 0) {
    fprintf(stderr, "Failed to create dataset\n");
    return;
  }

  // Define memory space for the dataset
  count[0] = num_cells_loc;
  count[1] = dims2[1];
  memspace_id = H5Screate_simple(2, count, NULL);
  if(memspace_id < 0) {
    fprintf(stderr, "Failed to create memory space\n");
    H5Dclose(dataset_id);
    return;
  }

  // Select hyperslab in the dataset
  dataspace_id = H5Dget_space(dataset_id);
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  if(status < 0) {
    fprintf(stderr, "Failed to select hyperslab\n");
    H5Dclose(dataset_id);
    H5Sclose(memspace_id);
    return;
  }

  // Create property list for collective dataset write
  acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
  if(acc_tpl1 < 0) {
    fprintf(stderr, "Failed to create property list\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    return;
  }

  status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
  if(status < 0) {
    fprintf(stderr, "Failed to set MPIO property\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(acc_tpl1);
    return;
  }

  // Write data to the dataset
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, acc_tpl1, topo_dset);
  if(status < 0) {
    fprintf(stderr, "Failed to write data to dataset\n");
  }

  // Close all HDF5 objects to release resources
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(acc_tpl1);
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
  if(status < 0) {
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

/** 
### create_xdmf_file(): Helper to open/create the HDF5 file
*/
hid_t create_xdmf_file(char *name){
  hid_t file_id;
  hid_t acc_tpl1; // File access template

#if _MPI
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);
  H5Pset_file_locking(acc_tpl1, 0, 0);

  // Enable collective metadata operations for better parallel I/O performance (HDF5 1.10.0+)
  #if (H5_VERS_MAJOR > 1) || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 10)
    H5Pset_coll_metadata_write(acc_tpl1, 1);
    H5Pset_all_coll_metadata_ops(acc_tpl1, 1);
  #endif

  // Create a new HDF5 file collectively
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
  H5Pclose(acc_tpl1);
#else
  // Create a new HDF5 file without parallel I/O
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_file_locking(acc_tpl1, 0, 0);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
  H5Pclose(acc_tpl1);
#endif

  return file_id;
}

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_DATA_H
