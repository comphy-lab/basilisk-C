#ifndef OUTPUT_HDF5_HELPERS_H
#define OUTPUT_HDF5_HELPERS_H

#ifdef HAVE_HDF5

#define HDF5_CONTIGUOUS 0
#define HDF5_CHUNKED 1
#define HDF5_ZFP 2

/** ### compute_chunk_size(): Centralized calculation for HDF5 chunking 

Targeting a balance between metadata overhead and I/O performance.
Ensures chunk dimensions are less than 2^32.

*/
hsize_t compute_chunk_size(hsize_t num_cells) {

  // Adaptively increase nchunks_per_proc so that chunks don't grow too large.
  // Start at 1 chunk per process; double until chunk_size fits within max_chunk_rows.
  const hsize_t max_chunk_rows = 1ULL << 20; // 1M rows (~8 MB for double, ~8 MB for long)

  hsize_t per_proc = (npe() > 0) ? num_cells / npe() : num_cells;
  int nchunks_per_proc = 1;
  while (nchunks_per_proc < (1 << 20) && per_proc / nchunks_per_proc > max_chunk_rows)
    nchunks_per_proc <<= 1;

  hsize_t chunk_size = per_proc / nchunks_per_proc;

  // Hard upper bound: HDF5 requires total elements per chunk < 2^32
  if (chunk_size >= 4294967296ULL) chunk_size = 4294967295ULL;

  if (pid() == 0) fprintf(stderr, "compute_chunk_size: %llu, nchunks %d\n", (unsigned long long)chunk_size, nchunks_per_proc);

  return chunk_size;
}

/** ### create_hdf5_file(): Helper to open/create the HDF5 file */
hid_t create_hdf5_file(char *name){
  hid_t file_id;
  hid_t acc_tpl1; // File access template

#if _MPI
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);

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
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
  H5Pclose(acc_tpl1);
#endif

  return file_id;
}

/** ### create_contiguous_dataset(): Creates a contiguous dataset in an HDF5 file */
void create_contiguous_dataset(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                               long num_cells, int num_cells_loc, int num_dims, const void *data,
                               hid_t datatype)
{
  // HDF5 cannot create a dataset with a zero-length dimension
  // (H5D__cache_dataspace_info: unable to get the next power of 2)
  if (num_cells == 0) return;

  hid_t dataspace_id, dataset_id, memspace_id, acc_tpl1;
  hsize_t dims2[2];
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

  // Create the dataset
  dataset_id = H5Dcreate2(file_id, dataset_name, datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Error creating dataset\n");
    H5Sclose(dataspace_id);
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
    return;
  }

  // Select hyperslab in the dataset
  dataspace_id = H5Dget_space(dataset_id);
  if (dataspace_id < 0) {
    fprintf(stderr, "Error getting dataspace\n");
    H5Dclose(dataset_id);
    H5Sclose(memspace_id);
    return;
  }
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  if (status < 0) {
    fprintf(stderr, "Error selecting hyperslab\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    return;
  }

#if _MPI
  // Create property list for collective dataset write
  acc_tpl1 = H5Pcreate(H5P_DATASET_XFER);
  if (acc_tpl1 < 0) {
    fprintf(stderr, "Error creating property list for collective dataset write\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    return;
  }
  status = H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);
  if (status < 0) {
    fprintf(stderr, "Error setting collective dataset write property\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(acc_tpl1);
    return;
  }

  // Write data to the dataset
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, acc_tpl1, data);
  if (status < 0) {
    fprintf(stderr, "Error writing data to dataset\n");
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    H5Pclose(acc_tpl1);
    return;
  }

  // Close all HDF5 objects to release resources
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(acc_tpl1);
#else
  // Write data to the dataset (serial)
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, H5P_DEFAULT, data);
  if (status < 0)
    fprintf(stderr, "Error writing data to dataset\n");

  // Close all HDF5 objects to release resources
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
#endif
}

/** ### create_chunked_dataset(): Creates a chunked dataset in an HDF5 file */
void create_chunked_dataset(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                            long num_cells, int num_cells_loc, int num_dims, const void *data, 
                            hid_t datatype, hsize_t chunk_size, int compression_level)
{
  // HDF5 cannot create a dataset with a zero-length dimension
  // (H5D__cache_dataspace_info: unable to get the next power of 2)
  if (num_cells == 0) return;

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

  // HDF5 constraint 1: chunk rows must not exceed the dataset rows for fixed-size datasets
  if (chunk_dims[0] > dims2[0] && dims2[0] > 0)
    chunk_dims[0] = dims2[0];

  // HDF5 constraint 2: total elements per chunk must be < 2^32
  if (chunk_dims[1] > 0 && chunk_dims[0] > 4294967295ULL / chunk_dims[1])
    chunk_dims[0] = 4294967295ULL / chunk_dims[1];

  status = H5Pset_chunk(plist_id, 2, chunk_dims);
  if (status < 0) {
    fprintf(stderr, "Error setting chunking properties\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Set the shuffle and compression properties
  status = H5Pset_shuffle(plist_id);
  if (status < 0) {
    fprintf(stderr, "Error setting shuffle properties\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

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

#if _MPI
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
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, acc_tpl1, data);
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
#else
  // Write data to the dataset (serial)
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, H5P_DEFAULT, data);
  if (status < 0)
    fprintf(stderr, "Error writing data to dataset\n");

  // Close all HDF5 objects to release resources
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(plist_id);
#endif
}

/** ### create_attribute_type(): Helper for string attributes */
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

/** ### create_attribute(): Helper for integer/scalar attributes */
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

/** ### write_simple_dataset(): Helper for small 1D datasets */
herr_t write_simple_dataset(hid_t group_id, const char *dataset_name, const long *data, const hsize_t *dims) {
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
  status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
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

#endif // HAVE_HDF5
#endif // OUTPUT_HDF5_HELPERS_H
