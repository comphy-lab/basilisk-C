#ifndef OUTPUT_XDMF_HELPERS_ZFP_H
#define OUTPUT_XDMF_HELPERS_ZFP_H

#ifdef HAVE_HDF5

#ifndef H5Z_FILTER_ZFP
#define H5Z_FILTER_ZFP 32013
#endif

/** ### create_chunked_dataset_zfp(): Creates a chunked dataset using ZFP compression 
  *
  * (Moved to cold storage)
  */
void create_chunked_dataset_zfp(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                                int num_cells, int num_cells_loc, int num_dims, const void *data, 
                                hid_t datatype, int chunk_size, double precision)
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

  // Create the dataset creation property list
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  if (plist_id < 0) {
    fprintf(stderr, "Error creating dataset creation property list\n");
    H5Sclose(dataspace_id);
    return;
  }

  // Set the chunking properties (Required for ZFP)
  chunk_dims[0] = chunk_size;
  chunk_dims[1] = dims2[1];
  status = H5Pset_chunk(plist_id, 2, chunk_dims);
  if (status < 0) {
    fprintf(stderr, "Error setting chunking properties\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Set the ZFP filter
  unsigned int cd_values[10];
  for(int i=0; i<10; i++) cd_values[i] = 0;
  
  if (precision == 0) {
    cd_values[4] = 5; // H5Z_ZFP_MODE_REVERSIBLE (Lossless)
    status = H5Pset_filter(plist_id, H5Z_FILTER_ZFP, H5Z_FLAG_OPTIONAL, 5, cd_values);
  }
  else {
    cd_values[4] = 2; // ZFP_MODE_PRECISION
    cd_values[5] = (unsigned int)precision;
    status = H5Pset_filter(plist_id, H5Z_FILTER_ZFP, H5Z_FLAG_OPTIONAL, 6, cd_values);
  }

  if (status < 0) {
    fprintf(stderr, "Error setting ZFP filter. Ensure H5Z-ZFP plugin is installed.\n");
    H5Sclose(dataspace_id);
    H5Pclose(plist_id);
    return;
  }

  // Create the dataset
  dataset_id = H5Dcreate2(file_id, dataset_name, datatype, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Error creating dataset with ZFP\n");
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

  // Select hyperslab
  dataspace_id = H5Dget_space(dataset_id);
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
  H5Pset_dxpl_mpio(acc_tpl1, H5FD_MPIO_COLLECTIVE);

  // Write data to the dataset
  status = H5Dwrite(dataset_id, datatype, memspace_id, dataspace_id, acc_tpl1, data);
  if(status < 0) fprintf(stderr, "Error writing data to dataset\n");

  // Close all HDF5 objects
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  H5Pclose(plist_id);
  H5Pclose(acc_tpl1);
}

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_ZFP_H
