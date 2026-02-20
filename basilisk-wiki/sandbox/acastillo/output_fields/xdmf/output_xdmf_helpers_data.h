#ifndef OUTPUT_XDMF_HELPERS_DATA_H
#define OUTPUT_XDMF_HELPERS_DATA_H

/** 
## Functions to write heavy data 
*/

#ifdef HAVE_HDF5

#include "acastillo/output_fields/output_hdf5_helpers.h"

/** ### write_dataset(): Wrapper for the writing routines

Provides a single entry point to write datasets with different compression modes.

*/
trace
void write_dataset(hid_t file_id, hsize_t *count, hsize_t *offset, const char *dataset_name,
                   long num_cells, long num_cells_loc, int num_dims, const void *data, 
                   hid_t datatype, int mode = HDF5_CHUNKED, hsize_t chunk_size = 0, int compression_level = 6)
{
  if (mode == HDF5_CONTIGUOUS || chunk_size == 0) {
    create_contiguous_dataset(file_id, count, offset, dataset_name, num_cells, num_cells_loc, num_dims, data, datatype);
  }
  else if (mode == HDF5_CHUNKED) {
    create_chunked_dataset(file_id, count, offset, dataset_name, num_cells, num_cells_loc, num_dims, data, datatype, chunk_size, compression_level);
  }
}

#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_DATA_H
