/**
# XDMF File Format

These routines are compatible with the `XDMF Model and Format` which can be read
using Paraview or Visit. Data is split in two categories: Light data and Heavy
data. 

* Light data is stored using eXtensible Markup Language (`XML`) to describe the
  data model and the data format. 

* Heavy data is composed of large datasets stored using the Hierarchical Data
  Format HDF5. As the name implies, data is organized following a hierarchical
  structure. HDF5 files can be read without any prior knowledge of the stored
  data. The type, rank, dimension and other properties of each array are stored
  inside the file in the form of meta-data. Additional features include support
  for a large number of objects, file compression, a parallel I/O implementation
  through the MPI-IO or MPI POSIX drivers. Using this format requires the HDF5
  library, which is usually installed in most computing centers or may be
  installed locally through a repository. Linking is automatic but requires the
  environment variables `HDF5_INCDIR` and `HDF5_LIBDIR`, which are usually set
  when you load the module hdf5. You might also add something like this to your
  `Makefile`

```
    HDFLIBS=-I$(HDF5_INCDIR) -L$(HDF5_LIBDIR) -lhdf5 -lhdf5_hl
    LIBS+=$(HDFLIBS)
```

**To use HDF5 we need to "declare" the library in the qcc system**. This should 
not be the case, but HDF5 is a messy library that just won't pass through `qcc`

```
    echo "@include <hdf5.h>" > $BASILISK/ast/std/hdf5.h
    echo "@include <hdf5_hl.h>" > $BASILISK/ast/std/hdf5_hl.h
```

**We also need to declare the datatypes at the end of `$BASILISK/ast/defaults.h`**

```
    echo "typedef hid_t, hsize_t, herr_t, H5L_info_t;" >> $BASILISK/ast/defaults.h
```

An explanation can be found [here](https://groups.google.com/g/basilisk-fr/c/CM270hBSfWo) 
*/

// Check if HDF5 is available by testing if the header can be included
#if __has_include(<hdf5.h>)
  #define HAVE_HDF5 1
  #pragma autolink -lhdf5
  #include <hdf5.h>
#else
  #warning "HDF5 library not found. XDMF output functions will be disabled."
#endif

#include "output_xdmf_helpers.h"

/** 
## output_xmf(): Exports full 2D (or 3D) fields.

This function performs the following steps:

1. Constructs the HDF5 file name based on the input subname.
2. Counts the number of points and cells in the data and calculates offsets for parallel I/O.
3. Initializes a marker for topology reconstruction.
4. Writes the light data to the XDMF file.
5. Sets up file access template with parallel I/O access and creates a new HDF5 file collectively.
6. Populates and writes the topology, points, and cell data datasets to the HDF5 file.
7. Closes the HDF5 resources.

The arguments and their default values are:

 * **slist** : Pointer to an array of scalar data.
 * **vlist** : Pointer to an array of vector data.
 * **subname** : String used to construct the HDF5 file name.
 * **mode** : Writing mode (HDF5_CONTIGUOUS or HDF5_CHUNKED).
 * **compression_level** : Compression level for GZIP or rate for ZFP.
  
 */
trace void output_xmf(scalar *slist, vector *vlist, char *subname, int mode = HDF5_CHUNKED, int compression_level = 6){
#ifdef HAVE_HDF5
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  // Construct the HDF5 file name
  char name[260];  // Buffer for file name construction
  snprintf(name, sizeof(name), "%s.h5", subname); 

  // Define a scalar mask for periodic conditions
  scalar per_mask[];
  foreach () {
    per_mask[] = 1.;
  }

  // Obtain the number of points and cells and get a marker to reconstruct the topology
  long num_points = 0, num_cells = 0;
  long num_points_loc = 0, num_cells_loc = 0;
  count_points_and_cells(&num_points, &num_cells, &num_points_loc, &num_cells_loc, per_mask);

  // Calculate offsets for parallel I/O
  long offset_points[npe()], offset_cells[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  // Initialize marker for topology reconstruction
  vertex scalar marker[];
  initialize_marker(marker, offset);

  // Write the light data
  if (pid() == 0) {
    write_xdmf_light_data(slist, vlist, name, subname, num_cells, num_points);
  }

  // Write the heavy data
  // Create HDF5 file using helper
  file_id = create_hdf5_file(name);
  if (file_id < 0) return; // Exit if file creation failed

  // Centralized chunk size calculation
  hsize_t chunk_size = compute_chunk_size(num_cells);

  // Populate and write the topology dataset
  long *topo_dset;
  populate_topo_dset_xdmf(&topo_dset, num_cells_loc, offset_cells, count, offset, per_mask, marker);
  write_dataset(file_id, count, offset, "/Topology", num_cells, num_cells_loc, pow(2, dimension), topo_dset, H5T_NATIVE_LONG, mode, chunk_size, compression_level);
  free(topo_dset);

  // Create group for mesh geometry data
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset(&points_dset, num_points_loc, offset_points, count, offset);
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
    populate_scalar_dset(s, scalar_dset, num_cells_loc, offset_cells, count, offset, per_mask);
    write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, mode, chunk_size, compression_level);
  }
  free(scalar_dset);

  // Allocate memory and write vector datasets
  double *vector_dset = (double *)malloc(num_cells_loc * 3 * sizeof(double));
  for (vector v in vlist) {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", v.x.name);
    populate_vector_dset(v, vector_dset, num_cells_loc, offset_cells, count, offset, per_mask);
    write_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 3, vector_dset, H5T_NATIVE_DOUBLE, mode, chunk_size, compression_level);
  }
  free(vector_dset);
  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
#else
  // HDF5 not available - print warning and return
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_xmf() called but HDF5 is not available. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#if dimension > 2
#include "output_xdmf_slice.h"
#endif
#include "output_xdmf_box.h"
#include "output_xdmf_facets.h"
