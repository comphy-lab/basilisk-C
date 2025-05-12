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

#pragma autolink -lhdf5
#include <hdf5.h>

/** 
## preamble: define some useful macros
*/
#define shortcut_slice(n, _alpha)                                \
  double alpha = (_alpha - n.x * x - n.y * y - n.z * z) / Delta; \
  if (fabs(alpha) > 0.87)                                        \
    continue;

#include "output_xdmf_helpers.h"

/** 
## write_xdmf_light_data(): write an `.xdmf` file containing the light data
*/

void write_xdmf_light_data(scalar *slist, vector *vlist, char *file_name, char *subname, int num_cells = 0, int num_points = 0, int dim = dimension){

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  // Open a file for writing
  char name[111];
  sprintf(name, "%s.xmf", subname);
  FILE *fp = fopen(name, "w");

  // Write the header of the xdmf file.
  write_xdmf_header(fp, file_name);
  write_xdmf_topology(fp, dim, num_cells, num_points, t);
  write_xdmf_attributes(fp, num_cells, slist, vlist);
  write_xdmf_footer(fp);

  // Close the file
  fflush(fp);
  fclose(fp);

#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

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
 * **compression_level** : Level of compression to use when writing data to the HDF5 file.
 
*/
trace void output_xmf(scalar *slist, vector *vlist, char *subname, int compression_level = 9){
  hid_t acc_tpl1;    // File access template
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  char name[111];                  // Buffer for file name construction
  sprintf(name, "%s.h5", subname); // Construct the HDF5 file name

  // Define a scalar mask for periodic conditions
  scalar per_mask[];
  foreach () {
    per_mask[] = 1.;
  }

  // Obtain the number of points and cells and get a marker to reconstruct the topology
  int num_points = 0, num_cells = 0, num_points_loc = 0, num_cells_loc = 0;
  count_points_and_cells(&num_points, &num_cells, &num_points_loc, &num_cells_loc, per_mask);

  // Calculate offsets for parallel I/O
  int offset_points[npe()], offset_cells[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  // Initialize marker for topology reconstruction
  vertex scalar marker[];
  initialize_marker(marker, offset);

  // Write the light data
  if (pid() == 0) {
    write_xdmf_light_data(slist, vlist, name, subname, num_cells, num_points);
  }

  // Write the heavy data
  // Setup file access template with parallel I/O access
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);

  // Create a new HDF5 file collectively
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);

  // Release file-access template
  H5Pclose(acc_tpl1);

  // Define chunk size for parallel I/O
  hsize_t chunk_size = num_cells / npe() / 8;

  // Populate and write the topology dataset
  long *topo_dset;
  populate_topo_dset(&topo_dset, num_cells_loc, offset_cells, count, offset, per_mask, marker);
  create_chunked_dataset(file_id, count, offset, "/Topology", num_cells, num_cells_loc, pow(2, dimension), topo_dset, H5T_NATIVE_LONG, chunk_size, compression_level);
  free(topo_dset);

  // Create group for mesh geometry data
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset(&points_dset, num_points_loc, offset_points, count, offset);
  create_chunked_dataset(group_id, count, offset, "/Geometry/Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
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
    create_chunked_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
  }
  free(scalar_dset);

  // Allocate memory and write vector datasets
  double *vector_dset = (double *)malloc(num_cells_loc * 3 * sizeof(double));
  for (vector v in vlist) {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", v.x.name);
    populate_vector_dset(v, vector_dset, num_cells_loc, offset_cells, count, offset, per_mask);
    create_chunked_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 3, vector_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
  }
  free(vector_dset);
  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fclose(file_id);
}

/** 
## output_xmf_slice(): Exports a 2D slice of a 3D field along x, y, or z.

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
 * **n** : Vector defining the plane
 * **alpha** : Intersect defining the plane
 * **compression_level** : Level of compression to use when writing data to the HDF5 file. 

*/

#if dimension == 3
trace void output_xmf_slice(scalar *slist, vector *vlist, char *subname, coord n = {0, 0, 1}, double _alpha = 0, int compression_level = 9){
  hid_t acc_tpl1;    // File access template
  hid_t file_id;     // HDF5 file ID
  hid_t group_id;    // HDF5 group ID
  hsize_t count[2];  // Hyperslab selection parameters
  hsize_t offset[2]; // Offset for hyperslab

  char name[111];                  // Buffer for file name construction
  sprintf(name, "%s.h5", subname); // Construct the HDF5 file name

  // Define a scalar mask to deal with solids and periodic conditions
  scalar per_mask[];
  foreach (){
    per_mask[] = 0.;
    shortcut_slice(n, _alpha);
    if (alpha > 0.){
#if EMBED
      per_mask[] = cs[];
#else
      per_mask[] = 1.;
#endif
    }
  }

  // Obtain the number of points and cells and get a marker to reconstruct the topology
  int num_points = 0, num_cells = 0, num_points_loc = 0, num_cells_loc = 0;
  count_points_and_cells_slice(&num_points, &num_cells, &num_points_loc, &num_cells_loc, per_mask, n, _alpha);

  // Calculate offsets for parallel I/O
  int offset_points[npe()], offset_cells[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  // Initialize marker for topology reconstruction
  vertex scalar marker[];
  initialize_marker_slice(marker, offset, n, _alpha);

  // Write the light data to an XDMF file
  if (pid() == 0)
    write_xdmf_light_data(slist, vlist, name, subname, num_cells, num_points, dim = dimension - 1);

  // Setup file access template with parallel I/O access
  acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);

  // Create a new HDF5 file collectively
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);

  // Release file-access template
  H5Pclose(acc_tpl1);

  // Define chunk size for parallel I/O
  hsize_t chunk_size = num_cells / npe() / 8;

  // Populate and write the topology dataset
  long *topo_dset;
  populate_topo_dset_slice(&topo_dset, num_cells_loc, offset_cells, count, offset, per_mask, marker, n, _alpha);
  create_chunked_dataset(file_id, count, offset, "/Topology", num_cells, num_cells_loc, pow(2, dimension - 1), topo_dset, H5T_NATIVE_LONG, chunk_size, compression_level);
  free(topo_dset);

  // Create group for mesh geometry data
  group_id = H5Gcreate(file_id, "Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Populate and write the points dataset
  double *points_dset;
  populate_points_dset_slice(&points_dset, num_points_loc, offset_points, count, offset, n, _alpha);
  create_chunked_dataset(group_id, count, offset, "/Geometry/Points", num_points, num_points_loc, 3, points_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
  free(points_dset);
  H5Gclose(group_id);

  // Create group for cell data
  group_id = H5Gcreate(file_id, "Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Allocate memory and write scalar datasets
  double *scalar_dset = (double *)malloc(num_cells_loc * sizeof(double));
  for (scalar s in slist)
  {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", s.name);
    populate_scalar_dset_slice(s, scalar_dset, num_cells_loc, offset_cells, count, offset, per_mask, n, _alpha);
    create_chunked_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 1, scalar_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
  }
  free(scalar_dset);

  // Allocate memory and write vector datasets
  double *vector_dset = (double *)malloc(num_cells_loc * 3 * sizeof(double));
  for (vector v in vlist)
  {
    char substamp[1024];
    sprintf(substamp, "/Cells/%s", v.x.name);
    populate_vector_dset_slice(v, vector_dset, num_cells_loc, offset_cells, count, offset, per_mask, n, _alpha);
    create_chunked_dataset(group_id, count, offset, substamp, num_cells, num_cells_loc, 3, vector_dset, H5T_NATIVE_DOUBLE, chunk_size, compression_level);
  }
  free(vector_dset);
  H5Gclose(group_id);

  // Close HDF5 resources
  H5Fclose(file_id);
}
#endif

/** ## postamble: delete macros */
#undef shortcut_slice
