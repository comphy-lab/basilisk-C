#ifndef OUTPUT_CATALYST_H
#define OUTPUT_CATALYST_H

/**
# Modernized Paraview Catalyst Integration

This wrapper completely replaces `deprecated/CatalystAdaptor.h` and `deprecated/FEDataStructures.h`.
It seamlessly translates the `output_fields` array generation logic directly into ParaView's
Conduit nodes with zero-copy overhead.
*/

#pragma autolink -I${CATALYST_INCDIR} -L${CATALYST_LIBDIR} -lcatalyst -DPARAVIEW_IMPL_DIR=${CATALYST_LIBDIR}/catalyst
@include <catalyst.h>



#if __has_include(<hdf5.h>)
  #define HAVE_HDF5 1
  #pragma autolink -lhdf5
  #include <hdf5.h>
#else
  #warning "HDF5 library not found. VTKHDF output functions will be disabled."
#endif

#include "acastillo/output_fields/output_common_helpers.h"
#include "acastillo/output_fields/vtkhdf/output_vtkhdf_helpers_populate.h"

//-----------------------------------------------------------------------------
/**
 * Initialize Catalyst using command line argument scripts.
 */
//-----------------------------------------------------------------------------
void output_catalyst_init(int argc, char* argv[]){
  conduit_node* catalyst_init_params = conduit_node_create();
  // pass scripts pass on command line.
  for (int cc = 1; cc < argc; ++cc){
    if (strcmp(argv[cc], "--output") == 0 && (cc + 1) < argc) {
      conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/type", "io");
      conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/filename", argv[cc + 1]);
      conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/channel", "grid");
      ++cc;
    }
    else {
      char buf[256];
      snprintf(buf, 256, "catalyst/scripts/script%d/filename", (cc - 1));
      conduit_node_set_path_char8_str(catalyst_init_params, buf, argv[cc]);
    }
  }
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/implementation", "paraview");
  
  const char* pv_libdir = getenv("CATALYST_LIBDIR");
  char impl_dir[1024];
  if (pv_libdir) {
      snprintf(impl_dir, sizeof(impl_dir), "%s/catalyst", pv_libdir);
  } else {
      snprintf(impl_dir, sizeof(impl_dir), "/usr/lib/catalyst");
  }
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/search_paths/paraview", impl_dir);

#if _MPI
  // Critical for cluster rendering! Hands off the Basilisk MPI communicator to ParaView
  conduit_node_set_path_int64(catalyst_init_params, "catalyst/mpi_comm", MPI_Comm_c2f(MPI_COMM_WORLD));
#endif

  enum catalyst_status err = catalyst_initialize(catalyst_init_params);
  conduit_node_destroy(catalyst_init_params);
  if (err != catalyst_status_ok)
    printf("Failed to initialize Catalyst: %d\n", err);
}

//-----------------------------------------------------------------------------
/**
 * Finalize Catalyst.
 */
//-----------------------------------------------------------------------------
void output_catalyst_finalize(){
  conduit_node* catalyst_fini_params = conduit_node_create();
  enum catalyst_status err = catalyst_finalize(catalyst_fini_params);
  if (err != catalyst_status_ok)
    printf("Failed to execute Catalyst finalizing: %d\n", err);
  conduit_node_destroy(catalyst_fini_params);
}

//-----------------------------------------------------------------------------
/**
 * Execute Catalyst processing for the current cycle.
 */
//-----------------------------------------------------------------------------
trace
void output_catalyst(scalar *slist, vector *vlist, int cycle, double time) {
#ifdef HAVE_HDF5
  conduit_node* catalyst_exec_params = conduit_node_create();
  if (pid() == 0) {
    fprintf(stderr, "catalyst_execute: Simulating cycle %d at time %g\n", cycle, time);
    fflush(stderr);
  }

  conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", cycle);
  conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", time);

  // declare the type of the channel
  conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh");
  conduit_node* mesh = conduit_node_create();

  // Define mask for valid domains (including periodicity treatment)
  scalar per_mask[];
  foreach () {
    per_mask[] = 1.;
    #if dimension == 2
      if (Period.y && y + Delta > Y0 + L0) per_mask[] = 0.;
      if (Period.x && x + Delta > X0 + L0) per_mask[] = 0.;
    #elif dimension == 3
      if (Period.z && z + Delta > Z0 + L0) per_mask[] = 0.;
      if (Period.y && y + Delta > Y0 + L0) per_mask[] = 0.;
      if (Period.x && x + Delta > X0 + L0) per_mask[] = 0.;
    #endif
  }

  long num_points = 0, num_points_loc = 0;
  long num_cells = 0, num_cells_loc = 0;
  count_points_and_cells(&num_points, &num_cells, &num_points_loc, &num_cells_loc, per_mask);

  hsize_t count[2], offset[2];
  long offset_points[npe()], offset_cells[npe()];
  calculate_offsets(offset_points, offset_cells, num_points_loc, num_cells_loc, offset);

  vertex scalar marker[];
  initialize_marker(marker, offset, 0);

  // Populate points using zero-copy allocation wrapper mapping
  double *points_dset;
  populate_points_dset(&points_dset, num_points_loc, offset_points, count, offset);

  // Populate topology (Connectivity)
  long *topo_dset;
  populate_topo_dset_vtkhdf(&topo_dset, num_cells_loc, offset_cells, count, offset, per_mask, marker);

  // Add Coordinate sets to Conduit
  conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "explicit");
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/x",
    /*data=*/points_dset, /*num_elements=*/num_points_loc, /*offset=*/0,
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/y",
    /*data=*/points_dset, /*num_elements=*/num_points_loc, /*offset=*/sizeof(double),
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/z",
    /*data=*/points_dset, /*num_elements=*/num_points_loc, /*offset=*/2 * sizeof(double),
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);

  // Add Topologies to Conduit
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "unstructured");
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords");
  #if dimension == 2
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/elements/shape", "quad");
  #else
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/elements/shape", "hex");
  #endif
  conduit_node_set_path_external_int64_ptr(mesh, "topologies/mesh/elements/connectivity", topo_dset, num_cells_loc * pow(2, dimension));

  // Build Field buffers
  int max_fields = 100;
  double *scalar_dsets[max_fields];
  double *vector_dsets[max_fields];
  
  // Scalar mapping
  int s_index = 0;
  for (scalar s in slist) {
    if (s_index >= max_fields) break;
    scalar_dsets[s_index] = (double *)malloc(num_cells_loc * sizeof(double));
    populate_scalar_dset(s, scalar_dsets[s_index], num_cells_loc, offset_cells, count, offset, per_mask);
    
    char path[1024];
    sprintf(path, "fields/%s/association", s.name);
    conduit_node_set_path_char8_str(mesh, path, "element");
    sprintf(path, "fields/%s/topology", s.name);
    conduit_node_set_path_char8_str(mesh, path, "mesh");
    sprintf(path, "fields/%s/volume_dependent", s.name);
    conduit_node_set_path_char8_str(mesh, path, "false");
    sprintf(path, "fields/%s/values", s.name);
    
    conduit_node_set_path_external_float64_ptr_detailed(
      mesh, path,
      /*data=*/scalar_dsets[s_index], /*num_elements=*/num_cells_loc,
      /*offset=*/0, /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
      /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
    s_index++;
  }

  // Vector mapping
  int v_index = 0;
  for (vector v in vlist) {
    if (v_index >= max_fields) break;
    vector_dsets[v_index] = (double *)malloc(num_cells_loc * 3 * sizeof(double));
    populate_vector_dset(v, vector_dsets[v_index], num_cells_loc, offset_cells, count, offset, per_mask);
    
    char path[1024];
    sprintf(path, "fields/%s/association", v.x.name);
    conduit_node_set_path_char8_str(mesh, path, "element");
    sprintf(path, "fields/%s/topology", v.x.name);
    conduit_node_set_path_char8_str(mesh, path, "mesh");
    sprintf(path, "fields/%s/volume_dependent", v.x.name);
    conduit_node_set_path_char8_str(mesh, path, "false");
    
    // Interleaved component mapping
    sprintf(path, "fields/%s/values/x", v.x.name);
    conduit_node_set_path_external_float64_ptr_detailed(mesh, path,
      /*data=*/vector_dsets[v_index], /*num_elements=*/num_cells_loc, /*offset=*/0,
      /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
      /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
      
    sprintf(path, "fields/%s/values/y", v.x.name);
    conduit_node_set_path_external_float64_ptr_detailed(mesh, path,
      /*data=*/vector_dsets[v_index], /*num_elements=*/num_cells_loc, /*offset=*/sizeof(double),
      /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
      /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
      
    sprintf(path, "fields/%s/values/z", v.x.name);
    conduit_node_set_path_external_float64_ptr_detailed(mesh, path,
      /*data=*/vector_dsets[v_index], /*num_elements=*/num_cells_loc, /*offset=*/2 * sizeof(double),
      /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
      /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
      
    v_index++;
  }

  conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh);

  // Execute
  enum catalyst_status err = catalyst_execute(catalyst_exec_params);
  if (err != catalyst_status_ok && pid() == 0) {
    fprintf(stderr, "catalyst_execute: WARNING! Catalyst Engine Failed Execution Check (Error Code %d)\n", err);
    fflush(stderr);
  }

  // Teardown API handles
  conduit_node_destroy(catalyst_exec_params);
  conduit_node_destroy(mesh);
  
  // Teardown allocated dynamic buffers
  free(points_dset);
  free(topo_dset);
  for (int i=0; i<s_index; i++) free(scalar_dsets[i]);
  for (int i=0; i<v_index; i++) free(vector_dsets[i]);
  
#else
  static int warning_printed = 0;
  if (!warning_printed && pid() == 0) {
    fprintf(stderr, "Warning: output_catalyst() called but HDF5 definitions are missing. Output skipped.\n");
    warning_printed = 1;
  }
#endif
}

#endif // OUTPUT_CATALYST_H
