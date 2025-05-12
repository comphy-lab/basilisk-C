/**

 echo "@include <catalyst.h>" > $BASILISK/ast/std/catalyst.h

 and declade the datatypes at the end of $BASILISK/ast/defaults.h

 echo "typedef conduit_node;" >> $BASILISK/ast/defaults.h

 see https://groups.google.com/g/basilisk-fr/c/CM270hBSfWo
 */

#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h
#pragma autolink -I${CATALYST_INCDIR} -L${CATALYST_LIBDIR} -lcatalyst -DPARAVIEW_IMPL_DIR=\"${CATALYST_LIBDIR}/catalyst\"
@include <catalyst.h>

//-----------------------------------------------------------------------------
/**
 * Initialize Catalyst.
 */
//-----------------------------------------------------------------------------
void do_catalyst_initialization(int argc, char* argv[]){
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
      snprintf(buf, 256, "catalyst/scripts/script%d", (cc - 1));
      conduit_node_set_path_char8_str(catalyst_init_params, buf, argv[cc]);
    }
  }
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/implementation", "paraview");
  conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/search_paths/paraview", PARAVIEW_IMPL_DIR);
  enum catalyst_status err = catalyst_initialize(catalyst_init_params);
  conduit_node_destroy(catalyst_init_params);
  if (err != catalyst_status_ok)
    printf("Failed to initialize Catalyst: %d\n", err);
}

//-----------------------------------------------------------------------------
/**
 * Execute per cycle
 */
//-----------------------------------------------------------------------------
void do_catalyst_execute(int cycle, double time, CGrid* cgrid, CAttributes* attribs) {
  conduit_node* catalyst_exec_params = conduit_node_create();
  conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", cycle);
  conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", time);

  // the data must be provided on a named channel. the name is determined by the
  // simulation. for this one, we're calling it "grid".

  // declare the type of the channel; we're using Conduit Mesh Blueprint
  // to describe the mesh and fields.
  conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh");

  // now, create the mesh.
  conduit_node* mesh = conduit_node_create();

  // add coordsets
  conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "explicit");
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/x",
    /*data=*/cgrid->points, /*num_elements=*/cgrid->no_points, /*offset=*/0,
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/y",
    /*data=*/cgrid->points, /*num_elements=*/cgrid->no_points, /*offset=*/1 * sizeof(double),
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/z",
    /*data=*/cgrid->points, /*num_elements=*/cgrid->no_points, /*offset=*/2 * sizeof(double),
    /*stride=*/3 * sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);

  // add topologies
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "unstructured");
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords");
  conduit_node_set_path_char8_str(mesh, "topologies/mesh/elements/shape", "hex");
  conduit_node_set_path_external_int64_ptr(mesh, "topologies/mesh/elements/connectivity", cgrid->cells, cgrid->no_cells * pow(2, dimension));

  //add velocity (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/velocity/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/velocity/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/velocity/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "fields/velocity/values/x",
    /*data=*/attribs->Velocity, /*num_elements=*/cgrid->no_cells, /*offset=*/0,
    /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "fields/velocity/values/y",
    /*data=*/attribs->Velocity, /*num_elements=*/cgrid->no_cells,
    /*offset=*/cgrid->no_cells * sizeof(double),
    /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_float64_ptr_detailed(mesh, "fields/velocity/values/z",
    /*data=*/attribs->Velocity, /*num_elements=*/cgrid->no_cells,
    /*offset=*/2 * cgrid->no_cells * sizeof(double),
    /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);

  //add pressure (cell-field)
  conduit_node_set_path_char8_str(mesh, "fields/pressure/association", "element");
  conduit_node_set_path_char8_str(mesh, "fields/pressure/topology", "mesh");
  conduit_node_set_path_char8_str(mesh, "fields/pressure/volume_dependent", "false");
  conduit_node_set_path_external_float64_ptr_detailed(
    mesh, "fields/pressure/values",
    /*data=*/attribs->Pressure, /*num_elements=*/cgrid->no_cells,
    /*offset=*/0 * cgrid->no_cells * sizeof(double),
    /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    /*endianness*/ CONDUIT_ENDIANNESS_DEFAULT_ID);
  conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh);

  // print for debugging purposes, if needed
  // conduit_node_print(catalyst_exec_params);

  // print information with details about memory allocation
  // conduit_node* info = conduit_node_create();
  // conduit_node_info(catalyst_exec_params, info);
  // conduit_node_print(info);
  // conduit_node_destroy(info);

  enum catalyst_status err = catalyst_execute(catalyst_exec_params);
  if (err != catalyst_status_ok)
    printf("Failed to execute Catalyst: %d\n", err);

  conduit_node_destroy(catalyst_exec_params);
  conduit_node_destroy(mesh);
}

//-----------------------------------------------------------------------------
/**
 * Finalize Catalyst.
 */
//-----------------------------------------------------------------------------
void do_catalyst_finalization(){
  conduit_node* catalyst_fini_params = conduit_node_create();
  enum catalyst_status err = catalyst_finalize(catalyst_fini_params);
  if (err != catalyst_status_ok)
    printf("Failed to execute Catalyst: %d\n", err);
  conduit_node_destroy(catalyst_fini_params);
}

#endif
