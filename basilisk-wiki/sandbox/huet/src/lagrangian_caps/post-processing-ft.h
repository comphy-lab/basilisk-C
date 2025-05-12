/**
# Post-processing functions for triangulated surfaces
*/

/** ## Output membrane position in plain text for external post-processing */
void dump_plain_nodes_pos(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nln; i++) {
      fprintf(file, ",%g %g %g", mesh->nodes[i].pos.x, mesh->nodes[i].pos.y,
      mesh->nodes[i].pos.z);
    }
    fprintf(file, "\n");
    fclose(file);
    }
}

void dump_plain_triangles(lagMesh* mesh, char* filename) {
  if (pid() == 0) {
    FILE* file = fopen(filename, "a");
    assert(file);
    fprintf(file, "%g", t);
    for(int i=0; i<mesh->nlt; i++) {
      fprintf(file, ",%d %d %d %g", mesh->triangles[i].node_ids[0],
      mesh->triangles[i].node_ids[1], mesh->triangles[i].node_ids[2],
      mesh->triangles[i].area);
    }
    fprintf(file, "\n");
    fclose(file);
  }
}

/** ## Visualization in paraview */
struct _pv_output_ascii {
  char* name;
  FILE* fp;
};

int pv_timestep = 0;
void pv_output_ascii(struct _pv_output_ascii p) {
  if (pid() == 0) {
    char name[128];
    char default_name[5] = "caps\0";
    char* prefix = p.name ? p.name : default_name;
    char suffix[64];
    sprintf(suffix, "_T%d.vtk", pv_timestep);
    sprintf(name, "%s%s", prefix, suffix);
    FILE* file = p.fp ? p.fp : fopen(name, "w");
    assert(file);

    /* Populate the header and other non-data fields */
    fprintf(file, "# vtk DataFile Version 4.2\n");
    fprintf(file, "Capsules at time %g\n", t);
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    
    /**
    First, we find the number of nodes and triangles to display
    */
    int nbpts_tot = 0;
    int nbtri_tot = 0;
    for(int j=0; j<NCAPS; j++) {
      if (CAPS(j).isactive) {
        nbpts_tot += CAPS(j).nln;
        nbtri_tot += CAPS(j).nlt;
      }
    }

    /* Populate the coordinates of all the Lagrangian nodes */
    fprintf(file, "POINTS %d double\n", nbpts_tot);
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nln; k++) {
        coord node_pos = correct_periodic_node_pos(
          CAPS(j).nodes[k].pos, CAPS(j).centroid);
        fprintf(file, "%g %g %g\n", node_pos.x, node_pos.y, node_pos.z);
      }
    }

    /* Populate the connectivity of the triangles */
    fprintf(file, "TRIANGLE_STRIPS %d %d\n", nbtri_tot, 
    4*nbtri_tot);
    int node_offset = 0;
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++) {
      fprintf(file, "%d %d %d %d\n", 3,
        CAPS(j).triangles[k].node_ids[0] + node_offset,
        CAPS(j).triangles[k].node_ids[1] + node_offset,
        CAPS(j).triangles[k].node_ids[2] + node_offset);
      }
      node_offset += CAPS(j).nln;
    }
    fprintf(file, "CELL_DATA %d\n", nbtri_tot);
    fprintf(file, "SCALARS T1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {  
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS T2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].tension[1]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tmax double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", max(CAPS(j).triangles[k].tension[0], 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Tavg double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", .5*(CAPS(j).triangles[k].tension[0] + 
          CAPS(j).triangles[k].tension[1]));
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_1 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[0]);
    }
    fprintf(file, "\n");
    fprintf(file, "SCALARS Lambda_2 double 1 \n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for(int j=0; j<NCAPS; j++) {
      for(int k=0; k<CAPS(j).nlt; k++)
        fprintf(file, "%g ", CAPS(j).triangles[k].stretch[1]);
    }
    fprintf(file, "\n");
    fclose(file);
  }
  pv_timestep++;
}