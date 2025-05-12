/**
# Dump and restore functions for triangulated membranes
The functions below write/read the lagrangian mesh to/from a file in order
to restart simulations.

At the moment these functions are only valid for three-dimensional simulations.
*/

void dump_lagnode(FILE* fp, lagNode* node) {
  foreach_dimension() fwrite(&(node->pos.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->lagVel.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->normal.x), sizeof(double), 1, fp);
  fwrite(&(node->curv), sizeof(double), 1, fp);
  fwrite(&(node->gcurv), sizeof(double), 1, fp);
  fwrite(&(node->ref_curv), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->lagForce.x), sizeof(double), 1, fp);
  fwrite(&(node->nb_neighbors), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fwrite(&(node->neighbor_ids[j]), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fwrite(&(node->edge_ids[j]), sizeof(int), 1, fp);
  fwrite(&(node->nb_triangles), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fwrite(&(node->triangle_ids[j]), sizeof(int), 1, fp);
  fwrite(&(node->nb_fit_iterations), sizeof(int), 1, fp);
}

void restore_lagnode(FILE* fp, lagNode* node) {
  foreach_dimension() fread(&(node->pos.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->lagVel.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->normal.x), sizeof(double), 1, fp);
  fread(&(node->curv), sizeof(double), 1, fp);
  fread(&(node->gcurv), sizeof(double), 1, fp);
  fread(&(node->ref_curv), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->lagForce.x), sizeof(double), 1, fp);
  fread(&(node->nb_neighbors), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fread(&(node->neighbor_ids[j]), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fread(&(node->edge_ids[j]), sizeof(int), 1, fp);
  fread(&(node->nb_triangles), sizeof(int), 1, fp);
  for(int j=0; j<6; j++)
    fread(&(node->triangle_ids[j]), sizeof(int), 1, fp);
  fread(&(node->nb_fit_iterations), sizeof(int), 1, fp);
}

void dump_edge(FILE* fp, Edge* edge) {
  fwrite(&(edge->node_ids[0]), sizeof(int), 1, fp);
  fwrite(&(edge->node_ids[1]), sizeof(int), 1, fp);
  fwrite(&(edge->triangle_ids[0]), sizeof(int), 1, fp);
  fwrite(&(edge->triangle_ids[1]), sizeof(int), 1, fp);
  fwrite(&(edge->l0), sizeof(double), 1, fp);
  fwrite(&(edge->length), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(edge->normal.x), sizeof(double), 1, fp);
}

void restore_edge(FILE* fp, Edge* edge) {
  fread(&(edge->node_ids[0]), sizeof(int), 1, fp);
  fread(&(edge->node_ids[1]), sizeof(int), 1, fp);
  fread(&(edge->triangle_ids[0]), sizeof(int), 1, fp);
  fread(&(edge->triangle_ids[1]), sizeof(int), 1, fp);
  fread(&(edge->l0), sizeof(double), 1, fp);
  fread(&(edge->length), sizeof(double), 1, fp);
  foreach_dimension() fread(&(edge->normal.x), sizeof(double), 1, fp);
}

void dump_triangle(FILE* fp, Triangle* triangle) {
  for(int k=0; k<3; k++) {
    fwrite(&(triangle->node_ids[k]), sizeof(int), 1, fp);
    fwrite(&(triangle->edge_ids[k]), sizeof(int), 1, fp);
  }
  fwrite(&(triangle->area), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(triangle->normal.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(triangle->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() {
    fwrite(&(triangle->refShape[0].x), sizeof(double), 1, fp);
    fwrite(&(triangle->refShape[1].x), sizeof(double), 1, fp);
  }
  for(int k=0; k<3; k++)
    for(int l=0; l<2; l++)
      fwrite(&(triangle->sfc[k][l]), sizeof(double), 1, fp);
  fwrite(&(triangle->stretch[0]), sizeof(double), 2, fp);
  fwrite(&(triangle->tension[0]), sizeof(double), 2, fp);
}

void restore_triangle(FILE* fp, Triangle* triangle) {
  for(int k=0; k<3; k++) {
    fread(&(triangle->node_ids[k]), sizeof(int), 1, fp);
    fread(&(triangle->edge_ids[k]), sizeof(int), 1, fp);
  }
  fread(&(triangle->area), sizeof(double), 1, fp);
  foreach_dimension() fread(&(triangle->normal.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(triangle->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() {
    fread(&(triangle->refShape[0].x), sizeof(double), 1, fp);
    fread(&(triangle->refShape[1].x), sizeof(double), 1, fp);
  }
  for(int k=0; k<3; k++)
    for(int l=0; l<2; l++)
      fread(&(triangle->sfc[k][l]), sizeof(double), 1, fp);
  fread(&(triangle->stretch[0]), sizeof(double), 2, fp);
  fread(&(triangle->tension[0]), sizeof(double), 2, fp);
}

void dump_lagmesh(FILE* fp, lagMesh* mesh) {
  fwrite(&(mesh->nln), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nln; i++) dump_lagnode(fp, &(mesh->nodes[i]));
  fwrite(&(mesh->nle), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nle; i++) dump_edge(fp, &(mesh->edges[i]));
  fwrite(&(mesh->nlt), sizeof(int), 1, fp);
  for(int i=0; i<mesh->nlt; i++) 
    dump_triangle(fp, &(mesh->triangles[i]));
  foreach_dimension() fwrite(&(mesh->centroid.x), sizeof(double), 1, fp);
  fwrite(&(mesh->volume), sizeof(double), 1, fp);
  fwrite(&(mesh->initial_volume), sizeof(double), 1, fp);
  int tmp;
  tmp = mesh->updated_stretches ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = mesh->updated_normals ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = mesh->updated_curvatures ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = mesh->isactive ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
}

void restore_lagmesh(FILE* fp, lagMesh* mesh) {
  fread(&(mesh->nln), sizeof(int), 1, fp);
  mesh->nodes = malloc(mesh->nln*sizeof(lagNode));
  for(int i=0; i<mesh->nln; i++) restore_lagnode(fp, &mesh->nodes[i]);
  fread(&(mesh->nle), sizeof(int), 1, fp);
  mesh->edges = malloc(mesh->nle*sizeof(Edge));
  for(int i=0; i<mesh->nle; i++) restore_edge(fp, &mesh->edges[i]);
  fread(&(mesh->nlt), sizeof(int), 1, fp);
  mesh->triangles = malloc(mesh->nlt*sizeof(Triangle));
  for(int i=0; i<mesh->nlt; i++) restore_triangle(fp, &mesh->triangles[i]);
  foreach_dimension() fread(&(mesh->centroid.x), sizeof(double), 1, fp);
  fread(&(mesh->volume), sizeof(double), 1, fp);
  fread(&(mesh->initial_volume), sizeof(double), 1, fp);
  int tmp;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_stretches = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_normals = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_curvatures = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->isactive = (tmp == 0) ? false : true;
}

/** If the simulation contains several membranes, we dump and read one mesh at
a time, but all meshes are stored in the same file. */
struct _dump_capsules {
    char* name;
    FILE* fp;
};

void dump_capsules(struct _dump_capsules p) {
    char default_name[10] = "caps.dump\0";
    char* name = p.name ? p.name : default_name;
    FILE* file = p.fp ? p.fp : fopen(name, "w");
    assert(file);
    for(int i=0; i<NCAPS; i++) dump_lagmesh(file, &CAPS(i));
    fclose(file);
}

void restore_capsules(char* filename) {
  FILE* file = fopen(filename, "r");
  assert(file);
  for(int i=0; i<NCAPS; i++) restore_lagmesh(file, &CAPS(i));
  fclose(file);
  initialize_all_capsules_stencils();
  generate_lag_stencils(no_warning = true);
}
