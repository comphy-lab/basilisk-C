/**
# Toolbox to perform operations on a triangulated meshes

From defining geometric computations such as normal vectors, volume and 
centroid, useful macros, subdividing triangles, below is a collection of helpful 
functions to deal with triangular meshes.
*/

/**
## Geometric computations

The function below computes the length of an edge. It takes as arguments
a pointer to the mesh as well as the ID of the edge of interest.
*/
double edge_length(lagMesh* mesh, int i) {
  double length = 0.;
  int v1, v2;
  v1 = mesh->edges[i].node_ids[0];
  v2 = mesh->edges[i].node_ids[1];
  foreach_dimension() {
    length += sq(GENERAL_1DIST(mesh->nodes[v1].pos.x, mesh->nodes[v2].pos.x));
  }
  return sqrt(length);
}

/**
The function ```compute_lengths``` below computes the lengths of all edges. It
takes as an argument a pointer to the mesh. If the optional argument
```force``` is set to ```true```, the edges' lengths are computed no matter the value of ```updated_stretches```.
*/
struct _compute_lengths{
  lagMesh* mesh;
  bool force;
};

void compute_lengths(struct _compute_lengths p) {
  lagMesh* mesh = p.mesh;
  bool force = (p.force) ? p.force : false;
  if (force || !mesh->updated_stretches) {
    for(int i=0; i < mesh->nle; i++)
      mesh->edges[i].length = edge_length(mesh, i);
    mesh->updated_stretches = true;
  }
}

#if dimension < 3
/**
The two functions below compute the outward normal vector to all the edges of
a Lagrangian mesh, for 2D simulations.
*/
void comp_edge_normal(lagMesh* mesh, int i) {
  int node_id[2];
  for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].node_ids[j];
  mesh->edges[i].normal.y = GENERAL_1DIST(mesh->nodes[node_id[0]].pos.x,
    mesh->nodes[node_id[1]].pos.x);
  mesh->edges[i].normal.x = GENERAL_1DIST(mesh->nodes[node_id[1]].pos.y,
    mesh->nodes[node_id[0]].pos.y);
  double normn = sqrt(sq(mesh->edges[i].normal.x)
    + sq(mesh->edges[i].normal.y));
  foreach_dimension() mesh->edges[i].normal.x /= normn;
}

void comp_edge_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nle; i++) comp_edge_normal(mesh, i);
}
#else // dimension > 2
/** In 3D simulations, the function below assumes that the Lagrangian mesh contains the origin and
is convex, and swaps the order of the nodes in order to compute an outward
normal vector. This only need to be performed at the creation of the mesh since
the outward property of the normal vectors won't change through the simulation.
*/
void comp_initial_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) {
    int nid[3]; // node ids
    coord centroid; /** Note: the centroid is only valid if the triangle is not
    across periodic boundaries, which is fine for this function since it is
    assumed the center of the membrane is at the origin. */
    foreach_dimension() centroid.x = 0.;
    for(int j=0; j<3; j++) {
      nid[j] = mesh->triangles[i].node_ids[j];
      foreach_dimension() centroid.x += mesh->nodes[nid[j]].pos.x/3;
    }
    coord normal, e[2];
    for(int j=0; j<2; j++)
      foreach_dimension()
        e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
          mesh->nodes[nid[j+1]].pos.x);
    foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
    double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
    double dp = 0.; // dp for "dot product"
    foreach_dimension() dp += normal.x*centroid.x;
    /** If the dot product is negative, the computed normal is inward and we
    need to swap two nodes of the triangle.*/
    if (dp < 0) {
      mesh->triangles[i].node_ids[1] = nid[2];
      mesh->triangles[i].node_ids[2] = nid[1];
      foreach_dimension() normal.x *= -1;
    }
    foreach_dimension() {
      mesh->triangles[i].centroid.x = centroid.x;
      mesh->triangles[i].normal.x = normal.x/norm;
    }
    mesh->triangles[i].area = norm/2;
  }
}

/** The two functions below compute the outward normal vector to all the
triangles of the mesh, for 3D simulations. */
void comp_triangle_area_normal(lagMesh* mesh, int i) {
  int nid[3]; // node ids
  for(int j=0; j<3; j++) nid[j] = mesh->triangles[i].node_ids[j];
  /** The next 15 lines compute the centroid of the triangle, making sure it is
  valid when the triangle lies across periodic boundaries. */
  foreach_dimension() mesh->triangles[i].centroid.x = 0.;
  for(int j=0; j<3; j++) {
    foreach_dimension() {
      mesh->triangles[i].centroid.x +=
        ACROSS_PERIODIC(mesh->nodes[nid[j]].pos.x,
        mesh->nodes[nid[0]].pos.x) ? mesh->nodes[nid[j]].pos.x/3 - L0 :
        mesh->nodes[nid[j]].pos.x/3;
    }
  }
  coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
  foreach_dimension() {
    if (fabs(mesh->triangles[i].centroid.x - origin.x) > L0/2.) {
      if (mesh->triangles[i].centroid.x - origin.x > 0)
        mesh->triangles[i].centroid.x -= L0;
      else mesh->triangles[i].centroid.x += L0;
    }
  }
  coord normal, e[2];
  for(int j=0; j<2; j++)
    foreach_dimension()
      e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
        mesh->nodes[nid[j+1]].pos.x);
  foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
  double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
  foreach_dimension() mesh->triangles[i].normal.x = normal.x/norm;
  mesh->triangles[i].area = norm/2.;
}

void comp_triangle_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) comp_triangle_area_normal(mesh, i);
}
#endif

/**
If a Lagrangian node falls exactly on an edge or a vertex of the Eulerian
mesh, some issues arise when checking for periodic boundary conditions. As a
quick fix, if this is the case we shift the point position by $10^{-10}$, as
is done in the two functions below.
*/
bool on_face(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) return true;
  else return false;
}

void correct_node_pos(coord* node) {
  coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
  foreach_dimension() {
    if (on_face(node->x, N, L0))
      node->x += 1.e-10;
    //FIXME: the nodes should not be sent to the other side of the domain if the boundary is not periodic...
    if (node->x > origin.x + L0/2)
      node->x -= L0;
    else if (node->x < origin.x - L0/2)
      node->x += L0;
  }
}

void correct_lag_pos(lagMesh* mesh) {
  for(int i=0; i < mesh->nln; i++) {
    correct_node_pos(&mesh->nodes[i].pos);
  }
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}

/**
The function below computes the centroid of the capsule as the average of the
coordinates of all its nodes. The centroid is stored as an attribute of the
caps structure.
*/
void comp_centroid(lagMesh* mesh) {
  coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
  foreach_dimension() mesh->centroid.x = 0.;
  for(int i=0; i<mesh->nln; i++)
    foreach_dimension() {
      double tentative_pos = mesh->nodes[i].pos.x - mesh->nodes[0].pos.x;
      mesh->centroid.x += (tentative_pos < origin.x - L0/2) ?
        tentative_pos + L0 : 
          ((tentative_pos > origin.x + L0/2) ? tentative_pos - L0 :
          tentative_pos);
    }
  foreach_dimension() 
    mesh->centroid.x = mesh->centroid.x/mesh->nln + mesh->nodes[0].pos.x;
  correct_node_pos(&mesh->centroid);
}

trace
void comp_volume(lagMesh* mesh) {
  coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
  comp_centroid(mesh);
  double volume = 0;
  for(int i=0; i<mesh->nlt; i++) {
    coord nodes[3];
    for(int j=0; j<3; j++)
      foreach_dimension() {
        double tentative_pos = mesh->nodes[mesh->triangles[i].node_ids[j]].pos.x
          - mesh->centroid.x;
        nodes[j].x = (tentative_pos < origin.x - L0/2) ?
          tentative_pos + L0 : 
          ((tentative_pos > origin.x + L0/2) ? tentative_pos - L0 :
          tentative_pos);
      }
    for(int j=0; j<3; j++) {
      coord cross_product;
      foreach_dimension() 
        cross_product.x = nodes[(j+1)%3].y*nodes[(j+2)%3].z - 
          nodes[(j+1)%3].z*nodes[(j+2)%3].y;
      volume += cdot(nodes[j],cross_product);
    }
  }
  mesh->volume = volume/18;
}


/** The function below updates the normal vectors on all the nodes as well as
the lengths and midpoints of all the edges (in 2D) or the area and centroids of
all the triangles (in 3D). */
void comp_normals(lagMesh* mesh) {
  if (!mesh->updated_normals) {
    #if dimension < 3
    compute_lengths(mesh);
    for(int i=0; i<mesh->nln; i++) {
      coord n[2];
      double l[2];
      double normn;
      for(int j=0; j<2; j++) {
        int edge_id;
        edge_id = mesh->nodes[i].edge_ids[j];
        l[j] = mesh->edges[edge_id].length;
        comp_edge_normal(mesh, edge_id);
        foreach_dimension() n[j].x = mesh->edges[edge_id].normal.x;
      }
      /** the normal vector at a node is the weighted average of the normal
      vectors of its edges. The average is weighted by the distance of the node
      to each of the edges' centers. */
      double epsilon = l[1]/(l[0] + l[1]);
      normn = 0.;
      foreach_dimension() {
        mesh->nodes[i].normal.x = epsilon*n[0].x + (1. - epsilon)*n[1].x;
        normn += sq(mesh->nodes[i].normal.x);
      }
      normn = sqrt(normn);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #else // dimension == 3
    comp_triangle_area_normals(mesh);
    for(int i=0; i<mesh->nln; i++) {
      foreach_dimension() mesh->nodes[i].normal.x = 0.;
      double sw = 0.; // sum of the weights
      for(int j=0; j<mesh->nodes[i].nb_triangles; j++) {
        int tid = mesh->nodes[i].triangle_ids[j];
        double dist = 0.;
        foreach_dimension()
          dist += sq(mesh->nodes[i].pos.x - mesh->triangles[tid].centroid.x);
        dist = sqrt(dist);
        sw += dist;
        foreach_dimension()
          mesh->nodes[i].normal.x += mesh->triangles[tid].normal.x*dist;
      }
      foreach_dimension()
        mesh->nodes[i].normal.x /= sw;
      double normn = cnorm(mesh->nodes[i].normal);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #endif
    mesh->updated_normals = true;
  }
}

#if dimension > 2

/**
## Useful macros

The macros below are useful to define an icosahedron
*/
#define GET_LD(NODE) ((fabs(fabs(NODE.pos.x) - ll) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - ll) < 1.e-8 ? 1 : 2)))
#define GET_LD_SIGN(NODE) ((GET_LD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_LD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_SD(NODE) ((fabs(fabs(NODE.pos.x) - sl) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - sl) < 1.e-8 ? 1 : 2)))
#define GET_SD_SIGN(NODE) ((GET_SD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_SD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_ZD(NODE) ((fabs(NODE.pos.x) < 1.e-8) ? 0 : \
  ((fabs(NODE.pos.y) < 1.e-8 ? 1 : 2)))

/**
## Operations on edges

The function below returns true if the two nodes $i$ and $j$ are neighbors
*/
bool is_neighbor(lagMesh* mesh, int i, int j) {
  for(int k=0; k<mesh->nodes[i].nb_neighbors; k++) {
    if (mesh->nodes[i].neighbor_ids[k] == j) return true;
  }
  return false;
}

/** The function below returns true if there is an edge connecting nodes i and
j */
bool edge_exists(lagMesh* mesh, int j, int k) {
  for(int i=0; i<mesh->nle; i++) {
    if ((mesh->edges[i].node_ids[0] == j && mesh->edges[i].node_ids[1] == k)
      || (mesh->edges[i].node_ids[0] == k
      && mesh->edges[i].node_ids[1] == j)) return true;
  }
  return false;
}

/** The function below returns true if the edge $i$ is across a priodic
boundary. */
bool is_edge_across_periodic(lagMesh* mesh, int i) {
  int n[2];
  for(int k=0; k<2; k++) n[k] = mesh->edges[i].node_ids[k];
  if (ACROSS_PERIODIC(mesh->nodes[n[0]].pos.x, mesh->nodes[n[1]].pos.x)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.y, mesh->nodes[n[1]].pos.y)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.z, mesh->nodes[n[1]].pos.z))
      return true;
  return false;
}

// foreach_dimension()
// bool is_edge_across_periodic_x(lagMesh* mesh, int i) {
//   int n[2];
//   for(int k=0; k<2; k++) n[k] = mesh->edges[i].node_ids[k];
//     if (ACROSS_PERIODIC(mesh->nodes[n[0]].pos.x, mesh->nodes[n[1]].pos.x))
//       return true;
//   return false;
// }


struct _write_edge {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  bool new_mesh;
  bool overwrite;
};

/** The function below writes edge i, connecting nodes j and k. If the edge
exists, the function returns false (no edge creation), true otherwise (edge
creation) */
bool write_edge(struct _write_edge p) {
  lagMesh* mesh = p.mesh;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool new_mesh = (p.new_mesh) ? p.new_mesh : false;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && edge_exists(mesh, j, k)) return false;
  else {
    mesh->edges[i].node_ids[0] = j;
    mesh->edges[i].node_ids[1] = k;
    for(int ii=0; ii<2; ii++) mesh->edges[i].triangle_ids[ii] = -1;
    if (new_mesh) {
      mesh->nodes[j].neighbor_ids[mesh->nodes[j].nb_neighbors] = k;
      mesh->nodes[k].neighbor_ids[mesh->nodes[k].nb_neighbors] = j;
      mesh->nodes[j].edge_ids[mesh->nodes[j].nb_neighbors] = i;
      mesh->nodes[k].edge_ids[mesh->nodes[k].nb_neighbors] = i;
      mesh->nodes[j].nb_neighbors++;
      mesh->nodes[k].nb_neighbors++;
    }
    return true;
  }
}

/** The function below creates a new edge between nodes i and j, and updates the
connectivity information of its nodes (but not its triangles, since they
don't exist yet). */
void new_edge(lagMesh* mesh, int i, int j) {
  int eid = mesh->nle; // id of the new edge
  int nodes[2];
  nodes[0] = i; nodes[1] = j;
  for(int k=0; k<2; k++) {
    mesh->edges[eid].node_ids[k] = nodes[k];

    /** Add the edge id to the newly connected nodes */
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      if (mesh->nodes[nodes[k]].edge_ids[l] == -1) {
        mesh->nodes[nodes[k]].edge_ids[l] = eid;
        break;
      }
    }

    /** Update the neighbors' list of the newly connected nodes */
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      if (mesh->nodes[nodes[k]].neighbor_ids[l] == -1) {
        mesh->nodes[nodes[k]].neighbor_ids[l] = nodes[(k+1)%2];
        break;
      }
    }

    /** The newly created edge is not yet surrounded by any triangle */
    mesh->edges[eid].triangle_ids[k] = -1;
  }
  mesh->nle++;
}

/** The function below splits an edge in two smaller edges, creating a node
at its midpoint. */
void split_edge(lagMesh* mesh, int i) {
  int nid[2];
  for(int j=0; j<2; j++) nid[j] = mesh->edges[i].node_ids[j];

  /** Create new node */
  foreach_dimension()
    mesh->nodes[mesh->nln].pos.x =
      .5*(mesh->nodes[nid[0]].pos.x + mesh->nodes[nid[1]].pos.x);
  mesh->nodes[mesh->nln].nb_neighbors = 6;
  mesh->nodes[mesh->nln].nb_triangles = 6;
  mesh->nodes[mesh->nln].neighbor_ids[0] = nid[0];
  mesh->nodes[mesh->nln].neighbor_ids[1] = nid[1];
  mesh->nodes[mesh->nln].edge_ids[0] = i;
  mesh->nodes[mesh->nln].edge_ids[1] = mesh->nle;
  for(int j=0; j<6; j++) {
    mesh->nodes[mesh->nln].triangle_ids[j] = -1;
    if (j>1) {
      mesh->nodes[mesh->nln].neighbor_ids[j] = -1;
      mesh->nodes[mesh->nln].edge_ids[j] = -1;
    }
  }

  /** Create new edge and update current one */
  write_edge(mesh, i, nid[0], mesh->nln, overwrite = true);
  write_edge(mesh, mesh->nle, nid[1], mesh->nln);
  for (int j=0; j<2; j++) {
    mesh->edges[i].triangle_ids[j] = -1;
    mesh->edges[mesh->nle].triangle_ids[j] = -1;
  }

  /** Update node information: neighboring nodes, connecting edges */
  for(int j=0; j<mesh->nodes[nid[0]].nb_neighbors; j++)
    if (mesh->nodes[nid[0]].neighbor_ids[j] == nid[1])
      mesh->nodes[nid[0]].neighbor_ids[j] = mesh->nln;
  for(int j=0; j<mesh->nodes[nid[1]].nb_neighbors; j++)
    if (mesh->nodes[nid[1]].neighbor_ids[j] == nid[0])
      mesh->nodes[nid[1]].neighbor_ids[j] = mesh->nln;
  for(int j=0; j<mesh->nodes[nid[1]].nb_neighbors; j++)
    if (mesh->nodes[nid[1]].edge_ids[j] == i)
      mesh->nodes[nid[1]].edge_ids[j] = mesh->nle;

  mesh->nln++;
  mesh->nle++;
}

/**
## Operations on triangles

The function below returns true if the triangle connecting nodes i,j and k
already exists in the mesh. */
bool triangle_exists(lagMesh* mesh, int i, int j, int k) {
  for(int t=0; t<mesh->nlt; t++) {
    for(int a=0; a<3; a++) {
      if (mesh->triangles[t].node_ids[a] == i) {
        for(int b=0; b<3; b++) {
          if (b != a && mesh->triangles[t].node_ids[b] == j) {
            for(int c=0; c<3; c++) {
              if (c != a && c != b && mesh->triangles[t].node_ids[c] == k) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}

/** The function below returns true if the triangle_ids $i$ is across a periodic
boundary. */
bool is_triangle_across_periodic(lagMesh* mesh, int i) {
  int e[3];
  for(int k=0; k<3; k++) e[k] = mesh->triangles[i].edge_ids[k];
  if (is_edge_across_periodic(mesh, e[0])
    || is_edge_across_periodic(mesh, e[1])
    || is_edge_across_periodic(mesh, e[2]))
      return true;
  return false;
}

struct _write_triangle {
  lagMesh* mesh;
  int tid;
  int i;
  int j;
  int k;
  bool overwrite;
};

/** The function below writes a triangle at index location tid, connecting nodes
i, j and k. It updates the connectivity information of its nodes and edges.*/
bool write_triangle(struct _write_triangle p) {
  lagMesh* mesh = p.mesh;
  int tid = p.tid;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && triangle_exists(mesh, i, j, k)) return false;
  else {
    mesh->triangles[tid].node_ids[0] = i;
    mesh->triangles[tid].node_ids[1] = j;
    mesh->triangles[tid].node_ids[2] = k;
    int c = 0;
    for(int a=0; a<3; a++) {
      int va = mesh->triangles[tid].node_ids[a];
      int b=(a+1)%3;
      int vb = mesh->triangles[tid].node_ids[b];
      for(int m=0; m<mesh->nodes[va].nb_neighbors; m++) {
        for(int n=0; n<mesh->nodes[vb].nb_neighbors; n++) {
          if (mesh->nodes[va].edge_ids[m] == mesh->nodes[vb].edge_ids[n]) {
            mesh->triangles[tid].edge_ids[c] = mesh->nodes[va].edge_ids[m];
            c++;
            int p = (mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[0]
              == -1) ? 0 : 1;
            mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[p] = tid;
          }
        }
      }
    }
    mesh->nodes[i].triangle_ids[mesh->nodes[i].nb_triangles] = tid;
    mesh->nodes[j].triangle_ids[mesh->nodes[j].nb_triangles] = tid;
    mesh->nodes[k].triangle_ids[mesh->nodes[k].nb_triangles] = tid;
    mesh->nodes[i].nb_triangles++;
    mesh->nodes[j].nb_triangles++;
    mesh->nodes[k].nb_triangles++;
    return true;
  }
}

struct _new_triangle {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  int prev_tid;
};

/** The function below also writes a new triangle connecting nodes i, j and k,
but in the specific context of the refinement process. The variable $prev\_tid$
is used to update the connectivity information of the nodes and edges. */
void new_triangle(lagMesh* mesh, int i, int j, int k, int prev_tid) {
  assert((i < mesh->nln) && (j < mesh->nln) && (k < mesh->nln));
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;

  int tid = mesh->nlt;
  for(int k=0; k<3; k++) {
    mesh->triangles[tid].node_ids[k] = nodes[k];
    /** Specify the id of the new triangle for node k */
    bool replaced_tid = false;
    for(int l=0; l<mesh->nodes[nodes[k]].nb_triangles; l++)
      if (mesh->nodes[nodes[k]].triangle_ids[l] == prev_tid) {
        mesh->nodes[nodes[k]].triangle_ids[l] = tid;
        replaced_tid = true;
        break;
      }
    if (!replaced_tid) {
      for(int l=0; l<mesh->nodes[nodes[k]].nb_triangles; l++)
        if (mesh->nodes[nodes[k]].triangle_ids[l] == -1) {
          mesh->nodes[nodes[k]].triangle_ids[l] = tid;
          replaced_tid = true;
          break;
        }
    }
    /** Find the edges and: (i) add them to the list of edges of the new
    triangle; (ii) specify the id of the new triangle for the edges */
    /** First, we find the edge that connects node i to node i+1 */
    int ce = -1; // ce for "current edge
    for(int l=0; l<mesh->nodes[nodes[k]].nb_neighbors; l++) {
      ce = mesh->nodes[nodes[k]].edge_ids[l];
      int cn = (mesh->edges[ce].node_ids[0] == nodes[k]) ?
        mesh->edges[ce].node_ids[1] : mesh->edges[ce].node_ids[0];
      if (cn == nodes[(k+1)%3]) break;
    }
    /** Then, we add this edge to the list of edges of our new triangle */
    mesh->triangles[tid].edge_ids[k] = ce;
    /** And we have to update the triangle id of the edge */
    replaced_tid = false;
    for(int l=0; l<2; l++)
      if (mesh->edges[ce].triangle_ids[l] == prev_tid) {
        mesh->edges[ce].triangle_ids[l] = tid;
        replaced_tid = true;
        break;
      }
    if (!replaced_tid && mesh->edges[ce].triangle_ids[0] != tid) {
      for(int l=0; l<2; l++)
        if (mesh->edges[ce].triangle_ids[l] == -1) {
          mesh->edges[ce].triangle_ids[l] = tid;
          replaced_tid = true;
          break;
        }
    }
  }
  mesh->nlt++;
}

void overwrite_triangle(lagMesh* mesh, int tid, int i, int j, int k) {
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;
  for(int i=0; i<3; i++) {
    /** 1. Take care of the nodes of the triangle */
    mesh->triangles[tid].node_ids[i] = nodes[i];
    /** Update the list of triangles for the node $i$ */
    bool already_there = false;
    for(int j=0; j<mesh->nodes[nodes[i]].nb_triangles; j++) {
      if (mesh->nodes[nodes[i]].triangle_ids[j] == tid) already_there = true;
      else if (mesh->nodes[nodes[i]].triangle_ids[j] == -1 && !already_there) {
        mesh->nodes[nodes[i]].triangle_ids[j] = tid;
        break;
      }
    }

    /** 2. Take care of the egdes of the triangle */
    /** 2.1. Identify the edge id connecting nodes $i$, $i+1$ */
    int eid = -1; // eid for "edge id"
    for(int j=0; j<mesh->nodes[nodes[i]].nb_neighbors; j++) {
      eid = mesh->nodes[nodes[i]].edge_ids[j];
      int cn = (mesh->edges[eid].node_ids[0] == nodes[i]) ?
        mesh->edges[eid].node_ids[1] : mesh->edges[eid].node_ids[0];
      if (cn == nodes[(i+1)%3]) break;
    }
    /** 2.2 Add this edge to the triangle's list of edges */
    mesh->triangles[tid].edge_ids[i] = eid;
    /** 2.3 Add the triangle id to the edge's list of triangles */
    int index = (mesh->edges[eid].triangle_ids[0] > -1 &&
      mesh->edges[eid].triangle_ids[0] != tid) ? 1 : 0;
    mesh->edges[eid].triangle_ids[index] = tid;
  }
}

/** The function below returns true if node j is a vertex of triangle i */
bool is_triangle_vertex(lagMesh* mesh, int i, int j) {
  for(int k=0; k<3; k++) {
    if (mesh->triangles[i].node_ids[k] == j) return true;
  }
  return false;
}

/**
The function below returns the (positive) angle between the two vectors formed
by the nodes [*n1, *n2] and [*n1, *n3].
*/
double comp_angle(lagNode* n1, lagNode* n2, lagNode* n3) {
  double theta = 0.;
  foreach_dimension() {
    theta += GENERAL_1DIST(n1->pos.x, n2->pos.x)*
             GENERAL_1DIST(n1->pos.x, n3->pos.x);
  }
  double norm1 = GENERAL_SQNORM(n1->pos, n2->pos);
  double norm2 = GENERAL_SQNORM(n1->pos, n3->pos);
  theta /= sqrt(norm1)*sqrt(norm2);
  theta = acos(theta);
  return theta;
}

/**
The function below returns true if the angle of node $i$ ($0 < i < 2$) of
triangle $tid$ is greater than $\pi$ radians.
*/
bool is_obtuse_node(lagMesh* mesh, int tid, int i) {
  lagNode* n[3];
  for(int j=0; j<3; j++)
    n[j] = &(mesh->nodes[mesh->triangles[tid].node_ids[j]]);
  if (comp_angle(n[i], n[(i+1)%3], n[(i+2)%3]) > pi/2.) return true;
  else return false;
}

/**
The function below returns true if the triangle $tid$ is obtuse at any of its
three angles.
*/
bool is_obtuse_triangle(lagMesh* mesh, int tid) {
  for(int i=0; i<3; i++) if (is_obtuse_node(mesh, tid, i)) return true;
  return false;
}

/**
## Uniform refinement of a mesh by subdividing its triangles

The function below loops through all triangles in the mesh and divide them
in four smaller ones. Keeping the correct structure of the mesh, i.e. updating
the relationship between nodes, edges, triangles and their neighbors results
in the somewhat complicated implementation below. */
void refine_mesh(lagMesh* mesh) {
  /** Perform the loop subdivision algorithm: until we reach the desired number
  of nodes, we split each triangles into four smaller ones */
  int cnt = mesh->nlt; // current number of triangles
  for(int i=0; i<cnt; i++) {
    int mid_ids[3];
    /** If not done yet, we split each edge into two */
    for(int j=0; j<3; j++) {
      int edge_id = mesh->triangles[i].edge_ids[j];
      if (mesh->edges[edge_id].triangle_ids[1] > -1) {
        split_edge(mesh, edge_id);
        mid_ids[j] = mesh->nln-1;
      }
      else {
        mid_ids[j] = is_triangle_vertex(mesh, i,
          mesh->edges[edge_id].node_ids[0]) ?
          mesh->edges[edge_id].node_ids[1] :
          mesh->edges[edge_id].node_ids[0];
      }
    }

    /** Connect the three midpoints with edges, and create corner triangles */
    for(int j=0; j<3; j++) {
      /** create edges between midpoints */
      new_edge(mesh, mid_ids[j], mid_ids[(j+1)%3]);

      /** Create the corner triangle with the new edge */
      int corner_id = -1;
      for(int k=0; k<3; k++) { // Loop over the current triangle vertices
        corner_id = mesh->triangles[i].node_ids[k];
        if (is_neighbor(mesh, corner_id, mid_ids[j]) &&
          is_neighbor(mesh, corner_id, mid_ids[(j+1)%3])) {
          break;
        }
      }
      new_triangle(mesh, mid_ids[j], mid_ids[(j+1)%3], corner_id, i);
    }
    /** Shrink the original (big) triangle into the center smaller one */
    overwrite_triangle(mesh, i, mid_ids[0], mid_ids[1], mid_ids[2]);
  }
}


/**
## Periodicity helper functions

In some situations, for instance to compute the volume of a capsule, we need to 
take the dot and cross products of
neighboring nodes that can be across periodic boundaries. The next three 
functions implement ``periodic-friendly" versions of the dot and cross products
that do not take into account the coordinates jump across the periodic
boundaries. The implementation relies on the assumption that a capsule is
*always* smaller in the x, y and z directions that half the domain size L0/2.
*/

/**
The function below corrects the coordinates of one node $a$ in order to
ensure it is placed on the same side of a reference coordinate (in practice,
the centroid of the capsule). The function returns the corrected node 
coordinate.
*/
coord correct_periodic_node_pos(coord a, coord ref) {
    coord result;
    foreach_dimension() {
        result.x = (fabs(a.x - ref.x) < L0/2) ? a.x : 
            (a.x > ref.x) ? a.x - L0 : a.x + L0;
    }
    return result;
}

/**
The function below corrects the coordinates of two nodes $a$ and $b$ in order to
ensure they are placed on the same side of a reference coordinate (in practice,
the centroid of the capsule). The result is stored in an array of `coord` of
length 2.
*/
void correct_periodic_nodes_pos(coord* result, coord a, coord b, coord ref) {
  foreach_dimension() {
    result[0].x = (fabs(a.x - ref.x) < L0/2) ? a.x : 
      (a.x > ref.x) ? a.x - L0 : a.x + L0;
    result[1].x = (fabs(b.x - ref.x) < L0/2) ? b.x : 
      (b.x > ref.x) ? b.x - L0 : b.x + L0;
  }
}

/**
The function below computes the cross product of two coordinates $a$ and $b$
that potentially lie across periodic boundaries. The coordinates are 
temporarilly moved on the same side of the periodic boundary as a reference 
coordinate `ref`, in practice the centroid of the capsule.
*/
foreach_dimension()
double periodic_friendly_cross_product_x(coord a, coord b, coord ref) {
  coord nodes[2];
  correct_periodic_nodes_pos(nodes, a, b, ref);
  return nodes[0].y*nodes[1].z - nodes[0].z*nodes[1].y;
}

/**
The function below computes the dot product of two coordinates $a$ and $b$
that potentially lie across periodic boundaries. The coordinates are 
temporarilly moved on the same side of the periodic boundary as a reference 
coordinate `ref`, in practice the centroid of the capsule.
*/
double periodic_friendly_dot_product(coord a, coord b, coord ref) {
  coord nodes[2];
  correct_periodic_nodes_pos(nodes, a, b, ref);
  return cdot(nodes[0], nodes[1]);
}


#endif // dimension > 2