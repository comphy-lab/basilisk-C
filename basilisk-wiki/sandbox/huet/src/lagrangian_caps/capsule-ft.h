/**
# Triangulated mesh of a capsule

In this file, we implement a Lagrangian mesh for the front-tracking method,
meant to track the position and compute the stresses of an elasitc membrane.
*/

#ifndef dimension
  #define dimension 3
#endif

/**
## Structure of the mesh

In the Lagrangian mesh, each node is assigned several attributes:

* ```pos```, the node coordinates
* ```lagVel```, the node velocity vector
* ```normal```, the vector normal to the membrane at the node coordinates
* ```curv```, the (mean) curvature of the membrane at the node coordinates
* ```ref_curv```, the reference (mean) curvature of the membrane at the node coordinates (default is zero)
* ```lagForce```, an elastic and/or bending force exerted by the membrane on its surrounding fluid
* ```stencil```, a [Cache](http://basilisk.fr/src/grid/tree.h#82) structure used for averaging the neighboring velocities and spreading the Lagrangian force as a body force on the neighboring Eulerian nodes
* in case of MPI simulations, ```pid```, the rank of the processor owning the Eulerian cell which contains the lagNode
* ```edge_ids```, the IDs of its connecting edges: 2 in 2D, up to 6 in 3D (because every considered shape is derived by subdividing in icosahedron, leading to 5 or 6 neighbors only).

In case of 3D simulations, other attributes are introduced:

* ```nb_neighbors```, the number of neighbors of the node, which should only be 5 or 6
* ```neighbor_ids```, the IDs of the node neighbors
* ```nb_triangles```, the number of triangles having the node as a vertex
* ```triangle_ids```, the ids of the above triangles
* ```gcurv```, the Gaussian curvature of the membrane at the node coordinates
* ```nb_fit_iterations```, the number of iterations needed to compute the membrane curvature and normal vector to the desired convergence threshold.

*/
typedef struct lagNode {
  coord pos;
  coord lagVel;
  coord normal;
  double curv;
  double ref_curv;
  coord lagForce;
  Cache stencil;
  #if _MPI
    int pid;
  #endif
  #if dimension < 3
    int edge_ids[2];
  #else
    int nb_neighbors;
    int neighbor_ids[6];
    int edge_ids[6];
    int nb_triangles;
    int triangle_ids[6];
    double gcurv;
    int nb_fit_iterations;
  #endif
} lagNode;

/** We specify the size of the 5x5(x5) stencil in 2D or 3D. */
#if dimension < 3
  #define STENCIL_SIZE 25
#else
  #define STENCIL_SIZE 125
#endif

/** Similarly, the edges of the mesh are assigned:

* ```node_ids```, the IDs of the two nodes they connect
* In case of 3D simulations, ```triangle_ids```, the IDs of the two triangles they separate
* ```l0```, the length of the edge in the initial stress-free configuration
* ```length```, the current length of the edge
* ```normal```, a vector normal to the membrane at the edge midpoint.

*/
typedef struct Edge {
  int node_ids[2];
  #if dimension > 2
    int triangle_ids[2];
  #endif
  double l0;
  double length;
  coord normal;
} Edge;

/** In case of 3D simulations, we define triangle faces. Each ```Triangle``` structure has the following attributes:

* ```node_ids```, the IDs of the triangle vertices
* ```edge_ids```, the IDs of the triangle edges
* ```area```, the current area of the triangle
* ```normal```, the normal vector to the triangle pointing outside the membrane
* ```centroid```, the coordinates of the triangle centroid
* ```refShape```, the coordinates of the vertices in the Common Plane in the stress-free configuration. By convention the first vertex is always placed at $(0, 0, 0)$, so only the coordinates of the second and third vertex are stored in ```refShape```
* ```sfc```, the shape function coefficients used in the finite elements method and computed in [store_initial_configuration](elasticity-ft.h#store_initial_configuration). There are two coefficients per vertex, hence the $3 \times 2$ array structure of ```sfc```.

 */
#if dimension > 2
  typedef struct Triangle {
    int node_ids[3];
    int edge_ids[3];
    double area;
    coord normal;
    coord centroid;
    coord refShape[2];
    double sfc[3][2]; // sfc for "shape function coefficients"
    double stretch[2];
    double tension[2];
  } Triangle;
#endif

/** The ```lagMesh``` structure contains arrays of the previously introduced nodes, edges and triangles. It defines an unstructured mesh, the membrane of our capsule. Its attributes are:

*  ```nln``` the number of Lagrangian nodes
* ```nodes```, the array of Lagrangian nodes
* ```nle```, the number of Lagrangian edges
* ```edges```, the array of Lagrangian edges
* In case of 3D simulations:
    * ```nlt```, the number of Lagrangian triangles
    * ```triangles```, the array of Lagrangian triangles
* ```updated_stretches```, a boolean used to check if the current length of the edges has been updated since the last advection of the Lagrangian nodes
* ```updated_normals```, a similar boolean telling if the nodal normal vectors should be recomputed
* ```updated_curvatures```, a last boolean telling if the nodal curvatures should be recomputed.
* ```isactive```, a boolean indicating if the capsule exists in the flow (useful
when capsules are introduced during a simulation)
*/

typedef struct lagMesh {
  int nln;
  lagNode* nodes;
  int nle;
  Edge* edges;
  #if dimension > 2
    int nlt;
    Triangle* triangles;
  #endif
  coord centroid;
  double initial_volume;
  double volume;
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
  bool isactive;
} lagMesh;

/** We denote by ```NCAPS``` the number of Lagrangian meshes, or capsules, in
the simulation. It is one by default. */
#ifndef NCAPS
  #define NCAPS 1
#endif
#ifndef RESTART_CASE
  #define RESTART_CASE 0
#endif

/** The Lagrangian mesh is accessible in the code thanks to the structure
below, which is simply an array of Lagrangian meshes (useful when several of
them are considered). The macro $CAPS(k)$ can be used as a shortcut to access the
$k^{th}$ membrane. */
typedef struct Capsules {
  lagMesh caps[NCAPS];
  int nbcaps;
} Capsules;
Capsules allCaps;
#define CAPS(i) (allCaps.caps[i])


/**
## Initialization, memory management and useful macros.
*/
void initialize_empty_capsule(lagMesh* mesh) {
  mesh->nln = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  #if dimension > 2
    mesh->nlt = 0;
    mesh->triangles = NULL;
  #endif
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
  mesh->isactive = false;
}

void free_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) free(mesh->nodes[i].stencil.p);
  free(mesh->nodes);
  free(mesh->edges);
  #if dimension > 2
    free(mesh->triangles);
  #endif
}

void free_all_caps(Capsules* caps) {
  for(int i=0; i<caps->nbcaps; i++)
    if (CAPS(i).isactive)
      free_one_caps(&(caps->caps[i]));
}

/**
## Adding capsules after restart

In case we add membranes to a simulation after restart, we need to call
the function below. */
void initialize_capsules() {
  allCaps.nbcaps = NCAPS;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) initialize_empty_capsule(&CAPS(i));
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

void initialize_capsule_stencils(lagMesh* mesh) {
  for(int j=0; j<mesh->nln; j++) {
    mesh->nodes[j].stencil.n = STENCIL_SIZE;
    mesh->nodes[j].stencil.nm = STENCIL_SIZE;
    mesh->nodes[j].stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));
  }
}

void initialize_all_capsules_stencils() {
  for(int i=0; i<NCAPS; i++)
    if (CAPS(i).isactive)
      initialize_capsule_stencils(&CAPS(i));
}

void initialize_active_capsule(lagMesh* mesh) {
  initialize_empty_capsule(mesh);
  mesh->isactive = true;
  initialize_capsule_stencils(mesh);
}

/** The next few macros are useful to compute signed distances and averages
across periodic boundaries. We assume for this purpose that the length of
the edges are less that half the domain size, which in practice should always
be the case. */
#define ACROSS_PERIODIC(a,b) (fabs(a - b) > L0/2.)
#define PERIODIC_1DIST(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DIST(a,b) : a - b)
#define PERIODIC_1DAVG(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DAVG(a,b) : a + b)
#define GENERAL_SQNORM(a,b) (sq(GENERAL_1DIST(a.x, b.x)) + \
  sq(GENERAL_1DIST(a.y, b.y)) + sq(GENERAL_1DIST(a.z, b.z)))

#if dimension < 3
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y)
#else
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y + a.z*b.z)
#endif

/**
## Advection of the mesh

The advection have to be compatible with MPI, so we start by including the MPI
communication functions. We also include [reg-dirac.h](reg-dirac.h), which
implements the interpolation of the Eulerian velocities onto the Lagrangian
nodes.
*/

/** By default, the mesh is advected using a second-order two-step Runge Kutta
scheme. If the following macro is set to 0, a first-order forward Euler schme
is used instead. */
#ifndef ADVECT_LAG_RK2
  #define ADVECT_LAG_RK2 1
#endif

#if _MPI
  #include "capsule-ft-mpi.h"
#endif
#include "ibm-ft.h"
#ifndef CONSERVE_VOLUME
  #define CONSERVE_VOLUME 1
#endif

/**
A lot of triangulated mesh-related functions are defined in a separate header
file.
*/
#include "geometry-ft.h"

#if CONSERVE_VOLUME
  #include "volume-conservation-ft.h"
#endif

/**
The function below advects each Lagrangian node by
interpolating the velocities around the node of interest. By default, a
second-order Runge Kutta scheme is used. By setting the macro
```ADVECT_LAG_RK2``` to 0, a simple forward Euler scheme is used.
*/
trace
void advect_lagMesh(lagMesh* mesh) {
  eul2lag(mesh);
  #if !(ADVECT_LAG_RK2)
    for(int i=0; i < mesh->nln; i++) {
      foreach_dimension() {
        mesh->nodes[i].pos.x += dt*mesh->nodes[i].lagVel.x;
      }
    }
  #else
    lagMesh buffer_mesh;
    buffer_mesh.isactive = true;
    buffer_mesh.nln = mesh->nln;
    buffer_mesh.nodes = malloc(mesh->nln*sizeof(lagNode));
    for(int i=0; i<mesh->nln; i++) {
      // Step 1 of RK2
      foreach_dimension()
        buffer_mesh.nodes[i].pos.x = mesh->nodes[i].pos.x +
          .5*dt*mesh->nodes[i].lagVel.x;
    }
    correct_lag_pos(&buffer_mesh);
    for(int j=0; j<buffer_mesh.nln; j++) {
      buffer_mesh.nodes[j].stencil.n = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.nm = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.p = malloc(STENCIL_SIZE*sizeof(Index));
    }
    generate_lag_stencils_one_caps(&buffer_mesh);
    eul2lag(&buffer_mesh);
    for(int i=0; i<mesh->nln; i++) {
      // Step 2 of RK2
      foreach_dimension()
        mesh->nodes[i].pos.x += dt*buffer_mesh.nodes[i].lagVel.x;
    }
    for(int i=0; i<buffer_mesh.nln; i++) free(buffer_mesh.nodes[i].stencil.p);
    free(buffer_mesh.nodes);
  #endif
  correct_lag_pos(mesh);
  #if CONSERVE_VOLUME
    enforce_optimal_volume_conservation(mesh);
  #endif
  comp_centroid(mesh);
  comp_volume(mesh);
  generate_lag_stencils_one_caps(mesh);
}


/**
## Putting the pieces together

Below, we call the above functions at the appropriate time using the Basilisk
event syntax.

We start by creating empty Lagrangian meshes, and allocating an acceleration
field in case it isn't done yet by another Basilisk solver.
*/
event defaults (i = 0) {
  allCaps.nbcaps = NCAPS;
  // #if (RESTART_CASE == 0)
    for(int i=0; i<NCAPS; i++) {
      initialize_empty_capsule(&CAPS(i));
    }
  // #endif
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also use this loop as an opportunity to
re-initialize the Lagrangian forces to zero. */
event tracer_advection(i++) {
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) {
      advect_lagMesh(&CAPS(i));
      for(int j=0; j<CAPS(i).nln; j++)
        foreach_dimension() CAPS(i).nodes[j].lagForce.x = 0.;
    }
  }
}

/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */
vector forcing[];
event acceleration (i++) {
  face vector ae = a;
  foreach()
    if (cm[] > 1.e-20) foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) lag2eul(forcing, &CAPS(i));
  }
  foreach_face()
    if (fm.x[] > 1.e-20) ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  free_all_caps(&allCaps);
}

/**
## Additional functionalities
*/
#if dimension > 2
  #include "dump-ft.h"
  #include "post-processing-ft.h"
#endif

/**
## Tests
[advect_caps.c](../../tests/lagrangian_caps/advect_caps.c): Tests the
convergence of the advection scheme.


[curvature.c](../../tests/lagrangian_caps/curvature.c): Tests the computation of
the curvature at the Lagrangian nodes. Since the curvature depends on the
normals, this case also validates the computation of the normal vectors.
*/
