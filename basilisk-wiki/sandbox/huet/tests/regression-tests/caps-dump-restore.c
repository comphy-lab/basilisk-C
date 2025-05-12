/**
# Test case for the dump and restore of a capsule

We create a capsule, dump and restore it, and ensure the restored version is identical to the dumped one.

*/
#define L0 1.
#define RADIUS .2
#define LEVEL 5
#define LAG_LEVEL 3
#define NCAPS 2

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"

int main(int argc, char* argv[]) {
  origin(-.50001*L0, -.50001*L0, -.50001*L0);
  N = 1 << LEVEL;
  run();
}

event init (i = 0) {
  activate_spherical_capsule(&CAPS(0), radius = RADIUS, 
    level = LAG_LEVEL);

  comp_volume(&CAPS(0));
  CAPS(0).initial_volume = CAPS(0).volume;

  double c0, c1, c2;
  c0 = 0.2072; c1 = 2.0026; c2 = -1.1228;
  double radius = RADIUS;
  for(int i=0; i<CAPS(0).nln; i++) {
    double rho = sqrt(sq(CAPS(0).nodes[i].pos.x) +
    sq(CAPS(0).nodes[i].pos.z))/radius;
    rho = (rho > 1) ? 1 : rho;
    int sign = (CAPS(0).nodes[i].pos.y > 0.) ? 1 : -1;
    CAPS(0).nodes[i].pos.y = sign*.5*radius*sqrt(1 - sq(rho))*
    (c0 + c1*sq(rho) + c2*sq(sq(rho)));
  }
  correct_lag_pos(&CAPS(0));
  comp_normals(&CAPS(0));
  comp_centroid(&CAPS(0));
  comp_volume(&CAPS(0));
  store_initial_configuration(&CAPS(0));
  
  /** 
  We need to initialize quantities such as forces and stretches in order to 
  compare our results to a reference file in a deterministic way.
  */
  for(int i=0; i<CAPS(0).nln; i++) {
    CAPS(0).nodes[i].curv = 1.;
    CAPS(0).nodes[i].ref_curv = 2.;
    foreach_dimension() {
      CAPS(0).nodes[i].lagVel.x = 3.;
      CAPS(0).nodes[i].lagForce.x = 4.;
    }
  }
  for(int i=0; i<CAPS(0).nle; i++) {
    CAPS(0).edges[i].l0 = 5.;
    CAPS(0).edges[i].length = 6.;
    foreach_dimension() CAPS(0).edges[i].normal.x = 7.;
  }
  for(int i=0; i<CAPS(0).nlt; i++) {
    foreach_dimension() CAPS(0).triangles[i].refShape[0].x = 8.;
    foreach_dimension() CAPS(0).triangles[i].refShape[1].x = 9.;
    CAPS(0).triangles[i].stretch[0] = 10.;
    CAPS(0).triangles[i].stretch[1] = 11.;
    CAPS(0).triangles[i].tension[0] = 12.;
    CAPS(0).triangles[i].tension[1] = 13.;
  }

  FILE* file = fopen("caps0.dump", "w");
  dump_lagmesh(file, &CAPS(0));
  fclose(file);
  file = fopen("caps0.dump", "r");
  restore_lagmesh(file, &CAPS(1));
  fclose(file);

  /**
  We now compare the dumped and restored capsules and output their difference
  which should be close to zero.
  */
  fprintf(stderr, "%d\n", CAPS(0).nln - CAPS(1).nln);
  fprintf(stderr, "%d\n", CAPS(0).nle - CAPS(1).nle);
  fprintf(stderr, "%d\n", CAPS(0).nlt - CAPS(1).nlt);
  foreach_dimension()
    fprintf(stderr, "%.5g\n", CAPS(0).centroid.x - CAPS(1).centroid.x);
  fprintf(stderr, "%.5g\n", CAPS(0).initial_volume - CAPS(1).initial_volume);
  fprintf(stderr, "%.5g\n", CAPS(0).volume - CAPS(1).volume);
  fprintf(stderr, "%d\n", 
          (CAPS(0).updated_stretches^CAPS(1).updated_stretches));
  fprintf(stderr, "%d\n", (CAPS(0).updated_normals^CAPS(1).updated_normals));
  fprintf(stderr, "%d\n", 
          (CAPS(0).updated_curvatures^CAPS(1).updated_curvatures));
  fprintf(stderr, "%d\n", (CAPS(0).isactive^CAPS(1).isactive));

  fprintf(stderr, "\n\t--- nodes ---\n\n");
  for(int i=0; i<CAPS(0).nln; i++) {
    foreach_dimension() {
      fprintf(stderr, "%.5g ", 
              CAPS(0).nodes[i].pos.x - CAPS(1).nodes[i].pos.x);
      fprintf(stderr, "%.5g ",
              CAPS(0).nodes[i].lagVel.x - CAPS(1).nodes[i].lagVel.x);
      fprintf(stderr, "%.5g ",
              CAPS(0).nodes[i].normal.x - CAPS(1).nodes[i].normal.x);
      fprintf(stderr, "%.5g ",
              CAPS(0).nodes[i].lagForce.x - CAPS(1).nodes[i].lagForce.x);
    }
    fprintf(stderr, "%.5g ", CAPS(0).nodes[i].curv - CAPS(1).nodes[i].curv);
    fprintf(stderr, "%.5g ", CAPS(0).nodes[i].gcurv - CAPS(1).nodes[i].gcurv);
    fprintf(stderr, "%.5g ", 
            CAPS(0).nodes[i].ref_curv - CAPS(1).nodes[i].ref_curv);
    fprintf(stderr, "%d ", 
            CAPS(0).nodes[i].nb_neighbors - CAPS(1).nodes[i].nb_neighbors);
    fprintf(stderr, "%d ", 
            CAPS(0).nodes[i].nb_triangles - CAPS(1).nodes[i].nb_triangles);
    for(int j=0; j<CAPS(1).nodes[i].nb_neighbors; j++) {
      fprintf(stderr, "%d ", 
        CAPS(0).nodes[i].neighbor_ids[j] - CAPS(1).nodes[i].neighbor_ids[j]);
      fprintf(stderr, "%d ", 
        CAPS(0).nodes[i].edge_ids[j] - CAPS(1).nodes[i].edge_ids[j]);
    }
    for(int j=0; j<CAPS(1).nodes[i].nb_triangles; j++) {
      fprintf(stderr, "%d ", 
        CAPS(0).nodes[i].triangle_ids[j] - CAPS(1).nodes[i].triangle_ids[j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n\t--- edges ---\n\n");
  for(int i=0; i<CAPS(0).nle; i++) {
    fprintf(stderr, "%d ", 
      CAPS(0).edges[i].node_ids[0] - CAPS(1).edges[i].node_ids[0]);
    fprintf(stderr, "%d ", 
      CAPS(0).edges[i].node_ids[1] - CAPS(1).edges[i].node_ids[1]);
    fprintf(stderr, "%d ", 
      CAPS(0).edges[i].triangle_ids[0] - CAPS(1).edges[i].triangle_ids[0]);
    fprintf(stderr, "%d ", 
      CAPS(0).edges[i].triangle_ids[1] - CAPS(1).edges[i].triangle_ids[1]);
    fprintf(stderr, "%d ", 
      CAPS(0).edges[i].node_ids[0] - CAPS(1).edges[i].node_ids[0]);
    fprintf(stderr, "%.5g ", CAPS(0).edges[i].l0 - CAPS(1).edges[i].l0);
    fprintf(stderr, "%.5g ", 
      CAPS(0).edges[i].length - CAPS(1).edges[i].length);
    foreach_dimension()
      fprintf(stderr, "%.5g ", 
              CAPS(0).edges[i].normal.x - CAPS(1).edges[i].normal.x);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n\t--- triangles ---\n\n");
  for(int i=0; i<CAPS(0).nlt; i++) {
    for(int j=0; j<3; j++) {
      fprintf(stderr, "%d ", 
        CAPS(0).triangles[i].node_ids[j] - CAPS(1).triangles[i].node_ids[j]);
      fprintf(stderr, "%d ", 
        CAPS(0).triangles[i].edge_ids[j] - CAPS(1).triangles[i].edge_ids[j]);
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].sfc[j][0] - CAPS(1).triangles[i].sfc[j][0]);
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].sfc[j][1] - CAPS(1).triangles[i].sfc[j][1]);
    }
    foreach_dimension() {
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].normal.x - CAPS(1).triangles[i].normal.x);
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].centroid.x - CAPS(1).triangles[i].centroid.x);
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].refShape[0].x - 
        CAPS(1).triangles[i].refShape[0].x);
      fprintf(stderr, "%.5g ", 
        CAPS(0).triangles[i].refShape[1].x - 
        CAPS(1).triangles[i].refShape[1].x);
    }
    fprintf(stderr, "%.5g ", 
      CAPS(0).triangles[i].area - CAPS(1).triangles[i].area);
    fprintf(stderr, "%.5g ", 
      CAPS(0).triangles[i].stretch[0] - CAPS(1).triangles[i].stretch[0]);
    fprintf(stderr, "%.5g ", 
      CAPS(0).triangles[i].stretch[1] - CAPS(1).triangles[i].stretch[1]);
    fprintf(stderr, "%.5g ", 
      CAPS(0).triangles[i].tension[0] - CAPS(1).triangles[i].tension[0]);
    fprintf(stderr, "%.5g ", 
      CAPS(0).triangles[i].tension[1] - CAPS(1).triangles[i].tension[1]);
    fprintf(stderr, "\n");
  }

  exit(0);
}

