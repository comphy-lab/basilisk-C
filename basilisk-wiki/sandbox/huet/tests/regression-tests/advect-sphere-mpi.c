/**
# Advection of a front-tracking mesh in MPI

In this file we test the advection of a capsule described by a Lagrangian mesh
across a periodic boundary using 4 processors.

This is a test for [capsule-ft-mpi.h](../../src/lagrangian_caps/capsule-ft-mpi.h). */

#define LEVEL 6
#define LAG_LEVEL 4
#define RADIUS .125
#define L0 1.
#define T_END 1.

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  init_grid(N);
  periodic(left);
  TOLERANCE = HUGE;
  stokes = true;
  DT = 5.e-3;
  run();
}

coord* ref_data = NULL;
event init (i = 0) {
  activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL, radius = RADIUS);
  ref_data = malloc(CAPS(0).nln*sizeof(coord));
  for(int i=0; i<CAPS(0).nln; i++)
    foreach_dimension() ref_data[i].x = CAPS(0).nodes[i].pos.x;
}

event impose_u (i++) {
  foreach() {
    double x0 = x + L0 -t*L0*T_END;
    double y0 = y + L0;
    double a = -2*sq(sin(pi*x0))*sin(pi*y0)*cos(pi*y0)*cos(pi*t/T_END);
    double b = -2*sin(pi*x0)*cos(pi*x0)*sq(cos(pi*y0))*cos(pi*t/T_END);
    double theta = pi/4.;
    u.x[] = a*cos(theta)/2 + L0*T_END;
    u.y[] = b/2;
    u.z[] = a*sin(theta)/2;
  }
}

event adapt (i++) {
  tag_ibm_stencils(&CAPS(0));
  adapt_wavelet({stencils}, (double []){1.e-2}, maxlevel = LEVEL);
  generate_lag_stencils(&CAPS(0));
}

/** We compute the time evolutions of the normalized area and volume 
of the capsule */
event progress_output (i++) {
  if (pid() == 0) {
    comp_triangle_area_normals(&CAPS(0));
    double narea = 0;
    for(int j=0; j<CAPS(0).nlt; j++) narea += CAPS(0).triangles[j].area;
    narea /= 4*pi*sq(RADIUS);
    comp_volume(&CAPS(0));
    double nvolume = CAPS(0).volume/CAPS(0).initial_volume;
    fprintf(stderr, "%d, %.5g, %.5g, %.5g, %.5g %.5g\n", i, t, narea, nvolume,
    CAPS(0).centroid.x, CAPS(0).initial_volume);
    fflush(stderr);
  }
}

event movie (i+=5) {
  view(fov = 25, bg = {1,1,1}, theta = 5*pi/6, psi = 0., phi = pi/8);
  clear();
  draw_lag(&CAPS(0), lw = .5, edges = true, facets = true);
  cells(n = {0,0,1});
  save("advected_sphere_mpi.mp4");
}

/** At the end of the simulation, we also compare the position of the
membrane to its initial position. With an asymptotically fine space
and time resolutions, they would coincide. */
event output (t = T_END) {
  if (pid() == 0) {
    double avg_err, max_err;
    avg_err = 0.; max_err = -HUGE;
    for(int i=0; i < CAPS(0).nln; i++){
    double err = 0.;
    foreach_dimension() err += sq(GENERAL_1DIST(ref_data[i].x,
        CAPS(0).nodes[i].pos.x));
    err = sqrt(err);
    avg_err += err;
    if (err > max_err) max_err = err;
    }
    avg_err /= CAPS(0).nln;
    fprintf(stderr, "%.3g, %.3g\n", avg_err, max_err);
    fflush(stderr);
  }
}

event end (t = T_END) {
  free(ref_data);
  return 0;
}

/**
##Results

![Movie of an advected sphere](advect-sphere-mpi/advected_sphere_mpi.mp4)
*/
