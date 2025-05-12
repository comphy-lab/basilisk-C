/**
# Capsule migration in a helicoidal channel

![Helicoidal channel geometry](https://damienhuet.github.io/images/basilisk_sandbox_files/helix/helix_mesh.png)(width="20%")

![Movie of the capsule in the helix at Re=100](https://damienhuet.github.io/images/basilisk_sandbox_files/helix/ux_cells.mp4)(width="15%")

We consider the inertial migration of an elastic capsule in an helicoidal
channel for three Reynolds number: 10, 50 and 100.

## How to run
This case needs to be compiled and run three times:
  - Once with the flag ```-DINIT_GEOM=1``` where the STL geometry is read and the mesh and volume fraction field are dumped in ```geometry_generation.dump```. This specific run has to be done on a single processor, and extra RAM memory is allocated (about 30GB are needed).
  - Once with the flag ```-DINITIAL_FLOW_FIELD=1``` where the flow field is developped for 20 iterations.
  - Once with none of the above flags where the capsule is inserted at the top of the helix and the physics happens.

## Source file
*/

#define L0 (13.)
#define HELIX_RADIUS (1.)
#define PIPE_RADIUS (.2)
#define RADIUS (.05)
#define U_AVG (20*RADIUS)
#ifndef RE
  #define RE 10.
#endif
#define RHO 1.
#define MU (RHO*U_AVG*(2*PIPE_RADIUS)/RE)
#define GRAD_P 100.
#ifndef CA
  #define CA 0.2
#endif
#define E_S (MU*U_AVG/CA)

#if INITIAL_FLOW_FIELD
  #define NCAPS 0
#else
  #if INIT_GEOM
    #define NCAPS 0
  #else
    #define NCAPS 1
  #endif
#endif

/** We define the following solver-ralated quantities:

* the minimum and maximum refinement levels of the Eulerian grid
* the non-dimensional duration time of the simulation TEND
* the number of Lagrangian points NLP on the capsule
* the maximum time-step DT_MAX, which appears to be Reynolds dependant
* the tolerance of the Poisson solver (maximum admissible residual)
* the tolerance of the wavelet adaptivity algorithm for the velocity
* the $stokes$ boolean, in order to ignore the convective term in the
Navier-Stokes solver
* the Jacobi preconditionner, which we switch on for this case.
* the output frequency: a picture is generated every OUTPUT_FREQ iterations
*/

#ifndef MINLEVEL
  #define MINLEVEL 2
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 11
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4
#endif
#ifndef DT_MAX
  #define DT_MAX 1.e-3
#endif
#ifndef END_COND
  #if INITIAL_FLOW_FIELD
    #define END_COND (i = NB_ITER_INITIAL_FLOW_FIELD)
  #else
    #define END_COND (t = 23./U_AVG)
  #endif
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-3
#endif
#ifndef U_TOL
  #define U_TOL (0.1*U_AVG)
#endif
#ifndef STOKES
  #define STOKES false
#endif
#define JACOBI 1
#define NB_ITER_INITIAL_FLOW_FIELD 20
#if (NCAPS > 0)
  #define INIT_GEOM 0
  #define INITIAL_FLOW_FIELD 0
  #define OUTPUT_FREQ 50
#else
  #define OUTPUT_FREQ 1
#endif


/**
## Simulation setup

We import the octree grid, the embedded boundaries, the centered Navier-Stokes
solver, two header files useful to read STL geometries, the Lagrangian
mesh, the neo-Hookean elastic law, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

scalar dist_embed[];
void distance_from_stl (scalar c, face vector f, FILE * fp, int maxlevel, int minlevel)
{
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1.e-30}, maxlevel = maxlevel, minlevel = minlevel).nf);

  foreach() dist_embed[] = d[];
}

FILE* foutput = NULL;
FILE* fperf = NULL;

scalar rhov[];
face vector muv[];
face vector alphav[];

int main(int argc, char* argv[]) {
  origin(-.5*L0 + 1.e-6, -.5*L0 + 1.e-6, -0.5);
  N = 1 << MINLEVEL;
  init_grid(N);
  mu = muv;
  alpha = alphav;
  rho = rhov;
  TOLERANCE = MY_TOLERANCE*U_AVG;
  TOLERANCE_MU = MY_TOLERANCE*MU;
  /** We don't need to compute the convective term in this case, so we set the
  boolean $stokes$ to false. However it is still important to choose $Re \ll 1$
  since we are solving the unsteady Stokes equation. */
  stokes = STOKES;
  DT = DT_MAX;

  run();
}
u.n[top] = dirichlet(-U_AVG);
u.n[bottom] = dirichlet(-U_AVG);
u.t[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[top] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.r[front] = dirichlet(0.);
u.r[back] = dirichlet(0.);
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

int restored_iteration = 0;
event init (i = 0) {
  #if NCAPS > 0
    foutput  = fopen("output.txt","w");
    fperf = fopen("perf.txt", "w");
  #else
    foutput  = fopen("output_initial_flow_field.txt","w");
    fperf = fopen("perf_initial_flow_field.txt", "w");
  #endif
  fprintf(fperf, "i, t, dt, total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");

  /** In this case we use a third order face-flux interpolation. */
  for (scalar s in {u})
    s.third = true;

  for (scalar s in {u, p, pf}) {
    s.gradient = minmod2;
  }

  #if INITIAL_FLOW_FIELD
    restore("geometry_generation.dump");
    restored_iteration = i;
    vertex scalar phi[];
    foreach_vertex()
      phi[] = (dist_embed[] + dist_embed[-1] + dist_embed[0,-1] +
        dist_embed[-1,-1] + dist_embed[0,0,-1] + dist_embed[-1,0,-1] +
        dist_embed[0,-1,-1] + dist_embed[-1,-1,-1])/8.;
    fractions(phi, cs, fs);
    fractions_cleanup(cs, fs, smin = 1.e-30);
  #else
    #if (NCAPS > 0)
    dist_embed.prolongation = refine_injection;
    restore("initial_flow_field.dump");
    initialize_capsules();
    restored_iteration = i;
    vertex scalar phi[];
    foreach_vertex()
      phi[] = (dist_embed[] + dist_embed[-1] + dist_embed[0,-1] +
        dist_embed[-1,-1] + dist_embed[0,0,-1] + dist_embed[-1,0,-1] +
        dist_embed[0,-1,-1] + dist_embed[-1,-1,-1])/8.;
    fractions(phi, cs, fs);
    fractions_cleanup(cs, fs, smin = 1.e-30);
    #endif
  #endif
}

#if INIT_GEOM
  event init_adapt (i = 0) {
#else
  #if INITIAL_FLOW_FIELD
    event init_adapt (i = 1) {
  #else
    event init_adapt (i = NB_ITER_INITIAL_FLOW_FIELD + 1) {
  #endif
#endif

  #if INIT_GEOM
    FILE* fp = fopen("helix2.stl", "r");
    distance_from_stl(cm, fm, fp, MAXLEVEL, MINLEVEL);
    fclose(fp);
    vertex scalar phi[];
    foreach_vertex()
      phi[] = (cs[] + cs[-1] + cs[0,-1] +
        cs[-1,-1] + cs[0,0,-1] + cs[-1,0,-1] +
        cs[0,-1,-1] + cs[-1,-1,-1])/8.;
    fractions(phi, cs, fs);
    fractions_cleanup(cs, fs, smin = 1.e-30);
  #else
    u.n[embed] = dirichlet(0.);
    u.t[embed] = dirichlet(0.);
    u.r[embed] = dirichlet(0.);

    #if (NCAPS > 0)
    /** We initialize a spherical membrane using the pre-defined function in
    [common-shapes-ft.h](../../src/lagrangian_caps/common-shapes-ft.h).*/
    activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL, radius = RADIUS);
    initialize_all_capsules_stencils();
    store_initial_configuration(&CAPS(0));

    for(int k=0; k<CAPS(0).nln; k++) {
      /** We place the capsule at the top of the helix, close to the inner
      (high-curvature) wall. */
      CAPS(0).nodes[k].pos.x += HELIX_RADIUS;
      CAPS(0).nodes[k].pos.y += 1.5;
      CAPS(0).nodes[k].pos.z += 12.;
    }
    generate_lag_stencils(no_warning = true);

    /** We refine the mesh around the embedded geometry and around the
    membrane*/
    astats ss;
    int ic = 0;
    do {
      ic++;
      tag_ibm_stencils();
      ss = adapt_wavelet ({stencils, cs, u}, (double[]) {1.e-30, 1.e-30,
        U_TOL, U_TOL, U_TOL}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
      generate_lag_stencils(no_warning = true);
    } while ((ss.nf || ss.nc) && ic < 100);

    #else
      foreach() foreach_dimension() u.x[] = 0.;
      foreach_face(x) uf.x[] = 0.;
    #endif
  #endif
}


/** Due to embedded boundaries, the viscosity and density fields are not
constant since they have to be multiplied by the volume fraction. */
event properties (i++) {
  foreach_face() {
    muv.x[] = MU*fm.x[];
    alphav.x[] = 1./RHO*fm.x[];
  }
  foreach() rhov[] = RHO*cm[];
}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm,
* second, we adapt according to the fluid volume fraction in order to refine
near the walls,
* lastly we adapt according to the velocity field.

When computing the initial flow field, there is no capsule and the adaptivity
crietrion is on the fluid volume fraction and the velocity only.
*/
event adapt (i++) {
  #if (NCAPS > 0)
    tag_ibm_stencils();
    adapt_wavelet({stencils, cs, u}, (double []){1.e-30,
      1.e-30, U_TOL, U_TOL, U_TOL}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    generate_lag_stencils();
  #else
    adapt_wavelet({cs, u}, (double []){1.e-30,
        U_TOL, U_TOL, U_TOL}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  #endif
}

event logperf (i++) {
  if (pid() == 0) {
    fprintf(fperf, "%d %g %g %ld %d %g %g %d %d %g %g %d\n", i, t, dt, grid->tn, mgu.i, mgu.resb, mgu.resa, mgu.nrelax, mgp.i, mgp.resb, mgp.resa, mgp.nrelax);
    fflush(fperf);
  }
}

/**
In the event below, we output the maximum length of the capsule in the x, y and
z directions, as well as the position of the centroid of the capsule.
*/
vector prev_u[];
event logfile (i += OUTPUT_FREQ) {
  #if INITIAL_FLOW_FIELD
    double max_deltau = -HUGE;
    double avg_deltau = 0;
    double flow_rate = 0.;
    int nbc = 0;
    foreach(reduction(max:max_deltau)
      reduction(+:avg_deltau) reduction(+:nbc) reduction(+:flow_rate)) {
      if (cm[] > 1.e-10) {
        nbc++;
        double nu = sqrt(sq(u.x[] - prev_u.x[]) + sq(u.y[] - prev_u.y[])
          + sq(u.z[] - prev_u.z[]));
        avg_deltau += nu;
        if (nu > max_deltau) max_deltau = nu;
        foreach_dimension() prev_u.x[] = u.x[];

        /* We also output the flow rate at the boundary of the domain */
        if (fabs(y - (L0 - Delta)/2.) < Delta/2.)
          flow_rate += fabs(u.y[])*sq(Delta)*fm.x[];
      }
    }
    if (nbc > 0) avg_deltau /= nbc;
    double avg_output_vel = flow_rate/(sq(PIPE_RADIUS)*3.1415);

    fprintf(foutput, "%d %d %g %g %g %g %g\n", i, nbc, avg_deltau, max_deltau, normf(u.x).max, flow_rate, avg_output_vel);
    fflush(foutput);
  #endif

  #if (NCAPS > 0)
    if (pid() == 0) {
      coord center, max, min;
      foreach_dimension() {
        center.x = 0.;
        max.x = -HUGE;
        min.x = HUGE;
      }
      for(int i=0; i<CAPS(0).nln; i++) {
        foreach_dimension() {
          center.x += CAPS(0).nodes[i].pos.x;
          if (CAPS(0).nodes[i].pos.x > max.x) max.x = CAPS(0).nodes[i].pos.x;
          if (CAPS(0).nodes[i].pos.x < min.x) min.x = CAPS(0).nodes[i].pos.x;
        }
      }
      foreach_dimension() center.x /= CAPS(0).nln;
      fprintf(foutput, "%d %g ", i, t);
      foreach_dimension() fprintf(foutput, "%g %g %g ", center.x, min.x, max.x);
      fprintf(foutput, "\n");
      fflush(foutput);
    }
  #endif
}

#if NCAPS > 0
  event save (t += 10*OUTPUT_FREQ) {
    dump(file = "flow.dump");
    dump_capsules("caps.dump");
  }
#endif

/** We also output a snapshot every OUTPUT_FREQ iteration. These four snapshots
are then glued side-by-side together and made into a movie using ffmpeg. */
event pictures (i += OUTPUT_FREQ) {
  char name[32];
  view(fov = .33*19.09, bg = {1,1,1}, height=600, width=1800, camera = "right");
  clear();
  translate (z = -6.) {
    squares("u.x", n = {1,0,0}, map = cool_warm, min = -2*U_AVG, max = 2*U_AVG);
    cells(n = {1, 0, 0});
    draw_lags(lw = 1., edges = true, facets = true);
  }
  sprintf(name, "ux_cells_%d.png", i/OUTPUT_FREQ);
  save(name);
}

#if INIT_GEOM
  event vof (i = 0) {
    dump (file = "geometry_generation.dump");
    exit(0);
  }
#endif

event end (END_COND) {
  fclose(foutput);
  #if INITIAL_FLOW_FIELD
    view(fov = .33*19.09, bg = {1,1,1}, height=600, width=1800, camera = "right");
    clear();
    translate (z = -6.) {
      squares("u.x", n = {1,0,0}, map = cool_warm, min = -2*U_AVG, max = 2*U_AVG);
      cells(n = {1, 0, 0});
    }
    save("ux_initial_field.png");
    dump (file = "initial_flow_field.dump");
  #endif
  return 0.;
}

/**
## Results

![Three-dimensional trajectory of the capsule](https://damienhuet.github.io/images/basilisk_sandbox_files/helix/helix_3d_trajectory.png)

![Distance of the capsule from the helix centerline](https://damienhuet.github.io/images/basilisk_sandbox_files/helix/helix_radial_migration.png)

![Trajectory of the capsule in a cross-section of the helical pipe](https://damienhuet.github.io/images/basilisk_sandbox_files/helix/helix_plane_migration.png)

*/
