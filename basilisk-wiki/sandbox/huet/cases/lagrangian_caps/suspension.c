/**
# Rheology of a suspension of deformable capsules in a shear flow

![Suspension of 50 capsules in a bi-periodic shear flow with volume fractions from 5\% to 30\%](https://damienhuet.github.io/images/basilisk_sandbox_files/suspension/suspension.png)(width="60%")

![](https://damienhuet.github.io/images/basilisk_sandbox_files/suspension/capsules_compressed.mp4)(width="20%")

## Definition of the relevant parameters

We define the following physical quantities:

* the size of the box $L_0 = 8$,
* the radius $a = 1$,
* the shear rate $\dot{\gamma} = 1$,
* the Reynolds number $Re \ll 1$,
* the viscosity $\mu = 1$,
* the density $\rho = \frac{Re \, \mu}{\dot{\gamma} a^2}$,
* the Capillary number defined $-$ in the case of capsules $-$ as the ratio of
viscous forces over elastic forces: $Ca = \frac{\mu a \dot{\gamma}}{E_S}$
* the elastic modulus $E_S$
*/

#define NCAPS 50
#define VOLUME_FRACTION 0.05
#define L0 1.
#define RADIUS (pow(NCAPS*3.1415926535*4/(3*VOLUME_FRACTION), -1./3)*L0)
#define SHEAR_RATE 1.
#ifndef RE
  #define RE 10.
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(RADIUS)))
#ifndef CA
  #define CA 0.1
#endif
#define E_S (MU*RADIUS*SHEAR_RATE/CA)

/**
We also define some solver-ralated quantities:

* the non-dimensional duration time of the simulation TEND
* the minimum and maximum refinement levels of the Eulerian grid
* the number of refinement levels LAGLEVEL of the Lagrangian mesh of the capsule
* the maximum time-step DT_MAX, which appears to be dependant on the Capillary number
* the tolerance of the Poisson solver (maximum admissible residual)
* the tolerance of the wavelet adaptivity algorithm for the velocity
* an frequency of post-processing output (here, every 10 iterations)
* the $stokes$ boolean, in order to ignore the convective term in the
Navier-Stokes solver
* the Jacobi preconditionner, which we switch on for this case.
*/

#ifndef TEND
  #define TEND (100./SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 4
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 8
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4
#endif
#ifndef DT_MAX
  #define DT_MAX 1.e-3
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-6
#endif
#ifndef U_TOL
  #define U_TOL 0.01
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 25
#endif
#ifndef OUTPUT_FREQ_PV_DUMP
  #define OUTPUT_FREQ_PV_DUMP 1000
#endif
#ifndef STOKES
  #define STOKES false
#endif
#define JACOBI 1
#define PARAVIEW_CAPSULE 1
#define PARAVIEW_FLOW_FIELD 0
#define OUTPUT_CAPS_NODE_TRI_INFO 0

#define MIN_DELTA (L0/(1 << MAXLEVEL))
#define MIN_INITIAL_GAP (3*MIN_DELTA)
#define MAX_POSITIONING_ATTEMPTS 1.e+6

#define rand_pos_periodic() (((double)rand() / (double)((unsigned)RAND_MAX + 1)\
  - 0.5)*L0)
#define rand_pos() (((double)rand() / (double)((unsigned)RAND_MAX + 1)\
  - 0.5)*(L0 - 2*(RADIUS + MIN_INITIAL_GAP)))

/**
## Simulation setup

We import the octree grid, the centered Navier-Stokes solver, the Lagrangian
mesh, the neo-Hookean elasticity, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"
#if PARAVIEW_FLOW_FIELD
  #include "lagrangian_caps/vtu_output.h"
#endif

bool no_intersection(lagMesh* caps, int k) {
  bool no_intersection = true;
  for(int i=0; i<k; i++) {
    foreach_dimension() {
      if (fabs(GENERAL_SQNORM(caps[i].centroid, caps[k].centroid)) < sq(2*RADIUS
        + MIN_INITIAL_GAP)) {
        no_intersection = false;
        break;
      }
    }
  }
  return no_intersection;
}

FILE* foutput = NULL;
FILE* fperf = NULL;

const scalar myrho[] = RHO;
const face vector mymu[] = {MU, MU, MU};
const face vector myalpha[] = {1./RHO, 1./RHO, 1./RHO};

int main(int argc, char* argv[]) {
  origin(-0.5*L0, -0.5*L0, -0.5*L0);
  /** We set periodic boundary conditions on the non-horizontal walls.*/
  periodic(right);
  periodic(front);
  N = 1 << MINLEVEL;
  mu = mymu;
  rho = myrho;
  alpha = myalpha;
  TOLERANCE = MY_TOLERANCE;
  /** We don't need to compute the convective term in this case, so we set the
  boolean $stokes$ to false. However it is still important to choose $Re \ll 1$
  since we are solving the unsteady Stokes equation. */
  stokes = STOKES;
  DT = DT_MAX;
  srand ( time(NULL) );
  run();
}

/** We impose shear-flow boundary conditions... */
u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] =  dirichlet(SHEAR_RATE*y);
uf.n[bottom] = dirichlet(0);
uf.n[top] = dirichlet(0);

event init (i = 0) {
  foutput  = fopen("output.txt","w");
  fperf = fopen("perf.csv", "w");
  fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");
  /** ... and we initialize the flow field to that of an undisturbed,
  fully-developed shear.*/
  foreach() {
    u.x[] = SHEAR_RATE*y;
    u.y[] = 0.;
    u.z[] = 0.;
  }

  for(int k=0; k<NCAPS; k++)
    activate_spherical_capsule(&CAPS(k), radius = RADIUS, level = LAG_LEVEL);

  if (pid() == 0) {
    for(int k=0; k<NCAPS; k++) {
      bool keep_drawing_positions = true;
      int nb_attempts = 0;
      while (keep_drawing_positions && nb_attempts < MAX_POSITIONING_ATTEMPTS) {
        CAPS(k).centroid.x = rand_pos_periodic();
        CAPS(k).centroid.y = rand_pos();
        CAPS(k).centroid.z = rand_pos_periodic();
        if (no_intersection(allCaps.caps, k)) {
          for(int i=0; i<CAPS(k).nln; i++)
            foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
          keep_drawing_positions = false;
          fprintf(stderr, "Number of attempts to insert capsule %d: %d\n", k,
            nb_attempts+1);
        }
        nb_attempts++;
      }
      if (nb_attempts == MAX_POSITIONING_ATTEMPTS) {
        fprintf(stderr, "Error: max number of attempts to insert capsule \
          %d reached.\n", k);
        return 1;
      }
    }
  }
  #if _MPI
    /** We now inform all the other processors of the positions of the capsules */
    double centroids[3*NCAPS];
    if (pid() == 0) {
      for(int k=0; k<NCAPS; k++) {
        centroids[3*k] = CAPS(k).centroid.x;
        centroids[3*k+1] = CAPS(k).centroid.y;
        centroids[3*k+2] = CAPS(k).centroid.z;
      }
    }
    MPI_Bcast(centroids, 3*NCAPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (pid() > 0) {
      for(int k=0; k<NCAPS; k++) {
        CAPS(k).centroid.x = centroids[3*k];
        CAPS(k).centroid.y = centroids[3*k+1];
        CAPS(k).centroid.z = centroids[3*k+2];
        for(int i=0; i<CAPS(k).nln; i++)
          foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
      }
    }
  #endif
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
  generate_lag_stencils(no_warning = true);
  astats ss;
  int ic = 0;
  do {
    ic++;
    tag_ibm_stencils();
    ss = adapt_wavelet ({stencils}, (double[]) {1.e-30}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
    generate_lag_stencils(no_warning = true);
    fprintf(stderr, "Refine initial mesh: step %d\n", ic);
  } while ((ss.nf || ss.nc) && ic < 100);
}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm. We also force the top and bottom walls to be at a fine level.
* second, we adapt according to the velocity field.
*/
event adapt (i++) {
  tag_ibm_stencils();
  adapt_wavelet({stencils, u}, (double []){1.e-2, U_TOL, U_TOL, U_TOL},
    minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  generate_lag_stencils();
}

int nb_dump = 0;
event output_and_dump_caps (i+=OUTPUT_FREQ_PV_DUMP) {
  #if OUTPUT_CAPS_NODE_TRI_INFO
    for(int k=0; k<NCAPS; k++) {
      char fposname[64];
      char ftriname[64];
      sprintf(fposname, "mb_%d_pos.csv", k);
      sprintf(ftriname, "mb_%d_tri.csv", k);
      dump_plain_nodes_pos(&CAPS(k), fposname);
      dump_plain_triangles(&CAPS(k), ftriname);
    }
  #endif
  #if PARAVIEW_CAPSULE
    pv_output_ascii();
  #endif
    #if PARAVIEW_FLOW_FIELD
      scalar * list = {p};
      vector * vlist = {u};
      save_data(list, vlist, t);
    #endif

    char fname[62];
    sprintf(fname, "flow_%d.dump", nb_dump);
    dump(file = fname);
    sprintf(fname, "caps_%d.dump", nb_dump);
    dump_capsules(fname);
    nb_dump++;
}

event logfile (i+=OUTPUT_FREQ) {
  double top_visc_stress = 0;
  int top_nb_cells = 0;
  foreach_boundary(top, reduction(+:top_visc_stress) reduction(+:top_nb_cells)) {
    top_nb_cells++;
    top_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  top_visc_stress *= MU/(sq(L0));

  double bottom_visc_stress = 0;
  int bottom_nb_cells = 0;
  foreach_boundary(bottom, reduction(+:bottom_visc_stress), reduction(+:bottom_nb_cells)) {
    bottom_nb_cells++;
    bottom_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  bottom_visc_stress *= MU/(sq(L0));

  if (pid() == 0) {
    double avg_ncaps_area = 0;
    double avg_ncaps_volume = 0;
    for(int k=0; k<NCAPS; k++) {
      avg_ncaps_volume += CAPS(k).volume;
      for(int i=0; i<CAPS(k).nlt; i++) avg_ncaps_area += CAPS(k).triangles[i].area;
    }
    avg_ncaps_area /= NCAPS*4*pi*sq(RADIUS);
    avg_ncaps_volume /= NCAPS*CAPS(0).initial_volume;

    fprintf(foutput, "%d %g %d %g %d %g %g %g\n", i, t, top_nb_cells,
      top_visc_stress, bottom_nb_cells, bottom_visc_stress, avg_ncaps_area,
      avg_ncaps_volume);
    fflush(foutput);

    fprintf(fperf, "%ld %d %g %g %d %d %g %g %d\n", grid->tn, mgu.i, mgu.resb,
      mgu.resa, mgu.nrelax, mgp.i, mgp.resb, mgp.resa, mgp.nrelax);
    fflush(fperf);
  }
}

/** We also output a movie frame every OUTPUT_FREQ iteration */
int nb_pic = 0;
event pictures (i+=OUTPUT_FREQ) {
  char fname[32];
  view(fov = 25, bg = {1,1,1}, camera = "front");
  clear();
  // cells(n = {0,0,1}, alpha = -.49*L0);
  cells(n = {0,0,1});
  // squares("u.x", n = {0,0,1}, alpha = -.49*L0, map = cool_warm);
  draw_lags(lw = .5, edges = true, facets = true);
  sprintf(fname, "ux_%d.png", nb_pic);
  save(fname);
  nb_pic++;
}

event end (t = TEND) {
  fclose(foutput);
  return 0.;
}

/**
## Results

~~~gnuplot
set multiplot

set xlabel "Non-dimensional time"
set ylabel "Normalized bulk viscosity μ_{bulk}"
set yrange [1:3.1]
set label 1 "Φ = 5%" at 75,1.35
set label 2 "Φ = 10%" at 75,1.62
set label 3 "Φ = 15%" at 75,1.8
set label 4 "Φ = 20%" at 75,2
set label 5 "Φ = 30%" at 75,2.42

stats 'phi0.05.txt' every ::1200 using (.5*($4 + $6)) prefix "A" nooutput
stats 'phi0.1.txt' every ::1200 using (.5*($4 + $6)) prefix "B" nooutput
stats 'phi0.15.txt' every ::1200 using (.5*($4 + $6)) prefix "C" nooutput
stats 'phi0.2.txt' every ::1200 using (.5*($4 + $6)) prefix "D" nooutput
stats 'phi0.3.txt' every ::1200 using (.5*($4 + $6)) prefix "E" nooutput


plot "phi0.05.txt" using 2:(.5*($4 + $6)) w l lc -1 dt 1 title "", \
"phi0.1.txt" using 2:(.5*($4 + $6)) w l lc -1 dt 1 title "", \
"phi0.15.txt" using 2:(.5*($4 + $6)) w l lc -1 dt 1 title "", \
"phi0.2.txt" using 2:(.5*($4 + $6)) w l lc -1 dt 1 title "", \
"phi0.3.txt" using 2:(.5*($4 + $6)) w l lc -1 dt 1 title "", \

set size .4, .3
set origin .07, .68
array phi[5]
array viscosity[5]
phi[1]=5
phi[2]=10
phi[3]=15
phi[4]=20
phi[5]=30
viscosity[1]=A_mean
viscosity[2]=B_mean
viscosity[3]=C_mean
viscosity[4]=D_mean
viscosity[5]=E_mean
set xlabel "Φ (\%)" offset -1,1
set ylabel "⟨μ_{bulk}⟩" offset 2.8,.5
set grid
set xrange [0:35]
set yrange [1:2.5]
set obj 1 rect from graph 0,0 to graph 1,1 fc rgb 0xffffff behind
set xtics 0,10,30 offset 0,.25
set ytics 1,.5,2.5 offset 0.5
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5
plot sample [i=1:5] '+' using (phi[i]):(viscosity[i]) w p lc -1 pt 1 title ""

unset multiplot
~~~
*/