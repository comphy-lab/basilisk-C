/**
# 3D deformation of an initially spherical capsule in a constricted channel

We reproduce a case from [Park \& Dimitrakopoulos](#park2013transient) where an
elastic capsule subject to the Skalak elastic law flows in a constricted square
channel. In this case we consider the capsule diameter to be equal to the
constriction height, and the viscosity ratio of the inner and outer fluids is
unity.

![Initially spherical capsule flowing through a constricted channel](https://damienhuet.github.io/images/basilisk_sandbox_files/constricted_channel/channel.mp4)(width="80%")

## Definition of the relevant parameters

We define the following solver-ralated quantities:

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
  #define MAXLEVEL 9
#endif
#ifndef TEND
  #define TEND (6.)
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4
#endif
#ifndef DT_MAX
  #define DT_MAX 1.e-4
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-3
#endif
#ifndef U_TOL
  #define U_TOL 0.01
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 100
#endif
#ifndef STOKES
  #define STOKES true
#endif
#define JACOBI 1

/**
We also define the following physical quantities:

* the length of the channel $L_0 = 20$,
* the half-height $H_c = 1$ of the constriction of the channel
* the radius of the capsule $a = H_c = 1$,
* the average x-component of the velocity at the inlet and outlet $u_{\text{avg}} = 1$
* the Reynolds number $Re \ll 1$,
* the viscosity $\mu = 1$,
* the density $\rho = \frac{Re \, \mu}{u_{\text{avg}} 2 a}$,
* the Capillary number defined $-$ in the case of capsules $-$ as the ratio of
viscous forces over elastic forces: $Ca = \frac{\mu u_{\text{avg}}}{E_s}$
* the elastic modulus $E_s$ of the membrane
* the area dilatation modulus: $C=1$ is a large value to enforce approximate
area conservation
* the prestress coefficient of the membrane $\alpha_p = 0.05$
*/

#define L0 20.
#define HC (1. + 1.e-7)
#define RADIUS (HC)
#define U_AVG 1.
#ifndef RE
  #define RE 0.01
#endif
#define MU 1.
#define RHO (RE*MU/(U_AVG*2*RADIUS))
#ifndef CA
  #define CA 0.1
#endif
#define E_S (MU*U_AVG/CA)
#define AREA_DILATATION_MODULUS 1.
#define ALPHA_P 0.05

/**
## Simulation setup

We import the octree grid, the centered Navier-Stokes solver, the Lagrangian
mesh, the Skalak elastic law, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/skalak-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

FILE* foutput = NULL;
FILE* fperf = NULL;

scalar rhov[];
face vector muv[];
face vector alphav[];

int main(int argc, char* argv[]) {
  origin(-0.5*L0, -0.5*L0, -0.5*L0);
  N = 1 << MINLEVEL;
  init_grid(N);
  mu = muv;
  alpha = alphav;
  rho = rhov;
  TOLERANCE = MY_TOLERANCE;
  /** We don't need to compute the convective term in this case, so we set the
  boolean $stokes$ to false. However it is still important to choose $Re \ll 1$
  since we are solving the unsteady Stokes equation. */
  stokes = STOKES;
  DT = DT_MAX;
  run();
}

/** We impose the boundary conditions: the velocity at the
inlet and outlet is set to $u_{\text{avg}}$ in the x-direction and to $0$ otherwise. */
u.n[left] = dirichlet(U_AVG);
u.n[right] = dirichlet(U_AVG);
u.t[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
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

event init (i = 0) {
  foutput  = fopen("output.txt","w");
  fperf = fopen("perf.txt", "w");
  fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");

  /** We initialize a spherical membrane using the pre-defined function in
  [common-shapes-ft.h](../../src/lagrangian_caps/common-shapes-ft.h).*/
  activate_spherical_capsule(&CAPS(0), level=LAG_LEVEL, radius=RADIUS/(1. + ALPHA_P));

  /** In this case we use a third order face-flux interpolation. */
  for (scalar s in {u})
    s.third = true;
}

event init_adapt (i = 0) {
  /**  The membrane is pre-inflated, so the stress-free configuration
  corresponds to a smaller radius $a_0$, and the current radius $a$ is related
  to $a_0$ by the prestress coefficient $\alpha_p$: $a = a_0 (1 + \alpha_p)$.*/
  for(int k=0; k<CAPS(0).nln; k++) {
    double cr = 0.;
    foreach_dimension() cr += sq(CAPS(0).nodes[k].pos.x);
    cr = sqrt(cr);
    foreach_dimension() CAPS(0).nodes[k].pos.x *= RADIUS/cr;
    CAPS(0).nodes[k].pos.x -= 4*HC;
  }
  generate_lag_stencils(&CAPS(0));

  /** We refine the mesh around the embedded geometry, around the membrane, and
  at the inlet and outlet.*/
  astats ss;
  int ic = 0;
  do {
    ic++;
    tag_ibm_stencils(&CAPS(0));
    foreach_boundary(left)
      if (cs[] > 1.e-10) stencils[] = noise();
    foreach_boundary(right)
      if (cs[] > 1.e-10) stencils[] = noise();
    solid (cs, fs, intersection(intersection(-fabs(y) + 2*HC,
      -fabs(z) + 2*HC), (fabs(x) < HC) ? -fabs(y) + HC :
      -fabs(y) + 2*HC));
    ss = adapt_wavelet ({stencils, cs}, (double[]) {1.e-30, 1.e-30},
      maxlevel = MAXLEVEL, minlevel = MINLEVEL);
    generate_lag_stencils(&CAPS(0));
  } while ((ss.nf || ss.nc) && ic < 100);

  /** The initial condition for the x-component of the velocity field is set to
  the average velocity everywhere, except in the constriction where the velocity
  is doubled. The y and z components are zero.*/
  foreach() u.x[] = (fabs(x) < HC) ? 2*U_AVG*cm[] : U_AVG*cm[];
  foreach_face() uf.x[] = (fabs(x) < HC) ? 2*U_AVG*fm.x[] : U_AVG*fm.x[];
}


/** Due to embedded boundaries, the viscosity and density fields are not
constant since they have to be multiplied by the volume fraction. */
event properties (i++) {
  foreach_face() {
    muv.x[] = MU*fm.x[];
    alphav.x[] = 1./RHO*fm.x[];
  }
  boundary((scalar*){muv});
  foreach() rhov[] = RHO*cm[];
}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm,
* second, we adapt according to the velocity field.
*/
event adapt (i++) {
  tag_ibm_stencils(&CAPS(0));
  foreach_boundary(left)
    if (cs[] > 1.e-10) stencils[] = noise();
  foreach_boundary(right)
    if (cs[] > 1.e-10) stencils[] = noise();
  adapt_wavelet({cs, stencils, u}, (double []){1.e-30,
    1.e-30, U_TOL, U_TOL, U_TOL}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  generate_lag_stencils(&CAPS(0));
}

/**
In the event below, we output the maximum length of the capsule in the x, y and
z directions, as well as the position of the centroid of the capsule.
*/
event logfile (i += OUTPUT_FREQ) {
  if (pid() == 0) {
    double min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = HUGE; max_x = -HUGE;
    min_y = HUGE; max_y = -HUGE;
    min_z = HUGE; max_z = -HUGE;
    coord xc;
    foreach_dimension() xc.x = 0.;
    for(int i=0; i<CAPS(0).nln; i++) {
      coord* cni = &(CAPS(0).nodes[i].pos);
      if (cni->x > max_x) max_x = cni->x;
      if (cni->x < min_x) min_x = cni->x;
      if (cni->y > max_y) max_y = cni->y;
      if (cni->y < min_y) min_y = cni->y;
      if (cni->z > max_z) max_z = cni->z;
      if (cni->z < min_z) min_z = cni->z;
      foreach_dimension() xc.x += cni->x;
    }
    double l_x = max_x - min_x;
    double l_y = max_y - min_y;
    double l_z = max_z - min_z;
    foreach_dimension() xc.x /= CAPS(0).nln;
    fprintf(foutput, "%g %g %g %g %g %g %g\n", t, xc.x, xc.y, xc.z,
      l_x, l_y, l_z);
    fflush(foutput);

    fprintf(fperf, "%ld %d %g %g %d %d %g %g %d\n", grid->tn, mgu.i, mgu.resb, mgu.resa, mgu.nrelax, mgp.i, mgp.resb, mgp.resa, mgp.nrelax);
  }
}

/** We also output a snapshot every OUTPUT_FREQ iteration. These four snapshots
are then glued side-by-side together and made into a movie using ffmpeg. */
event pictures (i += OUTPUT_FREQ) {
  if (i == 0) {
    view(fov = 19., bg = {1,1,1}, tx = 0.);
    clear();
    squares("u.x", n = {0,0,1}, map = cool_warm,
      min = 0., max = 3*U_AVG);
    cells(n = {0, 0, 1});
    draw_lag(&CAPS(0), lw = 1., edges = true, facets = true);
    save("initial_field.png");
  }
  char name[32];
  for(int k=0; k<4; k++) {
    view(fov = .25*19.09  , bg = {1,1,1}, tx = (3./8. - k*.25));
    clear();
    squares("u.x", n = {0,0,1}, map = cool_warm,
      min = 0., max = 3*U_AVG);
    cells(n = {0, 0, 1});
    draw_lag(&CAPS(0), lw = 1., edges = true, facets = true);
    sprintf(name, "ux_%d_%d.png", i, k);
    save(name);
  }
}

event end (t = TEND) {
  fclose(foutput);
  return 0.;
}

/**
## Results

This case ran on 96 processors for 3 days. The result
below show very good qualitative and quantitative agreement with the results
of [Park \& Dimitrakopoulos](#park2013transient).

~~~gnuplot
set title "Deformation of a capsule of radius a flowing through a constriction\
 of half-size H_c = a"
set encoding utf8
set xlabel "x_c/H_c"
set ylabel "l_{x,y,z}/2a"
set grid
set xrange [-4:6]

plot "../data/constricted_channel/output.txt" using 2:($5/2) w l lc -1 dt 1 title "", \
"../data/constricted_channel/output.txt" using 2:($6/2) w l lc -1 dt 1 title "", \
"../data/constricted_channel/output.txt" using 2:($7/2) w l lc -1 dt 1 title "", \
"../data/constricted_channel/lx_park.csv" every 8 w p lc -1 pt 6 ps 1.5 title "", \
"../data/constricted_channel/ly_park_even_xsampling.csv" w p lc -1 pt 8 ps 1.5 title "", \
"../data/constricted_channel/lz_park.csv" every 2 w p lc -1 pt 4 ps 1.5 title "", \
"https://raw.githubusercontent.com/DamienHuet/basilisk-sandbox-files/main/output.txt" using 2:($5/2) every 1000000 w p ps 0 title "Park \\& Dimitrakopoulos: lx ○, ly □, lz △", \
"https://raw.githubusercontent.com/DamienHuet/basilisk-sandbox-files/main/output.txt" using 2:($5/2) every 1000000 w p ps 0 title "Basilisk ―――"
~~~

![Qualitative comparison of the capsule shape](https://damienhuet.github.io/images/basilisk_sandbox_files/constricted_channel/constricted_channel_qualitative_comparison.png){ width="30%" }

## References
~~~bib
@Article{park2013transient,
  author    = {Park, Sun-Young and Dimitrakopoulos, P},
  journal   = {Soft matter},
  title     = {Transient dynamics of an elastic capsule in a microfluidic constriction},
  year      = {2013},
  number    = {37},
  pages     = {8844--8855},
  volume    = {9},
  file      = {:park2013transient - Transient Dynamics of an Elastic Capsule in a Microfluidic Constriction.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Royal Society of Chemistry},
}
~~~
*/
