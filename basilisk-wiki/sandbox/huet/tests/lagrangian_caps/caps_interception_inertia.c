/**
# Capsules interaction in the presence of inertia

![](https://damienhuet.github.io/images/basilisk_sandbox_files/caps_interception_inertia/uy.mp4)(width="80%")
The video above shows two elastic capsules interacting and reversing their motion at Reynolds 10. We borrowed this case from [Doddi \& Bagchi](#doddi2008effet)

## Definition of relevant parameters

*/

#define RADIUS 1.
#define DIAMETER (2*RADIUS)
#define L0 (DIAMETER/.16)
#define SHEAR_RATE 1.
#ifndef RE
  #define RE 10.
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(DIAMETER)))
#ifndef CA
  #define CA 0.05
#endif
#define E_S (MU*DIAMETER*SHEAR_RATE/CA)
#define NCAps 1.25
#define ALPHA_P 0.

#ifndef TEND
  #define TEND (40.*SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 5
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 8
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4
#endif
#ifndef DT_MAX
  #define DT_MAX (1.e-3)
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-3
#endif
#ifndef U_TOL
  #define U_TOL 0.05
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 100
#endif
#ifndef STOKES
  #define STOKES false
#endif
#define JACOBI 1

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

FILE* foutput = NULL;

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
  run();
}

/** We impose shear-flow boundary conditions. */
u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] =  dirichlet(SHEAR_RATE*y);
uf.n[bottom] = dirichlet(0);
uf.n[top] = dirichlet(0);

coord initial_centers[2] = {{-2.*DIAMETER, .1*DIAMETER, 0.},
  {2.*DIAMETER, -.1*DIAMETER, 0.}};
double initial_area;

event init (i = 0) {
  if (pid() == 0) foutput  = fopen("output.txt","w");

  /** We initialize two spherical membranes using the pre-defined functions in
  [common-shapes-ft.h](../../src/lagrangian_caps/common-shapes-ft.h) */
  activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL,
    radius = RADIUS, shift = {initial_centers[0].x,
    initial_centers[0].y, initial_centers[0].z});
  activate_spherical_capsule(&CAPS(1), level = LAG_LEVEL,
    radius = RADIUS, shift = {initial_centers[1].x,
    initial_centers[1].y, initial_centers[1].z});
}

event adapt_init (i=0) {
  /** For post-processing purposes, we compute the initial area of the
  membranes. */
  comp_triangle_area_normals(&CAPS(0));
  comp_triangle_area_normals(&CAPS(1));
  for(int j=0; j<CAPS(0).nlt; j++) initial_area += CAPS(0).triangles[j].area;

  /** We refine the mesh around the membranes.*/
  astats ss;
  int ic = 0;
  generate_lag_stencils();
  do {
    ic++;
    tag_ibm_stencils();
    ss = adapt_wavelet ({stencils}, (double[]) {1.e-30},
      minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    generate_lag_stencils();
  } while ((ss.nf || ss.nc) && ic < 100);

  /** ... and we initialize the flow field to that of an undisturbed,
  fully-developed shear.*/
  foreach() {
    u.x[] = SHEAR_RATE*y;
    u.y[] = 0.;
    u.z[] = 0.;
  }
}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm,
* second, we adapt according to the velocity field.
*/
event adapt (i++) {
  tag_ibm_stencils();
  adapt_wavelet({stencils, u}, (double []){1.e-2, U_TOL, U_TOL, U_TOL},
    minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  generate_lag_stencils();
}

/**
In the event below, we output the coordinates of the centroids of the two capsules as well as their areas and the total number of grid cells (to get an idea of the performance).
*/
event logfile (i += OUTPUT_FREQ) {
  if (pid() == 0) {
    coord current_centers[2];
    double current_areas[2];
    fprintf(foutput, "%g ", t);
    for(int k=0; k<NCAPS; k++) {
      current_areas[k] = 0.;
      foreach_dimension() current_centers[k].x = 0.;
      for(int l=0; l<CAPS(k).nln; l++)
        foreach_dimension()
          current_centers[k].x += CAPS(k).nodes[l].pos.x;
      foreach_dimension() current_centers[k].x /= CAPS(k).nln;
      for(int l=0; l<CAPS(k).nlt; l++)
        current_areas[k] += CAPS(k).triangles[l].area;
      fprintf(foutput, "%g %g %g %g ", current_centers[k].x,
        current_centers[k].y, current_centers[k].z, current_areas[k]);
    }
    fprintf(foutput, "%ld \n", grid->tn);
    fflush(foutput);
  }
}

/** We also output pictures of the x- and y-components of the velocity field and of the capsules at every OUTPUT_FREQ iterations */
event pictures (i += OUTPUT_FREQ) {
  char name[32];
  view(fov = .35*19, bg = {1,1,1}, camera = "front", height = 800, width = 2000);
  clear();
  cells(n = {0,0,1});
  squares("u.x", n = {0,0,1}, map = cool_warm, linear = true);
  draw_lags(lw = 1, edges = true, facets = true);
  sprintf(name, "ux_%d.png", i/OUTPUT_FREQ);
  save(name);

  view(fov = .35*19, bg = {1,1,1}, camera = "front", height = 800, width = 2000);
  clear();
  cells(n = {0,0,1});
  squares("u.y", n = {0,0,1}, map = cool_warm, linear = true);
  draw_lags(lw = 1, edges = true, facets = true);
  sprintf(name, "uy_%d.png", i/OUTPUT_FREQ);
  save(name);
}

event end (t = TEND) {
  if (pid() == 0) fclose(foutput);
  return 0.;
}

/**
## Results
~~~gnuplot
set ylabel "{/Symbol D}x_2/2a"
set xlabel "Non-dimensional time"
set xrange [0:30]
set yrange [-.25:.25]
set size ratio .5
set grid
set key top right

plot '../data/caps_interception_inertia/output_re3.txt' using 1:(($1 < 30) ? ($3)/2 : 1/0) w l lc "blue" lw 1.5 dt 1 title "Re = 3: Doddi \\& Bagchi ○, this study", \
'../data/caps_interception_inertia/output_re3.txt' using 1:(($1 < 30) ? ($7)/2 : 1/0) w l lc "blue" lw 1.5 dt 1 title "", \
'../data/caps_interception_inertia/output_re10.txt' using 1:(($1 < 30) ? ($3)/2 : 1/0) w l lc "black" lw 1.5 dt 2 title "Re = 10: Doddi \\& Bagchi □, this study", \
'../data/caps_interception_inertia/output_re10.txt' using 1:(($1 < 30) ? ($7)/2 : 1/0) w l lc "black" lw 1.5 dt 2 title "", \
'../data/caps_interception_inertia/output_re50.txt' using 1:(($1 < 30) ? ($3)/2 : 1/0) w l lc "dark-red" lw 1.5 dt 3 title "Re = 50: Doddi \\& Bagchi △, this study", \
'../data/caps_interception_inertia/output_re50.txt' using 1:(($1 < 30) ? ($7)/2 : 1/0) w l lc "dark-red" lw 1.5 dt 3 title "", \
'../data/caps_interception_inertia/doddi_fig8c_re3.csv' w p pt 6 ps 1.25 lc "blue" title "", \
'../data/caps_interception_inertia/doddi_fig8c_re3.csv' using 1:(-$2) w p pt 6 ps 1.25 lc "blue" title "", \
'../data/caps_interception_inertia/doddi_fig8c_re10.csv' w p pt 4 ps 1.25 lc -1 title "", \
'../data/caps_interception_inertia/doddi_fig8c_re10.csv' using 1:(-$2) w p pt 4 ps 1.25 lc -1 title "", \
'../data/caps_interception_inertia/doddi_fig8c_re50.csv' w p pt 8 ps 1.25 lc "dark-red" title "", \
'../data/caps_interception_inertia/doddi_fig8c_re50.csv' using 1:(-$2) w p pt 8 ps 1.25 lc "dark-red" title "", \
'../data/caps_interception_inertia/output_re3.txt' every 5000 w p ps 0 title ""

~~~

## References
~~~bib
@Article{doddi2008effect,
  author    = {Doddi, Sai K and Bagchi, Prosenjit},
  journal   = {International journal of multiphase flow},
  title     = {Effect of inertia on the hydrodynamic interaction between two liquid capsules in simple shear flow},
  year      = {2008},
  number    = {4},
  pages     = {375--392},
  volume    = {34},
  file      = {:1-s2.0-S0301932207001589-main.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Elsevier},
}
~~~
*/
