/**
# Capsules interception
We reprodue a case from [Lac and Barthès-Biesel](#lac2007hydrodynamic), where two elastic
capsules pass each other in a shear flow.

![Interception of two initially spherical capsules in a shear flow in Stokes conditions](https://damienhuet.github.io/images/basilisk_sandbox_files/caps_interception/uy_linear.mp4)(width="80%")

## Definition of relevant parameters

### Physical parameters
*/
#define L0 16.
#define RADIUS 1.
#define SHEAR_RATE 1.
#ifndef RE
  #define RE 0.01
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(RADIUS)))
#ifndef CA
  #define CA 0.45
#endif
#define E_S (MU*RADIUS*SHEAR_RATE/CA)
#define NCAPS 2
#define ALPHA_P 0.05
#ifndef TEND
  #define TEND (2.5e-4)
#endif

/**
### Numerical parameters
*/
#ifndef MINLEVEL
  #define MINLEVEL 5
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 9
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4
#endif
#ifndef DT_MAX
  #define DT_MAX (2.5e-4)
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-3
#endif
#ifndef U_TOL
  #define U_TOL 0.05
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 40
#endif
#ifndef STOKES
  #define STOKES true
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

/**
The initial position of the capsules are set according to figure 5 of Lu and
Barthès-Biesel. The capsules are slightly offset from the line $y=0$, meaning
their velocity is non-zero and of opposite sign: they will intercept.
*/
coord initial_centers[2] = {{-4.*RADIUS, .25*RADIUS, 0.},
  {4.*RADIUS, -.25*RADIUS, 0.}};
double initial_area;

event init (i = 0) {
  if (pid() == 0) foutput  = fopen("output.txt","w");

  /** We initialize two spherical membranes using the pre-defined function in
  [common-shapes-ft.h](../../src/lagrangian_caps/common-shapes-ft.h) */
  activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL,
    radius = RADIUS/(1. + ALPHA_P), shift = {initial_centers[0].x,
    initial_centers[0].y, initial_centers[0].z});
  activate_spherical_capsule(&CAPS(1), level = LAG_LEVEL,
    radius = RADIUS/(1. + ALPHA_P), shift = {initial_centers[1].x,
    initial_centers[1].y, initial_centers[1].z});
}

event adapt_init (i=0) {
  /**  The membrane are pre-inflated, so the stress-free configuration
  corresponds to a smaller radius $a_0$, and the current radius $a$ is related
  to $a_0$ by the prestress coefficient $\alpha_p$: $a = a_0 (1 + \alpha_p)$.*/
  for(int k=0; k<NCAPS; k++) {
    for(int l=0; l<CAPS(0).nln; l++) {
      double cr = 0.;
      foreach_dimension() cr += sq(CAPS(k).nodes[l].pos.x - initial_centers[k].x);
      cr = sqrt(cr);
      foreach_dimension()
        CAPS(k).nodes[l].pos.x = (CAPS(k).nodes[l].pos.x - initial_centers[k].x)*
          RADIUS/cr + initial_centers[k].x;
    }
  }

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
In the event below, we output the coordinates of the centroids of the two capsules, as well as their areas.
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
    fprintf(foutput, "\n");
    fflush(foutput);
  }
}

/** We also output a movie frame every OUTPUT_FREQ iteration */
event pictures (i += OUTPUT_FREQ) {
  char name[32];
  view(fov = .35*19, bg = {1,1,1}, camera = "front", height = 800, width = 2000);
  clear();
  cells(n = {0,0,1});
  squares("u.x", n = {0,0,1}, map = cool_warm);
  draw_lags(lw = 1, edges = true, facets = true);
  sprintf(name, "ux_%d.png", i);
  save(name);

  view(fov = .35*19, bg = {1,1,1}, camera = "front", height = 800, width = 2000);
  clear();
  cells(n = {0,0,1});
  squares("u.y", n = {0,0,1}, map = cool_warm);
  draw_lags(lw = 1, edges = true, facets = true);
  sprintf(name, "uy_%d.png", i);
  save(name);
}

event end (t = TEND) {
  if (pid() == 0) fclose(foutput);
  return 0.;
}

/**
## Results
~~~gnuplot
set xlabel "{/Symbol D}x_1/2a"
set ylabel "{/Symbol D}x_2/2a"
set xrange [-4:6]
set yrange [0:.875]
set size ratio .33
set grid

plot '../data/caps_interception/output.txt' using (($2-$6)/2): (($1 < 26) ? ($3-$7)/2 : 1/0) w l lc "blue" lw 1.5 title "This study", \
'../data/caps_interception/lac_fig4b.csv' using 1:2 w p pt 6 lc -1 ps 1.5 title "Lac et al. (2007)"

~~~

## References
~~~bib
@Article{lac2007hydrodynamic,
  author    = {Lac, Etienne and Morel, Arnaud and Barth{\`e}s-Biesel, Dominique},
  journal   = {Journal of Fluid Mechanics},
  title     = {Hydrodynamic interaction between two identical capsules in simple shear flow},
  year      = {2007},
  pages     = {149--169},
  volume    = {573},
  file      = {:hydrodynamic-interaction-between-two-identical-capsules-in-simple-shear-flow.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}
~~~


*/
