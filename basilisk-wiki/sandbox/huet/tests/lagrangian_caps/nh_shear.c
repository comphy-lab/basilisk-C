/**
# Deformation of an initially spherical capsule in a shear flow

![Initially spherical capsule in a shear flow at $Ca=0.2$](https://damienhuet.github.io/images/basilisk_sandbox_files/nh_shear/ux_zoom.mp4)(width="30%")

We reproduce a case from [Pozrikidis, 1995](#pozrikidis1995finite); [Ramanujan
& Pozrikidis, 1998](#ramanujan1998deformation) and  [Doddi and Bagchi, 2008](#doddi2008lateral).

In this case, an initially spherical capsule is deformed under the action of a
Stokes shear flow. The elastic law of the capsule is the neo-Hookean law.
To validate our results, we compute the Taylor deformation parameter
$D = \frac{r_{max} - r_{min}}{r_{max} + r_{min}}$ in the shear plane. The videos above correspond to Capillary numbers of 0.0125 (left) and 0.2 (right).


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

#define L0 8.
#define RADIUS 1.
#define SHEAR_RATE 1.
#ifndef RE
  #define RE 0.01
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(RADIUS)))
#ifndef CA
  #define CA 0.2
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
  #define TEND (8./SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 4
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 7
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
  #define OUTPUT_FREQ 10
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
  N = 1 << MAXLEVEL;
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
  /** ... and we initialize the flow field to that of an undisturbed,
  fully-developed shear.*/
  foreach() {
    u.x[] = SHEAR_RATE*y;
    u.y[] = 0.;
    u.z[] = 0.;
  }
  /** We initialize a spherical membrane using the pre-defined function in
  [common-shapes-ft.h](../../src/lagrangian_caps/common-shapes-ft.h) */
  activate_spherical_capsule(&CAPS(0), level = LAG_LEVEL, radius = RADIUS);
  foutput  = fopen("output.txt","w");
}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm. We also force the top and bottom walls to be at a fine level.
* second, we adapt according to the velocity field.
*/
event adapt (i++) {
  tag_ibm_stencils(&CAPS(0));
  foreach_boundary(top) stencils[] = noise();
  foreach_boundary(bottom) stencils[] = noise();
  adapt_wavelet({stencils, u}, (double []){1.e-2, U_TOL, U_TOL, U_TOL},
    minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  generate_lag_stencils(&CAPS(0));
}

/**
In the event below, we output the Taylor deformation parameter, as well as the angle between the largest radius $\bm{r_{max}}$
and the x-axis.
*/
event logfile (i++) {
  if (pid() == 0) {
    double rmax = -HUGE;
    double rmin = HUGE;
    double theta = 0.;
    for(int i=0; i<CAPS(0).nln; i++) {
      double x, y, z;
      x = CAPS(0).nodes[i].pos.x;
      y = CAPS(0).nodes[i].pos.y;
      z = CAPS(0).nodes[i].pos.z;
      double rad  = sqrt(sq(x) + sq(y) + sq(z));
      if (rad > rmax) {
        rmax = rad;
        theta = (fabs(x) < 1.e-14) ? (y>0. ? pi/2. : 3*pi/2.) : atan2(y,x);
      }
      if (rad < rmin)
        rmin = rad;
    }
    double D = (rmax - rmin)/(rmax + rmin);
    fprintf (foutput, "%g %g %g %g %g\n", t, rmin, rmax, D, theta);
    fflush(foutput);
  }
}

/** We also output a movie frame every OUTPUT_FREQ iteration */
event pictures (i+=OUTPUT_FREQ) {
  char fname[32];
  view(fov = 18.9, bg = {1,1,1}, camera = "front");
  clear();
  cells(n = {0,0,1});
  squares("u.x", n = {0,0,1}, map = cool_warm);
  draw_lag(&CAPS(0), lw = 1, edges = true, facets = true);
  sprintf(fname, "ux_%d.png", i/OUTPUT_FREQ);
  save(fname);
}

event end (t = TEND) {
  fclose(foutput);
  return 0.;
}

/**
## Results
### Qualitative
The command to generate the videos above from the simulation results:

```ffmpeg -y -framerate 100 -i ux_%d.png -c:v libx264 -vf "format=yuv420p, scale=1.1*iw:-1, crop=iw/1.1:ih/1.1" ux.mp4```

### Quantitative
~~~gnuplot
set grid
set xlabel "Non-dimensional time"
set ylabel "Taylor deformation parameter D"
set key reverse Left
set key top left
set xrange [0:5]
set yrange [0:.6]

set label "Ca 0.2" at 4.1,.525
set label "Ca 0.1" at 4.1,.41
set label "Ca 0.05" at 4.1,.3
set label "Ca 0.025" at 4.1,.19
set label "Ca 0.0125" at 4.1,0.105


plot '../data/nh_shear_3d/pozrikidis1998_ca_0.2.csv' every 4 w p pt 4 ps 1.25 lw 1.25 lc "royalblue" dt 3 title "", \
'../data/nh_shear_3d/pozrikidis1995_ca_0.2.csv' every 7 w p pt 8 ps 1.25 lw 1.5 lc "orange" dt 4 title "", \
'../data/nh_shear_3d/pozrikidis1998_ca_0.1.csv' every 4 w p pt 4 ps 1.25 lw 1.25 lc "royalblue" dt 3 title "", \
'../data/nh_shear_3d/pozrikidis1995_ca_0.1.csv' every 10 w p pt 8 ps 1.25 lw 1.25 lc "orange" dt 4 title "", \
'../data/nh_shear_3d/pozrikidis1998_ca_0.05.csv' every 4 w p pt 4 ps 1.25 lw 1.25 lc "royalblue" dt 3 title "", \
'../data/nh_shear_3d/pozrikidis1995_ca_0.05.csv' every 10 w p pt 8 ps 1.25 lw 1.25 lc "orange" dt 4 title "", \
'../data/nh_shear_3d/pozrikidis1998_ca_0.025.csv' every 3 w p pt 4 ps 1.25 lw 1.25 lc "royalblue" dt 3 title "", \
'../data/nh_shear_3d/pozrikidis1995_ca_0.025.csv' every 7 w p pt 8 ps 1.25 lw 1.25 lc "orange" dt 4 title "", \
'../data/nh_shear_3d/pozrikidis1995_ca_0.0125.csv' every 5 w p pt 8 ps 1.25 lw 1.25 lc "orange" dt 4 title "", \
'../data/nh_shear_3d/bagchi_ca0.2.csv' w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title "", \
'../data/nh_shear_3d/bagchi_ca0.1.csv' w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title "", \
'../data/nh_shear_3d/bagchi_ca0.05.csv' w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title "", \
'../data/nh_shear_3d/bagchi_ca0.025.csv' w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title"", \
'../data/nh_shear_3d/bagchi_ca0.0125.csv' w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title"", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4 every ::1::1 w l lw 1.25 lc "black" dt 1 title "This study", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4 every ::1::1 w l lw 1.25 lc "black" dt 2 title "This study - high resolution", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4 every ::1::1 w p pt 8 ps 1.25 lw 1.25 lc "orange" dt 4 title "Pozrikidis (1995)", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4 every ::1::1 w p pt 4 ps 1.25 lw 1.25 lc "royalblue" title "Ramanujan \\& Pozrikidis (1998)", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4 every ::1::1 w p pt 6 ps 1.25 lw 1.25 lc "dark-green" title "Doddi \\& Bagchi (2008)", \
'../data/nh_shear_3d/low_res_ca0.2.txt' using 1:4:(10.0) smooth acsplines w l lc -1 lw 2 title "", \
'../data/nh_shear_3d/low_res_ca0.1.txt' using 1:4:(10.0) smooth acsplines w l lc -1 lw 2 title "", \
'../data/nh_shear_3d/low_res_ca0.05.txt' using 1:4:(100.0) smooth acsplines w l lc -1 lw 2 title "", \
'../data/nh_shear_3d/low_res_ca0.025.txt' using 1:4:(100.0) smooth acsplines w l lc -1 lw 2 title "", \
'../data/nh_shear_3d/low_res_ca0.0125.txt' using 1:4:(100.0) smooth acsplines w l lc -1 lw 2 title "", \
'../data/nh_shear_3d/high_res_ca0.2.txt' using 1:4 w l lc -1 lw 2 dt 2 title "", \
'../data/nh_shear_3d/high_res_ca0.1.txt' using 1:4 w l lc -1 lw 2 dt 2 title "", \
'../data/nh_shear_3d/high_res_ca0.05.txt' using 1:4 w l lc -1 lw 2 dt 2 title "", \
'../data/nh_shear_3d/high_res_ca0.025.txt' using 1:4 w l lc -1 lw 2 dt 2 title "", \
'../data/nh_shear_3d/high_res_ca0.0125.txt' using 1:4 w l lc -1 lw 2 dt 2 title ""
~~~

## References
~~~bib
@Article{pozrikidis1995finite,
  author    = {Pozrikidis, C},
  journal   = {Journal of Fluid Mechanics},
  title     = {Finite deformation of liquid capsules enclosed by elastic membranes in simple shear flow},
  year      = {1995},
  pages     = {123--152},
  volume    = {297},
  file      = {:pozrikidis1995finite - Finite Deformation of Liquid Capsules Enclosed by Elastic Membranes in Simple Shear Flow.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}

@Article{ramanujan1998deformation,
  author    = {Ramanujan, S and Pozrikidis, C},
  journal   = {Journal of fluid mechanics},
  title     = {Deformation of liquid capsules enclosed by elastic membranes in simple shear flow: large deformations and the effect of fluid viscosities},
  year      = {1998},
  pages     = {117--143},
  volume    = {361},
  file      = {:deformation-of-liquid-capsules-enclosed-by-elastic-membranes-in-simple-shear-flow-large-deformations-and-the-effect-of-fluid-viscosities.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}

@Article{doddi2008lateral,
  author    = {Doddi, Sai K and Bagchi, Prosenjit},
  journal   = {International Journal of Multiphase Flow},
  title     = {Lateral migration of a capsule in a plane Poiseuille flow in a channel},
  year      = {2008},
  number    = {10},
  pages     = {966--986},
  volume    = {34},
  file      = {:bagchi2008.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Elsevier},
}
~~~
*/
