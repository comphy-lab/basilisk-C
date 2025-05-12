/**
# Deformation of an initially spherical capsule obeying the neo-Hookean elastic law and the Helfrich bending force in a shear flow

![Neo-Hookean capsule in a shear flow with bending stresses](https://damienhuet.github.io/images/basilisk_sandbox_files/bending_shear/ux.mp4)(width="30%")

## Definition of the relevant parameters

We define the following physical quantities:

* the size of the box $L_0 = 4$,
* the radius $a = 1$,
* the shear rate $\dot{\gamma} = 1$,
* the Reynolds number $Re \ll 1$,
* the viscosity $\mu = 1$,
* the density $\rho = \frac{Re \, \mu}{\dot{\gamma} a^2}$,
* the Capillary number defined $-$ in the case of capsules $-$ as the ratio of
viscous forces over elastic forces: $Ca = \frac{\mu a \dot{\gamma}}{E_S}$
* the elastic modulus $E_s$
* the bending modulus $E_b$, but it is more convenient to consider the non-dimensional bending modulus $\tilde{E}_b = \frac{E_b}{a^2 E_s}$
* the reference curvature, which is not considered in this case
*/

#define L0 4.
#define RADIUS 1.
#define SHEAR_RATE 1.
#ifndef RE
  #define RE 0.01
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(RADIUS)))
#ifndef CA
  #define CA 0.05
#endif
#define E_S (MU*RADIUS*SHEAR_RATE/CA)
#define ND_EB .0375
#define E_B (ND_EB*E_S*sq(RADIUS))
#define REF_CURV 0

/**
We also define some solver-ralated quantities:

* the non-dimensional duration time of the simulation TEND
* the minimum and maximum refinement levels of the Eulerian grid
* the number of refinement levels LAGLEVEL of the Lagrangian mesh of the capsule
* the maximum time-step DT_MAX, which appears to be limited by the bending modulus
* the tolerance of the Poisson solver (maximum admissible residual)
* the tolerance of the wavelet adaptivity algorithm for the velocity
* the output frequency: a picture is generated every OUTPUT_FREQ iterations
* the $stokes$ boolean, in order to ignore the convective term in the
Navier-Stokes solver
* the Jacobi preconditionner, which we turn on for this case.
*/

#ifndef TEND
  #define TEND (8./SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 5
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 7
#endif
#ifndef LAGLEVEL
  #define LAGLEVEL 5
#endif
#ifndef DT_MAX
  #define DT_MAX 2.5e-5
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-3
#endif
#ifndef U_TOL
  #define U_TOL 0.01
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 800
#endif
#ifndef STOKES
  #define STOKES true
#endif
#define JACOBI 1

/**
## Simulation setup

We import the octree grid, the centered Navier-Stokes solver, the Lagrangian
mesh, the neo-Hookean elasticity, the bending force in the front-tracking
context, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/bending-ft.h"
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
  activate_spherical_capsule(&CAPS(0), level = LAGLEVEL, radius = RADIUS);
  foutput  = fopen("output.txt","w");
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
  adapt_wavelet({stencils, u}, (double []){1.e-2, U_TOL, U_TOL, U_TOL},
    minlevel = MINLEVEL, maxlevel = MAXLEVEL);
  generate_lag_stencils(&CAPS(0));
}

/**
In the event below, we output the Taylor deformation parameter in the shear
plane $z = 0$, as well as the angle between the largest radius $\bm{r_{max}}$
and the x-axis.
*/
event logfile (i+=OUTPUT_FREQ) {
  if (pid() == 0) {
    double rmax = -HUGE;
    double rmin = HUGE;
    double max_edge_length = -HUGE;
    double theta = 0.;
    for(int i=0; i<CAPS(0).nle; i++) {
      double el = edge_length(&CAPS(0), i);
      if (el > max_edge_length) max_edge_length = el;
    }
    for(int i=0; i<CAPS(0).nln; i++) {
    /** The post-processing is only carried out if we are in the shear plane */
      if (sq(CAPS(0).nodes[i].pos.z) < .5*max_edge_length) {
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
     }
    double D = (rmax - rmin)/(rmax + rmin);
    fprintf (foutput, "%g %g %g %g %g\n", t, rmin, rmax, D, theta);
    fflush(foutput);
  }
}

/** We also output a movie frame every OUTPUT_FREQ iteration */
event pictures (i+=OUTPUT_FREQ) {
  view(fov = 18.9, bg = {1,1,1}, camera = "front");
  clear();
  cells(n = {0,0,1});
  squares("u.x", n = {0,0,1}, map = cool_warm);
  draw_lag(&CAPS(0), lw = 1, edges = true, facets = true);
  char name[32];
  sprintf(name, "ux_%d.png", i/OUTPUT_FREQ);
  save(name);
}

event end (t = TEND) {
  fclose(foutput);
  return 0.;
}

/**
### Qualitative
The command to generate the videos above from the simulation results:
```ffmpeg -y -r 25 -i ux_%d.png -c:v libx264 -vf "fps=25,format=yuv420p" ux.mp4```

### Quantitative
~~~gnuplot
set grid
set key bottom right
set xrange [0:5]
set ylabel "Taylor deformation parameter"
set xlabel "Non-dimensional time"

plot '../data/bending_shear/output.txt' using 1:4 w l lc "black" lw 2 dt 1 title "This study", \
'../data/bending_shear/pozrikidis_2001.csv' w l lc "web-blue" lw 1 dt 2 title "Pozrikidis 2001", \
'../data/bending_shear/methods_C_S.csv' w l lc "orange" lw 1 dt 3 title "Guckenberger et al., 2016", \
'../data/bending_shear/huang_2012.csv' w l lc "dark-red" lw 1 dt 4 title "Huang 2012", \
'../data/bending_shear/Le_2009.csv' w l lc "dark-spring-green" lw 1 dt 5 title "Le 2009", \
'../data/bending_shear/Le_2010.csv' w l lc "orchid4" lw 1 dt 6 title "Le 2010", \
'../data/bending_shear/zhu_2015.csv' w l lc "navy" lw 1 dt 7 title "Zhu 2015"
~~~

## References
~~~bib
@Article{pozrikidis2001effect,
  author    = {Pozrikidis, C},
  journal   = {Journal of Fluid Mechanics},
  title     = {Effect of membrane bending stiffness on the deformation of capsules in simple shear flow},
  year      = {2001},
  pages     = {269},
  volume    = {440},
  file      = {:pozrikidis2001effect - Effect of Membrane Bending Stiffness on the Deformation of Capsules in Simple Shear Flow.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}

@Article{le2010effect,
  author    = {Le, Duc Vinh},
  journal   = {Physical Review E},
  title     = {Effect of bending stiffness on the deformation of liquid capsules enclosed by thin shells in shear flow},
  year      = {2010},
  number    = {1},
  pages     = {016318},
  volume    = {82},
  file      = {:PhysRevE.82.016318.pdf:PDF},
  groups    = {Biological flows},
  publisher = {APS},
}

@Article{zhu2015motion,
  author    = {Zhu, Lailai and Brandt, Luca},
  journal   = {Journal of Fluid Mechanics},
  title     = {The motion of a deforming capsule through a corner},
  year      = {2015},
  pages     = {374--397},
  volume    = {770},
  file      = {:the-motion-of-a-deforming-capsule-through-a-corner.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}

@Article{guckenberger2016bending,
  author    = {Guckenberger, Achim and Schraml, Marcel P and Chen, Paul G and Leonetti, Marc and Gekle, Stephan},
  journal   = {Computer Physics Communications},
  title     = {On the bending algorithms for soft objects in flows},
  year      = {2016},
  pages     = {1--23},
  volume    = {207},
  file      = {:1-s2.0-S0010465516301072-main.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Elsevier},
}
~~~
*/
