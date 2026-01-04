/**
# Wake forming behind a porous cylinder
This is a simulation of the flow past a porous cylinder testing the 
implementation of the Darcy-Forchheimer model in Basilisk.
We consider a 2D domain with a cylinder of radius R0 placed at the symmetry 
axis. The flow is driven by a constant inflow velocity U0 at the left 
boundary. We set the fluid viscosity muG to achieve a specified Reynolds
number. This test case is inspired by the work of [Yu et al. (2011)](https://doi.org/10.1016/j.compfluid.2010.09.040).

We run several cases varying the permeability of the porous medium, 
represented by the Darcy number Da. The wake length is measured as the point at
which *u.x* along the centerline (y=0) becomes zero again after the cylinder.

![Streamlines behind a porous cylinder for Da=0.001](porous-cylinder/streamlines.png){width="100%"}

## Simulation parameters
*/
int maxlevel = 9;         // Maximum refinement level
double Re = 20;           // Reynolds number
double R0 = 0.5;          // Cylinder radius 
double U0 = 1.;           // Inflow velocity
double epsi0 = 0.7;       // Porosity
double side_length = 15.; // Domain length in terms of R0
double tend = 10.;        // End time

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "darcy.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"

/**
## Boundary conditions
Boundary conditions are inlet on the left, outlet on the right, symmetry
on top and on bottom.
*/
u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);
pf[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);

scalar porosity[];
double rhoG = 1., muG;

/**
## Simulation setup
We declare the list of Darcy numbers to simulate.
*/
double DaList[] = {1.00E-05, 5.00E-05, 1.00E-04, 5.00E-04, 1.00E-03, 2.50E-03, 5.00E-03};
unsigned int ii = 0; // Index of the current case

int main() {

  /**
  We set the fluid properties based on the Reynolds number.
  For both phases, we want to use the gas properties.
  */
  muG = R0*2*U0/Re;
  rho1 = rho2 = rhoG;
  mu1 = mu2 = muG;

  L0 = side_length*R0*2;
  origin (-L0/2, 0);
  init_grid (1 << maxlevel);
  f.tracers = {porosity};

  /**
  We download the reference data from Yu et al. (2011) for comparison.
  */
  system ("wget -q https://raw.githubusercontent.com/Riccaraccio/basilisk-sandbox-rcaraccio/refs/heads/master/data/porouscylinder/yuData");

  /**
  We run the cases for each Darcy number in DaList.
  */
  size_t n_cases = sizeof(DaList)/sizeof(DaList[0]);
  for (ii = 0; ii < n_cases; ii++) {
    Da = (coord) {DaList[ii], DaList[ii]};
    run();
  }
}

/**
# Initialization event
We initialize a cylinder of radius R0 at the center of the domain.
Exploiting symmetry, we only simulate the upper half of the domain.
*/

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
event init (i = 0) {
  fraction (f, circle(x, y, R0));
  foreach() {
    porosity[] = f[]*epsi0;
    u.x[] = U0;
  }
}

/**
# Stability and Adaptation events
To ensure numerical stability, we limit the CFL number to a maximum value.
*/
const double max_cfl = 0.1;
event stability (i++) {
  if (CFL > max_cfl)
    CFL = max_cfl;
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({u.x, u.y}, {f}, (double[]){1.e-3, 1.e-3}, maxlevel, 2, padding=2);
}

/**
We avoid transport of the interface by setting the velocity to zero in the 
VOF event. After the interface advection, we restore the original velocity
field.
 */
face vector ufsave[];
event vof (i++) {
  foreach_face() {
    ufsave.x[] = uf.x[];
    uf.x[] = 0.;
  }
}

event tracer_diffusion (i++) {
  foreach_face()
    uf.x[] = ufsave.x[];
}

/**
## Embed implementation (commented out)
Here we use the embed method to represent the cylinder as a solid object.
This is to obtain the wake length for comparison with the porous cylinder cases.
*/
/*
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "adapt_wavelet_leave_interface.h"

u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);

face vector muv[];

int main() {
  mu = muv;

  L0 = side_length*R0*2;
  origin(-L0/2, 0);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event properties (i++){
  foreach_face()
    muv.x[] = fm.x[]*R0*2*U0/Re;
}

event init (i = 0) {
  mask (y > L0/2 ? top : none);
  solid (cs, fs, -circle(x, y, R0));
  foreach()
    u.x[] = cs[] > 1.- 1.e-10 ? U0 : 0;
}

event adapt (i++) {
  adapt_wavelet_leave_interface ({u.x, u.y}, {cs}, (double[]){1.e-3, 1.e-3}, maxlevel, 3, padding=1);
}
*/

/**
## Movie event
We visualize the streamlines developing around the cylinder for case 4.
*/

scalar strline[], omega[];
strline[top] = dirichlet(0);
strline[bottom] = dirichlet(U0 * L0);
event movie (t = end) {
  if (ii == 4) {
    foreach ()
      omega[] = 0;

    vorticity(u, omega);

    poisson(strline, omega);
    boundary({strline});

    view(quat = {0.000, 0.000, 0.000, 1.000},
         fov = 30, near = 0.01, far = 1000,
         tx = -0.03, ty = -0.05, tz = -0.3,
         width = 1920, height = 1080);
    draw_vof(c = "f");
    isoline(phi = "strline", n = 100, min = U0*L0*0.99, max = U0*L0*1.01);
    save("streamlines.png");
  }
} 

/**
At the end of each simulation, we output the centerline velocity profile
to a file named "case-ii.dat", where "ii" is the index of the current case.
*/
event stop (t = tend) {
  char name[80];
  sprintf(name, "case-%d.dat", ii);
  FILE * fp = fopen(name, "w");
  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0/2; x+=step){
    fprintf(fp, "%g %g\n", x, interpolate(u.x, x, 0));
  }
  fflush(fp);
  fclose(fp);
}

/**
## Post-processing script
~~~gnuplot wake length
reset
set terminal svg size 450,400
set output "wake-length-comparison.svg"

array Da[7] = [1.00E-05, 5.00E-05, 1.00E-04, 5.00E-04, 1.00E-03, 2.50E-03, 5.00E-03]
array x_zero[7]
X0 = 0.5 # Cylinder radius
do for [i=1:7] {
# Filename for each case (note: using i-1 to get case-0.dat through case-6.dat)
  filename = sprintf("case-%d.dat", i-1)

  # Find exact zero crossing using linear interpolation
  stats filename using 1:2 nooutput

  # Initialize variables for finding zero crossing
  x_prev = 0
  y_prev = 0
  x_at_zero = 0

  # Read through the file to find sign change
  do for [j=0:STATS_records-1] {
    stats filename using (column(1)):(column(2)) every ::j::j nooutput
    x_curr = STATS_min_x
    y_curr = STATS_min_y

    # Check if y changes sign between previous and current point
    if (j > 0) {
      if ((y_prev * y_curr < 0) && (x_curr > 0.55)) {
        # Linear interpolation to find x where y=0
        x_at_zero = x_prev - y_prev * (x_curr - x_prev) / (y_curr - y_prev)
        break
      }
      # Also check if y is exactly 0
      if (y_curr == 0) {
        x_at_zero = x_curr
        break
      }
    }

    x_prev = x_curr
    y_prev = y_curr
  }

  # Store the result
  x_zero[i] = x_at_zero - X0

  if (x_zero[i] < 0) {
    x_zero[i] = 0
  }

  # Print the result
  print sprintf("Case %d (Da = %.2e): x at y=0 is %.6f", i-1, Da[i], x_at_zero)
}

set xlabel "Da"
set ylabel "Relative wake length"
set format x "10^{%T}"
set logscale x
set key bottom left
set xrange [5e-6:1e-2]
set yrange [0:1.]
set grid

Wake_Embed = 0.923063153 # Obtained from a separate simulation using the embed method

plot x_zero u  (Da[$1]):(x_zero[$1]) w lp notitle pt 4 lc "blue" lw 2,\
      "yuData" u 1:2 w p pt 6 lc "red" t "Yu et al. (2011)", \
      Wake_Embed w l lc "black" t "Embed"

~~~

## References

~~~bib
@article{yu2011steady,
  title={Steady flow around and through a permeable circular cylinder},
  author={Yu, Peng and Zeng, Yan and Lee, Thong See and Chen, Xiao Bing and Low, Hong Tong},
  journal={Computers \& Fluids},
  volume={42},
  number={1},
  pages={1--12},
  year={2011},
  publisher={Elsevier}
}
~~~
*/
