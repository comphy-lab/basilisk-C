/**
# Stefan Flow Problem
Section 4.1 of [Boyd et al, 2023](#boyd_consistent_2023).
*/

/**

![Stefan Flow](stefan_problem/temperature.mp4)(width="640" height="640")

*/
/**
## Setup the problem
### Domain size, initial vapor layer thickness, and simulation duration
*/
#define L 1.0e-3 
#define H0 0.1e-3
//
double start_time_stef = 0.0;
#define T_END (0.12-start_time_stef)

/**
### Fluid properties for saturatate water liquid and vapor
*/
#include "../src/properties/water_vapor.h"

/**
### Include header for NS, phase-change, gravity, etc
*/
#include "../src/01_vaporization/evap_include.h"

/**
### Boundary conditions for the high-temperature wall (bottom) and outflow (top)
*/

#define T_wall (T_sat + 10.0)
T_V[bottom] = dirichlet(T_wall);
fE_V[bottom] = dirichlet(rhocp_V * T_wall);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);

/**
### Analytical solution to the stefan problem
*/
#include "../src/theory_evap/stef_calc.h"

/**
### Main function.
*/
int main(int argc, char* argv[]) {
  MIN_LEVEL = 5;
  LEVEL = 5;
  MAX_LEVEL = 6;

  setup_evap();
  
  size(L);
  origin(0.0 * L, 0.0 * L);
  N = 1 << LEVEL;
  init_grid(N);

  G.x = 0.;
  Z.x = 0.;

  run();
}

/**
### Initialize the VOF, liquid and vapor temperature.
*/
#define plane(x, y, H) (H - y)
double beta_ = 0.;

event init(i = 0) {
  CFL = 0.2;

#if TREE
  refine(level < MAX_LEVEL && plane(x, y, (H0 - dR_refine)) < 0. &&
       plane(x, y, (H0 + dR_refine)) > 0.);
#endif
  fraction_LS(f, -plane(x, y, H0));

  beta_ = find_beta();
  start_time_stef = find_start_time(H0, beta_);

  foreach () {
    foreach_dimension() {
      u.x[] = 0.;
    }

    double profile = max(min(y / H0, 1.0), 0.0);
    double T_vap = profile * T_sat + (1.0 - profile) * T_wall;

    T_V[] = (1. - f[]) * T_vap + f[] * T_sat;
    T_L[] = (1. - f[]) * T_sat + f[] * T_sat;
  }

  init_Energy();
}


/**
### Adaptive mesh refinement.
*/
#define femax 1.0e-2
#define Temax 1.0
#define uemax 1.0e-1
#include "../src/01_vaporization/adapt_evap.h"

/**
### Outputs the interface location and the analytical interface location.
*/
event outputs(t = 0.; t += DELTA_T; t <= T_END) {

  double effective_height;
  scalar fc[];
  foreach ()
    fc[] = 1. - f[];

  effective_height = statsf(fc).sum / L0;

  char name[200];
  sprintf(name, "interface_location.dat");
  static FILE* fp = fopen(name, "w");

  double exact_location = interface_location_func(t+start_time_stef, beta_);
  fprintf(fp, "%g %g %g\n", t, 1.0e3*effective_height, 1.0e3*exact_location);
  fflush(fp);
}
/**
### Check the simulation progress.
*/
#include "../src/outputs/progress.h"
event progress(i++) { progress_check(i, t, T_END); }

/**
### Movie maker
*/
#include "view.h"
#include "../src/post_processing/movie_maker.h"
event movie_output(t += DELTA_T) {
  scalar temperature[];
  foreach () {
    temperature[] = f[] * T_L[] + (1. - f[]) * T_V[];
  }
  MESH_ON=0;
  BOX_ON=1;
  movie_maker_i("temperature", i, t);
}

/**
## Results
~~~gnuplot Evolution of the vapor layer thickness
set xlabel "time [s]"
set ylabel "Vapor Layer Thickness [mm]"
set key left top
set size square
set grid

plot "interface_location.dat" every 10 u 1:3 w p ps 2 t "Anal", \
     "interface_location.dat" u 1:2 w l lw 2 t "L6"
~~~
## References

~~~bib
@article{boyd_consistent_2023,
  title = {A consistent volume-of-fluid approach for direct numerical simulation of the aerodynamic breakup of a vaporizing drop},
	copyright = {All rights reserved},
	issn = {0045-7930},
	journal = {Computers \& Fluids},
	author = {Boyd, Bradley and Ling, Yue},
	year = {2023},
	keywords = {DNS, Drop breakup, Vaporization, Volume-of-fluid method},
	pages = {105807},
  doi = {https://doi.org/10.1016/j.compfluid.2023.105807},
  url = {https://www.sciencedirect.com/science/article/pii/S0045793023000324},
}
~~~
*/

