/**
# Sucking Problem
Section 4.2 of [Boyd et al, 2023](#boyd_consistent_2023).
*/
/**

![Sucking Flow](sucking_problem/temperature.mp4)(width="640" height="640")

*/
/**
## Setup the problem
### Domain size, initial vapor layer thickness, and simulation duration
*/
#define H0 0.02
#define L (1.0)
#define T_END 2.0

/**
### Fluid properties for saturatate water liquid and vapor
*/
#include "../src/properties/liq_gas.h"

#define T_inf (T_sat + 2.0)

/**
### Include header for NS, phase-change, gravity, etc
*/
#include "../src/01_vaporization/evap_include.h"

/**
### Analytical solution to the sucking problem
*/
double vel_i_max = 1.;
#include "../src/theory_evap/sucking_calc.h"

/**
### Boundary conditions for the wall (bottom) and outflow (top)
*/

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
T_L[top] = dirichlet(T_inf);
fE_L[top] = dirichlet(rhocp_L * T_inf);

u.n[bottom] = dirichlet(0.);
p[bottom] = neumann(0.);
T_V[bottom] = dirichlet(T_sat);
fE_V[bottom] = dirichlet(rhocp_V * T_sat);

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
double beta_GLOBAL = 0.;
double t_sucking_GLOBAL = 0.;
#define plane(x, y, H) (H - y)
event init(i = 0) {
  CFL = 0.1;

#if TREE
  refine(level < MAX_LEVEL && plane(x, y, (H0 - dR_refine)) < 0. &&
         plane(x, y, (H0 + dR_refine)) > 0.);
#endif
  fraction_LS(f, -plane(x, y, H0));

  double beta_ = find_beta();
  double start_time = find_start_time(H0, beta_);

  beta_GLOBAL = beta_;
  t_sucking_GLOBAL = start_time;

  // init interface velocity
  double uy = vel_interface_func(start_time, beta_);
  vel_i_max = uy;

  foreach () {
    foreach_dimension() {
      u.x[] = 0.;
    }
    u.y[] = uy;

    double T_vap = T_sat;
    double T_liq = T_inf;

    T_liq = T_liq_func(start_time, y, beta_);

    T_V[] = T_vap;
    T_L[] = T_liq;
  }

  init_Energy();
}

/**
### Adaptive mesh refinement.
*/
#define femax 1.0e-5
#define Temax (1.0e-5*(T_inf-T_sat))
#define uemax (1.0e-1*vel_i_max)
#include "../src/01_vaporization/adapt_evap.h"

/**
### Outputs the interface location and the analytical interface location.
*/
event outputs(t = 0.; t += DELTA_T; t < T_END) {
  double effective_height;
  scalar fc[];
  foreach ()
    fc[] = 1. - f[];

  effective_height = statsf(fc).sum / L0;

  char name[200];
  sprintf(name, "interface_location.dat");
  static FILE* fp = fopen(name, "w");

  double exact_location = interface_location_func(t+t_sucking_GLOBAL, beta_GLOBAL);
  fprintf(fp, "%g %g %g\n", t, effective_height, exact_location);
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
set ylabel "Vapor Layer Thickness [m]"
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


