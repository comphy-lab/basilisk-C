/**
# Bubble growth in a superheat liquid
Section 4.3 of [Boyd et al, 2023](#boyd_consistent_2023).
*/

/**

![Bubble cross-section](bubble_growth/temperature.mp4)(width="640" height="640")
![Bubble growth](bubble_growth/vof_movie.mp4)(width="640" height="640")

*/

/** 
### 3D test (alternatively could be 2D axi)
*/
//#define AXI 0
#include "grid/octree.h"

/**
## Setup the problem
### Initial bubble radius, domain size, and simulation duration
*/

#define R0 0.12
#define L (5.0 * R0)
#define T_END 2.0
#define BUBBLE_VOL

/**
### Fluid properties for saturatate water liquid and vapor
*/
#include "../src/properties/liq_gas.h"

/**
### Analytical solution to the bubble growth 
*/
#define T_inf (T_sat + 2.0)
#include "../src/theory_evap/bubble_growth_calc.h"
/**
### Include header for NS, phase-change, gravity, etc
*/
#include "../src/01_vaporization/evap_include.h"

/**
### Boundary conditions for the outflow (top, right, front)
*/

T_L[top] = dirichlet(T_inf);
T_L[right] = dirichlet(T_inf);
T_L[front] = dirichlet(T_inf);

fE_L[top] = dirichlet(rhocp_L * T_inf);
fE_L[right] = dirichlet(rhocp_L * T_inf);
fE_L[front] = dirichlet(rhocp_L * T_inf);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

u.n[front] = neumann(0.);
p[front] = dirichlet(0.);

/**
### Main function.
*/
int main(int argc, char* argv[]) {
  MIN_LEVEL = 5;
  LEVEL = 6;
  MAX_LEVEL = 6;

  setup_evap();

  size(L);
  origin(0.0, 0.0, 0.0);
  N = 1 << LEVEL;
  
  init_grid(N);

  G.x = 0.;
  Z.x = 0.;

  run();
}

/**
### Initialize the VOF, liquid and vapor temperature.
*/
#define sphere(x, y, z, R) (R - sqrt(sq(x) + sq(y) + sq(z)))

double beta_GLOBAL = 0.0;
double t_GLOBAL = 0.0;
event init(i = 0) {
  if (!restore (file = "restart")) {
    CFL = 0.5;
#if TREE
    refine(level < MAX_LEVEL && sphere(x, y, z, (0.5*R0)) < 0. && sphere(x, y, z, (2.0*R0)) > 0.);
#endif
    fraction_LS(f, -sphere(x, y, z, R0));

    double beta_ = find_beta();
    double time_ = time_given_radius(R0,beta_);

    beta_GLOBAL = beta_;
    t_GLOBAL = time_;

    if (pid()==0){
      printf("start_time_  %g\n", time_);
    }

    foreach() {
      foreach_dimension() {
        u.x[] = 0.;
      }

      double T_vap = T_sat;
      double radial_dist = sqrt(sq(x) + sq(y) + sq(z));
      double T_liq = T_liq_func(time_, radial_dist, beta_);

      T_V[] = T_vap;
      T_L[] = T_liq;
    }
  }

  init_Energy();
}

/**
### Adaptive mesh refinement.
*/
#define femax 1.0e-5
#define Temax 1.0e-2
#define uemax 1.0e-2
#include "../src/01_vaporization/adapt_evap.h"

/** 
### Bubble volume calcuation
*/
#ifdef BUBBLE_VOL
#define ffrac (1.0 - f[])
#else // droplet volume
#define ffrac f[]
#endif
double bubble_volume() {
     double vb = 0.;
     foreach (reduction(+:vb)){
         double dvb = ffrac * dv();
         vb += dvb; // volume of the bubble
     }
 #if AXI
     vb *= 2. * pi;
 #endif  // AXI
    // only modelling a 1/8 of a bubble
    return 8. * vb;
}

/**
### Outputs the interface location and the analytical interface location.
*/
event outputs(t = 0.; t += DELTA_T; t < T_END) {
  double effective_radius = pow(3.0/4.0*bubble_volume()/pi, 1.0/3.0);

  static FILE* fp_radius = fopen("radius_in_time.dat", "w");

  double exact_radius = bubble_radius(t+t_GLOBAL, beta_GLOBAL);
  fprintf(fp_radius, "%g %g %g\n", t, effective_radius, exact_radius);
  fflush(fp_radius);
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
event moive_output(t += DELTA_T) {
  scalar temperature[];
  foreach () {
    temperature[] = f[] * T_L[] + (1. - f[]) * T_V[];
  }
  TX = 0.0;
  TY = 0.0;
  MESH_ON=0;
  BOX_ON=1;
  vof_movie_3D_mirror_x_y_z();
  movie_maker_3D_mirror_x_y("temperature", "front", t);
}

/**
## Results
~~~gnuplot Evolution of the bubble radius
set xlabel "time [s]"
set ylabel "Bubble radius [m]"
set key left top
set size square
set grid

plot "radius_in_time.dat" every 10 u 1:3 w p ps 2 t "Anal", \
     "radius_in_time.dat" u 1:2 w l lw 2 t "L6"
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
