/**
# 2D Film boiling problem
Section 4.4 of [Boyd et al, 2023](#boyd_consistent_2023).
*/

/**

![Film boiling](film_boiling/temperature.mp4)(width="640" height="640")

*/

/**
## Setup the problem
### Domain size, and simulation duration
*/
double lambda_RT = 0.0;
#define L (lambda_RT)  // size of the box
//#define T_END (0.528)
#define T_END (0.4)
/**
### Fluid properties and parameters
*/
#include "../src/properties/film_boiling_case.h"
#define T_inf (5. + T_sat)
#define grav_mag 9.81

/**
### Include header for NS, phase-change, gravity, etc
*/
#include "../src/01_vaporization/evap_include.h"

/**
### Boundary conditions for the outflow (top) and wall (bottom)
*/
T_V[bottom] = dirichlet(T_inf);
T_L[top] = dirichlet(T_sat);

fE_V[bottom] = dirichlet(rhocp_V * T_inf);
fE_L[top] = dirichlet(rhocp_L * T_sat);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);

u.n[bottom] = dirichlet(0.);
p[bottom] = neumann(0.);

/**
### Main function.
*/
int main(int argc, char* argv[]) {

    lambda_RT = 2.0 * pi * sqrt(3. * sigma_lv / (grav_mag * (rho_L - rho_V)));
    fprintf(stdout, "lambda_RT %g\n", lambda_RT);

    MIN_LEVEL = 5;
    LEVEL = 8;
    MAX_LEVEL = 8;

    setup_evap();

    size(L);
    origin(-0.5*L0, 0.0);
    N = 1 << LEVEL;
    init_grid(N);

    G.y = -grav_mag;

    run();
}

/**
### Initialize the VOF, liquid and vapor temperature.
*/
#define RT_y_func(x) ((4. + cos(2. * pi * x / lambda_RT)) * lambda_RT / 128.)
#define RT_func(x, y) (y - RT_y_func(x))

double T_vap_func(double x, double y){
    double y_max = RT_y_func(x);
    double T_y = T_inf - y/y_max*(T_inf-T_sat);
    T_y = clamp(T_y, T_sat, T_inf);
    return T_y;
}

event init(i = 0) {
    if (!restore (file = "restart")) {
        CFL = 0.2;
#if TREE
        refine(level < MAX_LEVEL && RT_func(x, y - dR_refine) < 0. &&
               RT_func(x, y + dR_refine) > 0.);
#endif
        fraction_LS(f, RT_func(x, y));

        foreach () {
            foreach_dimension() {
                u.x[] = 0.;
            }

            double T_vap = T_vap_func(x, y);
            
            T_V[] = (1.0 - f[]) * T_vap + f[] * T_sat;
            T_L[] = T_sat;
        }
        foreach_face() uf.x[] = 0.;
    }

    init_Energy();
}

/**
### Adaptive mesh refinement.
*/
#define REFINE_VAPOR_REGION
#define femax 1.0e-5
#define Temax 1.0e-2
#define uemax 1.0e-2
#include "../src/01_vaporization/adapt_evap.h"


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
  TX = 0.05;
  MESH_ON=0;
  BOX_ON=1;
  movie_maker_i("temperature", i, t);
}

/**

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

