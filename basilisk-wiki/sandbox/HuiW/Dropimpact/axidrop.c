/**
Drop impact onto the deep pool of the same liquid at high speed. It is motivated by large raindrops impacting on the ocean surface near terminal speed, as experimentally studied by [David Murphy JFM 2015](https://doi.org/10.1017/jfm.2015.431). 

Here we numerically reproduce and investigate such high-energy splashing dynamics. You can find some recent numerical results in these publications: 

* Axisymmetric simulation: coming soon

* Three-dimensional simulation: Hui Wang, Shuo Liu, Annie-Claude Bayeul-Lainé, David Murphy, Joseph Katz, and Olivier Coutier-Delgosha. [Analysis of high-speed drop impact onto deep liquid pool](https://doi.org/10.1017/jfm.2023.701). Journal of Fluid Mechanics, 972:A31, October 2023.
*/


/**Header files included for this simulation. We solve the incompressible Navier-Stokes equations with surface tension in axisymmetric configurations.*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"                           
#include "reduced.h"                           
#include "view.h"     

/**Pre-defined macros for fluid properties and flow field. Here we use the air-water configuration.*/  
#define L_Ref 0.01  
#define DB (0.004/L_Ref)                                  
#define y_D (0.)                           
#define x_D (0.6*DB)                                              
#define z_D (0.)                                             
#define Drop(x,y) (sq((x-x_D)*4/3.8)+sq((y-y_D)*4/4.3)) 

#define DOMAIN (8.)                                
#define RHO_SW (1.)
#define RHO_A  (1.3/1018.3)
#define MU_SW (DB/Re)
#define MU_A  (MU_SW*1.8e-5/1.e-3)
#define SIGMA_SW (DB/We)                
#define Gravity (1./(36.*36.*DB)) 

#define drop(x,y) (sq(DB/2) - Drop(x,y))
#define pool(x) (-x)
#define epi (0.05*DB)
#define rfdrop(x,y) ((sq(DB/2-epi) < Drop(x,y)) && (sq(DB/2+epi) > Drop(x,y)))
#define rfpool(x) ((-epi < x) && (epi > x))

int MAXLEVEL = 15;
double fErr = 1e-4;          
double uErr = 1e-2;

/**Boundary conditions*/
u.n[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);


/**Main fuction. Gravity is included using reduced formula.*/
int main () {
  size (DOMAIN);
  origin (-DOMAIN/2., 0.);
  init_grid (pow(2.0, MINLEVEL));

  rho1 = RHO_SW;           
  rho2 = RHO_A;                                              
  mu1 = MU_SW;           
  mu2 = MU_A;
  f.sigma = SIGMA_SW;

  G.x = -Gravity;
  run();
}   

/**We define the initial flow field. The drop is initialsed above the pool surface, which allows the consideration of air disk entrapment.*/

event init (t = 0) {
    scalar f4[], f6[];

    refine ((rfdrop(x,y) || rfpool(x)) && (level < MAXLEVEL));
  
    fraction (f4, pool(x));
    fraction (f6, drop(x,y));
    foreach() {
      f[] = f4[] + f6[];
      u.x[] = -f6[];
    }
}

/**Adaptive mesh refinement*/
event adapt (i++) { 
  adapt_wavelet ((scalar *){f, u.x, u.y}, (double[]){fErr, uErr, uErr}, MAXLEVEL);

/**Simulation ends at MAXTIME*/
event end (t = MAXTIME) {                                         
  printf("i=%d t=%g\n",i,t);
}