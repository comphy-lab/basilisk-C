/**
# Forced isotropic turbulence in a triply-periodic box
*/

/**
## Parameters of the problem

Space, time and output parameters
*/
#define LEVEL 4             // Refinement level
#define Dx (Ls/(1<<LEVEL))  // Minimum mesh size
#define TDUMP (50*teddy)    // Time to save a dump
#define TMAX (100*teddy)    // Time of the simulation
#define OMOVIE 0            // Activate movie 
#define TMOVIES (0.1*teddy) // Movie framerate

/**
Forcing parameters
*/
#define TFABC 0    // Force with ABC approach (override linear forcing)
#define TFMETH 1   // Linear forcing method
#define TFRHSNUM 1 // Approximation of budget for control (TFMETH 4,5,6 or 7)
#define TFTAU 67.  // Relaxation time (TFMETH 4,5,6 or 7)
#define TFC1 1.    // Only for hybrid methods (TFMETH 6 or 7)
#define TFC2 1.    // Only for hybrid methods (TFMETH 6 or 7)

/**
Input parameters
*/
#define ReTaylor 40.      // Taylor scale
#define Ls (2.*pi)        // Domain length for periodic box
#define dLs 0.19          // Ratio between domain and integral length
#define rho_f 1.          // Carrier phase density
#define mu_f 0.005        // Carrier phase dynamic viscosity
#define nu_f (mu_f/rho_f) // Carrier phase kinematic viscosity

/**
Computation of turbulent quantities
*/ 
#define Lt (dLs*Ls)                     // Integral scale
#define A0 (sq(ReTaylor/Lt)*nu_f/45.)   // Forcing constant for HIT
#define k0 (27./2.*sq(Lt*A0))           // Turbulent kinetic energy
#define u0 sqrt(2./3.*k0)               // Turbulent velocity
#define eps0 (27.*sq(Lt)*pow(A0,3))     // Dissipation rate
#define ReTurb (sqrt(k0)*Lt/nu_f)       // Turbulent Reynolds
#define eta pow(pow(nu_f,3.)/eps0,0.25) // Kolmogorov scale
#define teddy (2./3.*k0/eps0)           // Eddy turnover time
#define LTaylor sqrt(10*nu_f*k0/eps0)   // Taylor scale                     

#define LtABC Lt                        // Forcing scale for ABC flow
#define uABC u0                         // Forcing amplitude for ABC flow

/**
## Solver
*/
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "../misc/algebra.h"
/**
We need to initialize two-phase flow quantities used in force_turbulence.h */
face vector av[];
scalar f[], * interfaces = {f};
double rho1, rho2, mu1, mu2;
#include "force_turbulence.h"

/**
We use the $\lambda_2$ criterion and Basilisk View for visualisation
of vortices. */

#include "lambda2.h"
#include "view.h"
#include "../misc/custom_cmap.h"

FILE * fParam;

int main (int argc, char * argv[])
{
  const face vector muc[] = {mu_f,mu_f,mu_f};
  mu = muc;
  rho1 = rho_f;
  mu1 = mu_f;
  rho2 = rho_f;
  mu2 = mu_f;
  a = av;
  L0 = Ls [0];
  DT = HUGE [0];
  init_grid(1<<(LEVEL));
  origin(-L0/2.,-L0/2.,-L0/2.);
  TOLERANCE = 1e-8 [*];
  foreach_dimension()
    periodic (right);
  run();
}

/**
## Initial conditions

The initial condition is an ABC flow with energy
and wavenumber coherent with the integral length 
of the target turbulence. */

event init (i = 0) {
  if (!restore(file = "../dump-hit")) {
    double kappa0 = 1., uf0 = 1.;
    kappa0 = 10.*pi/Ls;
    uf0 = u0;
    foreach() {
      f[] = 0.;
      u.x[] = uf0*(cos(kappa0*y) + sin(kappa0*z));
      u.y[] = uf0*(cos(kappa0*z) + sin(kappa0*x));
      u.z[] = uf0*(cos(kappa0*x) + sin(kappa0*y));
    }
  }
  fParam = fopen("params.csv","w");
  fprintf(fParam,"rho_f,mu_f,k0,eps0,Re0,eta,LTaylor,teddy\n");
  fprintf(fParam,"%g,%g,%g,%g,%g,%g,%g,%g",
                  rho_f,mu_f,k0,eps0,ReTaylor,eta,LTaylor,teddy);
  fclose(fParam);
}

/**
We generate a movie of the vortices. 
*/

#if OMOVIE
event movies(t<=TMAX; t+=TMOVIES){
  scalar l2[];
  lambda2(u, l2);
  view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
        tx = 0.0122768, ty = 0.0604286,
        width = 1000,
        height = 1000,
        bg = {1,1,1});
  squares ("u.z", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2.);
  squares ("u.x", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2., n = {1,0,0});
  squares ("u.y", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2., n = {0,1,0});
  isosurface ("l2", -20);
  save ("movie.mp4");
}
#endif

event dump_restart(t=TDUMP){
  dump ("dump-hit");
}

event stop(t=TMAX){
  return 1;
}
