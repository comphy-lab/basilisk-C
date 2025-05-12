/**
# Turbulent emulsion promoted by linear forcing in triply-periodic box
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
#define RESTORE 0           // Restart? Put 1

/**
Forcing parameters
*/
#define TFABC 0    // Force with ABC approach (override linear forcing)
#define TFMETH 1   // Linear forcing method
#define TFRHSNUM 1 // Approximation of budget for control (TFMETH 4,5,6 or 7)
#define TFTAU 67.  // Relaxation time (TFMETH 4,5,6 or 7)
#define TFC1 1.    // Only for hybrid methods (TFMETH 6 or 7)
#define TFC2 1.    // Only for hybrid methods (TFMETH 6 or 7)
#define UMEAN 2    // Way to deal with the kinetic energy exponential growth

/**
Input parameters
*/
#define ReTaylor 40.      // Taylor scale
#define WeCrit 25.        // Critical Weber
#define Ls (2.*pi)        // Domain length for periodic box
#define dLs 0.19          // Ratio between domain and integral length
#define rho_f 1.          // Carrier phase density
#define mu_f 0.005        // Carrier phase dynamic viscosity
#define nu_f (mu_f/rho_f) // Carrier phase kinematic viscosity
#define rhoRatio 1.       // Ratio between dispersed and carrier phase density
#define muRatio 1.        // Ratio  between dispersed and carrier phase viscosity
#define Nb 1              // Number of droplets
#define PHI 0.1           // Volume concentration

/**
Computation of turbulent quantities
*/ 
#define Lt (dLs*Ls)                    // Integral scale
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
Computation of dispersed phase properties
*/ 
#define rho_d (rho_f*rhoRatio)       // Dispersed phase density
#define mu_d (mu_f*muRatio)          // Dispersed phase dynamic viscosity
#define VOL (PHI*pow(Ls,3))          // Droplet volume
#define R1 pow(3./4./pi*VOL,1./3.)   // Radius 1
#define D1 (2*R1)                    // Droplet diameter
#define SIG (2./3.*rho_f*pi*k0/WeCrit)
#define etaRatio (D1/eta)            // Ratio between Kolmogorov scale and droplet diameter
#define SURF (4.*pi*pow(R1,2))       // Droplet surface

/**
## Solver
*/
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "../misc/algebra.h"
#include "force_turbulence.h"
#include "tension.h" // Needs to be after force_turbulence.h so acceleration events occur before
#include "../misc/fractions_init.h"

/**
We use the $\lambda_2$ criterion and Basilisk View for visualisation
of vortices. */

#include "lambda2.h"
#include "view.h"
#include "../misc/custom_cmap.h"
#include "../misc/output3d.h"

FILE * fParam;

int main(int argc, char** argv){  
  L0 = Ls [0];
  DT = HUGE [0];
  init_grid(1<<(LEVEL));
  origin(-L0/2.,-L0/2.,-L0/2.);
  rho1 = rho_f;
  mu1 = mu_f;
  rho2 = rho_f;
  mu2 = mu_f;
  f.sigma = SIG;
  TOLERANCE = 1e-5 [*];
  foreach_dimension()
    periodic(right);
  run();
}

/**
## Initial conditions

The initial condition is a restart from a [single-phase HIT](SP.c). 
A single droplet is generated at the center of the domain. */

event init (t=0)
{
  if (restore(file = "../dump-hit")) {
    coord Positions[1];
    double Radii[1];
    Radii[0] = R1;
    Positions[0] = coord_null; 
    foreach() {
      coord center = POS;
      f[] = refine_frac(center, Positions, Radii, 1, Delta, 1);
    }
  }
  else {
    coord Positions[1];
    double Radii[1];
    Radii[0] = R1;
    Positions[0] = coord_null; 
    foreach() {
      coord center = POS;
      f[] = refine_frac(center, Positions, Radii, 1, Delta, 1);
      u.x[] = (cos(y) + sin(z));
      u.y[] = (sin(x) + cos(z));
      u.z[] = (cos(x) + sin(y));
    }
  }
  fParam = fopen("params.csv","w");
  fprintf(fParam,"Nb,rho_f,rho_d,mu_f,mu_d,sig,R,V,S,Phi,");
  fprintf(fParam,"k0,eps0,Re0,eta,LTaylor,teddy\n");
  fprintf(fParam,"%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,",
                  Nb,rho_f,rho_d,mu_f,mu_d,SIG,R1,VOL,SURF,PHI);
  fprintf(fParam,"%g,%g,%g,%g,%g,%g",
                  k0,eps0,ReTaylor,eta,LTaylor,teddy);
  fclose(fParam);
}

/**
## Stop the simulation

Dump the solution and output the velocity and volume
fraction fields */

event stop(t=TMAX){
  dump ("dump-restart");
  output_field_3d ({u.x, u.y, u.z, f}, fopen ("u.dat", "w"), linear = true);
  return 1;
}