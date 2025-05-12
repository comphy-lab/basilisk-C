/**
# Droplet-laden HIT promoted by linear forcing in triply-periodic box
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
#define etaRatio 20.      // Ratio between Kolmogorov scale and droplet diameter
#define Nb 1              // Number of droplets

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
#define rho_d (rho_f*rhoRatio)      // Dispersed phase density
#define mu_d (mu_f*muRatio)         // Dispersed phase dynamic viscosity
#define D1 (eta*etaRatio)           // Droplet diameter
#define SIG (2.*rho_f*pow(D1,5./3.)*pow(eps0,2./3.)/WeCrit) // Surface tension
#define R1 (0.5*D1)                 // Radius 1
#define VOL (4./3.*pi*pow(R1,3))    // Droplet volume
#define SURF (4.*pi*pow(R1,2))      // Droplet surface
#define PHI (VOL*Nb/pow(Ls,3))      // Volume concentration

/**
## Solver
*/
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "../misc/algebra.h"
#include "force_turbulence.h"
#include "tension.h" // Needs to be after force_turbulence.h so acceleration events occur before
#include "no-coalescence.h"
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
Droplets are generated at random positions until the prescribed
number of droplets `Nb` is reached.

First, a routine is used to check if a new droplet does not 
overlap with existing droplets */

int check_all_dist(const coord * Position,const  double * Radii, int j){
  for (int k = 0; k < j; k++){
    if(dist_perio(Position[j],Position[k]) < ((Radii[j]+Radii[k])+5.*sqrt(dimension)*Dx))
      return 0;
  }
  return 1;
}

event init (t=0)
{
  coord Positions[Nb];
  double Radii[Nb];
  for (int j = 0; j < Nb; j++){
    Radii[j] = R1; 
  }
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  srand(time(NULL));
  int j = 0;
  int ntry = 0;
  fprintf (stderr, "Generating random droplet positions...\n");
  while (j < Nb) {
    ntry++;
    assert (ntry < 100000*Nb);
    double Rm = (double) RAND_MAX + 1.0;
    Positions[j].x = ((double) rand()/Rm) * Ls - Ls/2;
    Positions[j].y = ((double) rand()/Rm) * Ls - Ls/2;
    Positions[j].z = ((double) rand()/Rm) * Ls - Ls/2;
    if (check_all_dist(Positions, Radii, j)) {
      j++;
      ntry = 0;
    }
  }
  if (restore(file = "../dump-hit")) {
    foreach() {
      coord center = POS;
      f[] = refine_frac(center, Positions, Radii, Nb, Delta, 1);
    }
  }
  else {
    foreach() {
      coord center = POS;
      f[] = refine_frac(center, Positions, Radii, Nb, Delta, 1);
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
## Restoring number of fields for no-coalescence

If a restart is performed using no-coalescence.h, 
it is required to read the number of fields from 
a file (here ../restart/number_of_vof.txt). */

#ifndef RESTORE
  #define RESTORE 0
#endif

#if RESTORE
event defaults (i=0)
{
  FILE * fp = fopen("../restart/number_of_vof.txt", "r");
  assert (fp);
  assert(fscanf(fp, "%d", &length_of_interfaces) == 1);
}
#endif

/**
## Stop the simulation

When stopping the simulation, write nvof in a file and
dump the solution. Also output the velocity and volume
fraction fields */

event stop(t=TMAX){
  int nvof = list_len(interfaces);
  FILE * para = fopen("number_of_vof.txt","w");
  fprintf(para,"%d",nvof);
  fclose(para);
  dump ("dump-restart");
  output_field_3d ({u.x, u.y, u.z, f}, fopen ("u.dat", "w"), linear = true);
  return 1;
}