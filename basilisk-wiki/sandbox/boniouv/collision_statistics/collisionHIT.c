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
#define eps0 0.24         // Dissipation rate
#define WeCrit 0.5        // Critical Weber
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
#define A0 pow(eps0/27./sq(Lt),1./3.)   // Forcing constant for HIT
#define k0 (27./2.*sq(Lt*A0))           // Turbulent kinetic energy
#define u0 sqrt(2./3.*k0)               // Turbulent velocity
#define ReTurb (sqrt(k0)*Lt/nu_f)       // Turbulent Reynolds
#define eta pow(pow(nu_f,3.)/eps0,0.25) // Kolmogorov scale
#define teddy (2./3.*k0/eps0)           // Eddy turnover time
#define LTaylor sqrt(10*nu_f*k0/eps0)   // Taylor scale                     
#define ReTaylor (sqrt(45.*A0/nu_f)*Lt) // Taylor Reynolds

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
#define Nbsup 100

// Frequence of droplet tracking
#define DFREQ 40.
#define dtprint (D1/u0/DFREQ)

/**
## Solver
*/
// #include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "../misc/algebra.h"
#include "../linear_forcing/force_turbulence.h"
#include "tension.h" // Needs to be after force_turbulence.h so acceleration events occur before
#include "no-coalescence.h"
#include "tag.h"
#include "collision.h"
#include "../misc/fractions_init.h"

/**
We use the $\lambda_2$ criterion and Basilisk View for visualisation
of vortices. */

#include "lambda2.h"
#include "view.h"
#include "../misc/custom_cmap.h"
#include "../misc/output3d.h"

FILE * fParam;
Drop droplets[Nb+Nbsup]; // Creation of a global Drop instance (defined in Drop.h)
int idump = 0;
int ndrop = Nb;

/** The flow parameters are calculated based on the previously defined non-dimensional parameters. 
 * The tolerance is reduced to 1e-4 and the boundaries are set to be periodic*/
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

/**Initializes the domain with zero velocity and outputs the flow parameters*/
event init (t=0)
{
  int ndrop = 0;
  bool dfound = false;
  int isearch = 0;
  while (!dfound && (isearch<10)) {
    if (restore(file = "../dump-hit")) {
      coord Positions[Nb];
      double Radii[Nb];
      for (int j = 0; j < Nb; j++){
        Radii[j] = R1;
      }
      generate_drops_pos(Positions,Radii,Nb);
      foreach() {
        coord center = POS;
        f[] = refine_frac(center, Positions, Radii, Nb, Delta, 1);
      }
    }
    else if (restore(file = "../restart/dump-restart")) {
      idump = restore_drops(droplets,"../restart/drop_dump.csv");
    }
    else {
      coord Positions[Nb];
      double Radii[Nb];
      for (int j = 0; j < Nb; j++){
        Radii[j] = R1; 
      }
      generate_drops_pos(Positions,Radii,Nb);
      foreach() {
        coord center = POS;
        f[] = refine_frac(center, Positions, Radii, Nb, Delta, 1);
        u.x[] = (cos(y) + sin(z));
        u.y[] = (sin(x) + cos(z));
        u.z[] = (cos(x) + sin(y));
      }
    }
    ndrop = update_droplets (interfaces,droplets);
    dfound = (ndrop == Nb);
    isearch += 1;
    fprintf(stderr,"Loop %d of searching droplets...\n",isearch);
  }
  assert(dfound);
  fprintf(stderr,"Simulation can now start!\n");
}

event track_droplets(t <= TMAX; t += dtprint){
  ndrop = update_droplets (interfaces,droplets);
  // save the datas
  if (pid() == 0) {
    print_droplets (droplets,ndrop,t);
  }
}

event save_dump(t <= TMAX; t += 2*teddy){
  char name[80];
  sprintf(name, "dump-%d", idump);
  dump (name);
  sprintf(name, "../field/u-%d.dat", idump);
  output_field_3d ({f,u.x, u.y, u.z}, fopen (name, "w"), n = 128, linear = true);
  scalar l2[];
  lambda2(u, l2);
  view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
        tx = 0.0122768, ty = 0.0604286,
        width = 1000,
        height = 1000,
        bg = {1,1,1});
  for (scalar s in interfaces){
    draw_vof(s.name, color = "nsc_coloring", map = cool_warm, min = 0, max = 5, lw = 2.);
  }
  squares ("l2", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2.);
  squares ("l2", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2., n = {1,0,0});
  squares ("l2", linear = true, map = gray, min = -30, max = 0, alpha = -Ls/2., n = {0,1,0});
  sprintf (name, "../visu/pop_nsc_%d.png",idump);
  save (name);
  idump += 1;
}

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

#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  adapt_wavelet ((scalar *){f,u}, (double[]){1e-3,uemax,uemax,uemax}, LEVEL);
}
#endif

event save_restart(t <= TMAX; t += teddy){
  dump_drops (droplets,ndrop,t,"drop_dump.csv");
  int nvof = list_len(interfaces);
  FILE * para = fopen("number_of_vof.txt","w");
  fprintf(para,"%d",nvof);
  fclose(para);
  dump ("dump-restart");
}

event stop(t=TMAX){
  dump ("dump-restart");
  return 1;
}