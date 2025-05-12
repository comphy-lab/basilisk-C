/**
# Binary collision test case with collision statistics output
*/

/**
## Parameters of the problem

Space, time and output parameters
*/
#define LEVEL 7             // Refinement level
#define Dx (Ls/(1<<LEVEL))  // Minimum mesh size
#define teddy 0.1
#define TMAX (5.)    // Time of the simulation
#define OMOVIE 0            // Activate movie 
#define TMOVIES (0.1) // Movie framerate
#define RESTORE 0           // Restart? Put 1

/**
Input parameters
*/
#define WeCrit 0.5        // Critical Weber
#define Ls (2.*pi)        // Domain length for periodic box
#define u0 1.0            // Approaching velocity
#define rho_f 1.          // Carrier phase density
#define mu_f 0.005        // Carrier phase dynamic viscosity
#define nu_f (mu_f/rho_f) // Carrier phase kinematic viscosity
#define rhoRatio 1.       // Ratio between dispersed and carrier phase density
#define muRatio 1.        // Ratio  between dispersed and carrier phase viscosity
#define Nb 2              // Number of droplets

/**
Computation of dispersed phase properties
*/ 
#define rho_d (rho_f*rhoRatio)       // Dispersed phase density
#define mu_d (mu_f*muRatio)          // Dispersed phase dynamic viscosity
#define D1 (0.2*Ls)                  // Droplet diameter
#define SIG (rho_f*sq(u0)*D1/WeCrit) // Surface tension
#define R1 (0.5*D1)                  // Radius 1
#define VOL (4./3.*pi*pow(R1,3))     // Droplet volume
#define SURF (4.*pi*pow(R1,2))       // Droplet surface
#define PHI (VOL*Nb/pow(Ls,3))       // Volume concentration
#define Nbsup 100

// Frequence of droplet tracking
#define DFREQ 40.
#define dtprint (D1/u0/DFREQ)

/**
## Solver
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "../misc/algebra.h"
#include "tension.h"
#include "no-coalescence.h"
#include "tag.h"
#include "collision.h"
#include "../misc/fractions_init.h"

#include "view.h"
#include "../misc/custom_cmap.h"

Drop droplets[Nb+Nbsup];

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

event track_droplets(t <= TMAX; t += dtprint){
  int ndrop = update_droplets (interfaces,droplets);
  if (pid() == 0) {
    print_droplets (droplets,ndrop,t);
  }
  int nvof = list_len(interfaces);
  FILE * para = fopen("number_of_vof.txt","w");
  fprintf(para,"%d",nvof);
  fclose(para);
}

/**Initializes the domain with zero velocity and outputs the flow parameters*/
event init (t=0)
{
  coord Positions[2];
  double Radii[2];
  for (int j = 0; j < 2; j++){
    Radii[j] = R1; 
  }
  Positions[0].x = 0.75*D1;
  Positions[0].y = 0.;
  Positions[0].z = 0.;
  Positions[1].x = -0.75*D1;
  Positions[1].y = 0.;
  Positions[1].z = 0.;
  foreach() {
    coord center = POS;
    f[] = refine_frac(center, Positions, Radii, 2, Delta, 1);
    u.x[] = -f[]*0.5*u0*sign(x);
  }
}

/** Outputs videos of the velocity
*/

event output (t = TMOVIES; t <= TMAX; t += TMOVIES)
{
  view (fov = 44, 
        camera = "front", 
        width = 800, 
        height = 800,
        bg = {1,1,1});
  for (scalar s in interfaces){
    draw_vof(s.name, color = "nsc_coloring", map = cool_warm, min = 0, max = 3, lw = 2.);
  }
  squares(color = "u.x");
  save ("movie.mp4");
}