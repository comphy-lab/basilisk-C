/** 2D Droplet spreading simulation with Dynamic Contact Angle (DCA) models */
/** This code reproduces the experimental setup described in the paper:*/
/** Kensuke Yokoi, Damien Vadillo, John Hinch, and Ian Hutchings,*/
/** "Numerical studies of the influence of the dynamic contact angle on a drop impacting on a dry surface," Physics of Fluids 21, 072102 (2009). [Link to paper](https://www.researchgate.net/publication/228351471_Numerical_studies_of_the_influence_of_the_dynamic_contact_angle_on_a_drop_impacting_on_a_dry_surface)*/

// ---------------- Grid and Solver Setup ---------------- //
#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
// htg head file
#include "output_htg.h"
// Dynamic wetting headfiles
#include "contact_line_velocity.h"
#include "dca_models.h"

// ---------------- Simulation Parameters ----------------
#define ca 0 // Switch: 1 = fixed contact angle, 0 = none
#define dca 1 // Switch: 1 = dynamic contact angle model enabled

#define F1_ERR 1e-6 // Tolerance for identifying interfacial cells
#define WIDTH 5.0e-3 // Domain width (5 mm)
#define MAX_LEVEL 8 // Maximum AMR refinement level
#define MIN_LEVEL 5 // Minimum AMR refinement level

double TEND = 0.03; // End time [s]
#define Rb 1.14e-3 // Initial bubble radius [m]


// ---------------- Fluid Properties ----------------
// Gas (air)
#define rho_g 1.25 // Density [kg/m^3]
#define mu_g 1.82e-5 // Viscosity [Pa·s]


// Liquid (water)
#define rho_l 1000. // Density [kg/m^3]
#define mu_l 1.0e-3 // Viscosity [Pa·s]


// Surface tension
#define f_sigma_ST 0.072 // [N/m]


// ---------------- DCA Model Parameters ----------------

double theta_int = 0.; // Initial interface angle
double theta0 = 90.; // Equilibrium contact angle [deg]
double theta_adv = 110.; // Advancing contact angle [deg]
double theta_rec = 72.; // Receding contact angle [deg]

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential components. */
vector h[]; // Height function for curvature & boundary conditions
#if ca
h.t[left] = contact_angle(theta0 * pi / 180.); // Static contact angle if ca=1
#endif



/**********DCA models*****************/



/*---------------- Dynamic contact angle model----------------*/
#if dca
// Dynamic variables
double theta_var, ucl2, Ca2, factor; // Contact line velocity, capillary number, orientation
/** Event: Apply dynamic contact angle model at each time step */
event DCA_bc(i++){
  
double output_values[2];
/** Compute contact line velocity and orientation factor*/
contact_line_velocity_and_normal_vector(f, u, WIDTH, MAX_LEVEL , output_values);
ucl2 = output_values[0];
factor = output_values[1];

  /** Compute capillary number: $Ca = \mu u_{cl} / \gamma$*/
Ca2 = mu_l * ucl2 / f_sigma_ST;


// Adjust sign based on direction (advance vs recede)
if (factor >= 0.0)
Ca2 = fabs(Ca2); // Advancing
else
Ca2 = -fabs(Ca2); // Receding
  
  
  
/******Call the model parameters**********/
 ModelParams params;
params.Ratio =  1e6; // HT model: Length ratio: L/l_m
params.xi_mkt = 0.002; // MKT model: xi_mkt = kB_cst*Temp/kappa/lambda^3
params.Cst01 = 1e4; // Pinning slope coefficient
params.mu_liquid = mu_l; 

// Select DCA model (HT / MKT / COMBINED / COMBINED_PIN)
theta_var = dynamic_contact_angle(Ca2, theta0,  theta_adv, theta_rec, params, "COMBINED_PIN");
if (theta_var<0.)
theta_var = 0.0;

// Apply updated dynamic contact angle at wall
h.t[left] = contact_angle(theta_var *pi/180.);

}
#endif

int main()
{

  size(WIDTH); // length of domain;
  //origin(-1,0,-1);
  init_grid(pow(2,MIN_LEVEL));
  /**
  We use a constant viscosity. */
  rho1 = rho_l;       // liquid phase
  rho2 = rho_g; // gas phase
  mu1 = mu_l;      // liquid phase
  mu2 = mu_g;    // gas phase
  G.x = -9.8;   // gravity

// Associate height function with VOF tracer
  f.height = h;
  
  f.sigma = f_sigma_ST;

  DT=1e-7;
  TOLERANCE = 1e-5;
  CFL=0.25;
  NITERMAX = 500;
  run();
}


/**
The initial drop is a quarter of a sphere. */
event init (t = 0)
{
  if (!restore (file = "restart")) {
  double Hc = 0.99*Rb * cos(theta_int / 180. * pi);
  refine( (sq(x-Hc) + sq(y) < sq(1.2*Rb)) && level < MAX_LEVEL); 
  fraction (f,  -(sq(x-Hc) + sq(y) - sq(Rb)));

// Initial spreading velocity
  foreach(){
   u.x[] = -1.0*f[];
  }
  
  }
}



/* in htg file the domain is like this
                right
      x------------------------
      |                       |
      |                       |
      |                       |
      |                       |
bottom|                       |top
 axi  |                       |
      |                       |
      |                       |
      |                       |
      |                       |
      -------------------------y
                left
*/
/* No slip boundary condition. */
u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);



u.n[top] = neumann(0.); // zeroGradient of normal velocity
u.n[right] = neumann(0.); // zeroGradient of normal velocity

p[right] = dirichlet(0.);                     // zeroPressure at cell face (?)
pf[right] = dirichlet(0.);     


/**
We log statistics on the maximum velocity, curvature and volume. If
the standard deviation of curvature falls below $10^{-2}$, we assume
that the steady shape is reached and we stop the calculation. */

// adaptive mesh refinement
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-3, 1e-2,1e-2},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}





event log_simulation(i = i+100)
{
  double vol1 = 0.;
  scalar posy[], posx[], posxcl[]; // interface position
  position(f, posy, {0, 1});       // y-coordinate of interface
  position(f, posx, {1, 0});       // x-coordinate of interface
  position(f, posxcl, {0, 1});     // position of contact line

  foreach (reduction(+
                     : vol1))
  {
    if (f[] <= F1_ERR|| f[] >= 1. - F1_ERR)
    {
      posx[] = nodata;
      posy[] = nodata;
    }

    if (f[] <= F1_ERR || f[] >= 1. - F1_ERR || x > WIDTH / pow(2, MAX_LEVEL))
    {
      posxcl[] = nodata; // If vertical distance too far from surface then posycl nodata
    }

    vol1 += (1. - f[]) * dv();
  }

  fprintf(stderr, " %g %g %g %g %g %g %g %g\n", t, vol1, statsf(posx).max, statsf(posy).max, statsf(posxcl).max, ucl2, Ca2, theta_var);
}


event snapshot1(t = 0; t <= TEND; t = t + 1e-4)
{
  char name[80];
  sprintf(name, "dump/dump-%.7f", t);
  dump(file = name);
}

event snapshot2(t = 0; t <= TEND; t = t + 1e-4)
{
  char path[] = "htg"; // no slash at the end!!
  char prefix[80];

  sprintf(prefix, "data-%.7f", t);
  output_htg((scalar *){f,p,u.x,u.y}, (vector *){uf}, path, prefix, i, t);
}


event stop_simulation(t = TEND)
{
  return 1;
}


