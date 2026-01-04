/**
This setup is an axisymmetric simulation of an O₂ bubble in water growing at a solid boundary with fixed mass transfer rate $\dot{m}$.*/
/**
The phase-change solver is taken from -[Genari's sandbox](/sandbox/ggennari/README)
*/
/**
The contact line is pinned at a fixed possition $ap$ by using the pinning function from -[Cipriano's sandbox](/sandbox/ecipriano/src/pinning.h) 
*/


#define F_ERR 1e-10
#define MASS_TRANSFER_RATE (0.1)

#include "grid/quadtree.h"
#include "axi.h"
#define PHASE_CHANGE 1 // phase change
#include "solver/navier-stokes/centered.h"
#include "solver/phase_change/two-phase.h"
#include "solver/phase_change/diffusion.h"
#include "solver/phase_change/phase_change_pure_species.h"
#include "solver/interface/adapt_wavelet_leave_interface.h"
#include "fractions.h"
// contact angle & surface tension
#include "contact.h"
#include "tension.h"
#include "curvature.h"



// htg head file
#include "output_htg.h"

/*
  Two-phase system properties
*/

// Gas-O2
#define rho_g 1.43 // density 
#define mu_g 20.27e-6 // viscosity 
#define Henry_cst 34.375; // Henry's coeff
#define Molar_mass 0.032;// molar mass
// Liquid-H2O
#define rho_l 1000. // density 
#define mu_l 1.0016e-3 // viscosity 
// surface tension
#define f_sigma_ST 0.072;

#define zeta 64. // saturation ratio   // 20kPa=5.05; 35kPa=2.89; 50kPa=2.02.
#define Rb 20e-6 // initial bubble Radius in "m" unit


#define relax_tt 0     // 1: switch on the relaxation time



double theta_int = 130.;




double WIDTH = 0.0001;                 // domain size
double TEND = 1.; // simulation time

int LEVEL = 9; // max. refinement level
int basic_LEV = 5; // basic refinement level


/**
The concentration field is stored in the scalar *a_c*
*/
scalar a_c[];
scalar *species = {a_c};
scalar m_a[];
scalar *m_list = {m_a};
// curvature of bubble interface
//scalar kappa[];


vector h[];


/**
Known the position of the interface $a$ we compute the corresponding
heigth function at the neighbouring cell.*/

foreach_dimension()
static double line_x (Point point, scalar h, double a)
{
  return (h[] == nodata ? nodata : -h[] + 2.*(a-y)/Delta);
}

#define contact_line(theta) contact_line_ (point, neighbor, _s, theta)

double contact_line_ (Point point, Point neighbor, scalar h, double a)
{
  if (neighbor.i != point.i)
    return line_x (point, h, a);
  if (neighbor.j != point.j)
    return line_y (point, h, a);
  assert (false); // not reached
  return 0.;
}

double ap;
h.t[left] = y > 0 ? contact_line (ap) : neumann (0);



int main()
{
  size(WIDTH); // length of domain;

  init_grid (1 << basic_LEV);
  //init_grid(pow(2,basic_LEV));

  rho1 = rho_l;       // liquid phase
  rho2 = rho_g; // gas phase
  mu1 = mu_l;      // liquid phase
  mu2 = mu_g;    // gas phase

  //f.tracers = species; // switch on/off convection term in transport equation
  f.sigma = f_sigma_ST;

  a_c.inverse = false;          // compute the tracer in liquid phase (inverse== false 0); gas phase (inverse==true 1)
  //a_c.D = Diff_sim;             // Diffusion coefficient
  a_c.mol_mass = Molar_mass;         // molar mass
  //a_c.He = Henry_cst;              // Henry's coeff
  a_c.cd = rho2 / a_c.mol_mass; // uniform concentration inside the bubble
  a_c.diff = false;              // we turn on diffusion of the tracer

  ap = Rb * sin(theta0 / 180. * pi);;



  f.height = h;


#if relax_tt
  //double par1, beta1;
  //par1 = (zeta - 1.) * a_c.cd / a_c.He * a_c.mol_mass / rho2;
  //beta1 = (par1 + sqrt(par1 * par1 + 2. * par1)) / 2.;
  //tt = pow(Rb / 2. / beta1, 2) / a_c.D;
  tt=0.01*Diff_exp/Diff_sim;
  //tt=0.01*TEND;
#else
  tt = 0.0; // we switch off the volume change (the bubble must have a constant radius)
#endif

  // solution setup
  TOLERANCE = 1e-7;
  DT = 1.0e-7;
  CFL = 0.5;
  NITERMAX = 5000;
  run();
}

/* in htg file the domain is like this
                right
      -------------------------
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
      -------------------------
                left
*/
/**
 * No slip boundary condition. */
u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);


/**
  We need outflow boundary conditions to let the liquid exit the domain as the bubble grows.
  */
// Bcs, allow outflow at right and top
u.n[right] = neumann(0.); // zeroGradient of normal velocity
p[right] = dirichlet(0.);                       // zeroPressre at cell center
pf[right] = dirichlet(0.);                      // zeroPressure at cell face (?)
u.n[top] = neumann(0.);



//**************2D coordinates**************//
// x_simu =y;
// y_simu =x;
//******************************************//
event init(t = 0)
{
   if (!restore (file = "restart")) {
  double Hc = Rb * cos(theta_int / 180. * pi);
  refine(sq(x - Hc) + sq(y) < sq(1.2 * Rb) && level < LEVEL); // under-saturation
  fraction(f, sq(x - Hc) + sq(y) - sq(1.0 * Rb));
  
}
}


event acceleration(i++)
{
  face vector av = a;
  foreach_face(x)
      av.x[] -= 9.8; //gravity closed here
}


// adaptive mesh refinement
event adapt(i++)
{
  double uemax = 1e-3;
  adapt_wavelet_leave_interface({a_c, u.x, u.y}, {f}, (double[]){1e-3, uemax, uemax, 1e-3}, LEVEL, basic_LEV, 1);
}


event snapshot1(t = 0; t <= TEND; t = t + 5e-6)
//event snapshot1(i = i + 100)
{
  char name[80];
  sprintf(name, "dump/dump-%.7f", t);
  dump(file = name);
}

event snapshot2(t = 0; t <= TEND; t = t + 5e-6)
//event snapshot2(i = i + 100)
{
  char path[] = "htg"; // no slash at the end!!
  char prefix[80];

  sprintf(prefix, "data-%.7f", t);
  output_htg((scalar *){f, u.x, u.y, a_c, p}, (vector *){uf}, path, prefix, i, t);
}

event stop_simulation(t = TEND)
{
  return 1;
}

