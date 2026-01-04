
/**
This setup is an axisymmetric simulation of an O₂ bubble in oversaturated water growing at a solid boundary with dynamic wetting*/
/**
The phase-change solver is taken from -[Genari's sandbox](/sandbox/ggennari/README)
*/

#define F_ERR 1e-10

#include "grid/quadtree.h"
#include "axi.h"
#define PHASE_CHANGE 1 // phase change
#include "solver/navier-stokes/centered.h"
#include "solver/phase_change/two-phase.h"
#include "solver/phase_change/diffusion.h"
#include "solver/phase_change/phase_change_pure_species.h"
#include "solver/interface/adapt_wavelet_leave_interface.h"
// contact angle & surface tension
#include "contact.h"
#include "tension.h"
#include "curvature.h"
/** Dynamic wetting headfiles */
#include "contact_line_velocity.h"
#include "dca_models.h"


// htg head file
#include "output_htg.h"

/*
  Two-phase system properties
*/

// O2
//#define RHOR 699.3007 // density ratio
//#define MUR 71.9178 // viscosity ratio
// #define Sc 0.0525 //Schmidt number Sc = mu_l / rho_l / Diff
#define zeta 8. // saturation ratio


#define Rb 10e-6 // initial bubble Radius in "m" unit

#define ca 1           // is 1 when need to apply contact angle
#define dca 1          // dynamic contact angle
#define relax_tt 1     // 1: switch on the relaxation time
#define laplace_term 0 // 1 switch on the effect of laplace pressure on saturation concentration

double theta0 = 90.; // equilibrium contact angle
double theta_int = 90.;

double theta_var;
double Ls = 500e-9; // slip length

double WIDTH = 20. * Rb;
double TEND = 1.0e-3;

int LEVEL = 10; // max. refinement level
int basic_LEV = 6;

double t1;
double t2;

double bubble_mass = 0.;

/**
The concentration field is stored in the scalar *a_c*
*/
scalar a_c[];
scalar *species = {a_c};
scalar m_a[];
scalar *m_list = {m_a};


/** ---------------- DCA Model Parameters ---------------- */

double theta_int = 90.; // Initial interface angle
double theta0 = 90.; // Equilibrium contact angle [deg]
double theta_adv = 110.; // Advancing contact angle [deg]
double theta_rec = 72.; // Receding contact angle [deg]



//** Define  BC for height function */
#if ca
vector h[];
h.t[left] = contact_angle(theta0 * pi / 180.);
#endif

/**  Dynamic contact angle */
#if dca
double CL_0;
double ucl2, Ca2;
double vol = 0.;

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
  
  
  
/** ****Call the model parameters********* */
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
  init_grid(pow(2,basic_LEV));

#if ca
  f.height = h;
#endif

  rho1 = 1000.;       // liquid phase
  rho2 = 1.43; // gas phase
  mu1 = 1.0016e-3;      // liquid phase
  mu2 = 20.27e-6;    // gas phase

  f.tracers = species; // switch on/off convection term in transport equation
  f.sigma = 0.072;

  a_c.inverse = false;          // compute the tracer in liquid phase (inverse== false 0); gas phase (inverse==true 1)
                                // a_c.D = mu1/(rho1*Sc); //Diffusion coefficient
  a_c.D = 2. * 1e-5;             // Diffusion coefficient
  a_c.mol_mass = 0.032;         // molar mass
  a_c.He = 34.375;              // Henry's coeff
  a_c.cd = rho2 / a_c.mol_mass; // uniform concentration inside the bubble
  a_c.diff = true;              // we turn on diffusion of the tracer

#if relax_tt
  tt=0.001*TEND;
#else
  tt = 0.0; // we switch off the volume change (the bubble must have a constant radius)
#endif
  // TEND = TEND + tt;

  // solution setup
  TOLERANCE = 1e-7;
  DT = 1e-10;
  CFL = 0.25;
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
  No slip boundary condition. */
//u.t[left] = dirichlet(0.);
//u.n[left] = dirichlet(0.);

/**
 Navier slip boundary condition. */
u.t[left] = robin(0,Ls);
u.n[left] = dirichlet(0);

// Bcs, allow outflow at right and top
u.n[right] = neumann(0.); // zeroGradient of normal velocity
p[right] = dirichlet(0.);                       // zeroPressre at cell center
pf[right] = dirichlet(0.);                      // zeroPressure at cell face (?)
a_c[right] = dirichlet(zeta * a_c.cd / a_c.He); // concentration

u.n[top] = neumann(0.);
a_c[top] = dirichlet(zeta * a_c.cd / a_c.He); // concentration

/**
  We initialize the bubble with a radius $r=Rb$ and set the initial concentration in the liquid domain (according to the saturation ratio $\zeta$).
*/


event init(t = 0)
{
   if (!restore (file = "restart")) {
  double Hc = Rb * cos(theta_int / 180. * pi);
  refine(sq(x - Hc) + sq(y) < sq(1.25 * Rb) && level < LEVEL); // under-saturation
  fraction(f, sq(x - Hc) + sq(y) - sq(1.0 * Rb));
  // Init tr_c
  foreach ()
    a_c[] = f[] * zeta * a_c.cd / a_c.He; // bulk concentration
}
}

/** --------Laplace pressure inside the bubble----------------*/
double Temp = 293.15, Ru_cst = 8.314;
event change_saturation_concentration(i++)
{
#if laplace_term
  double vol = 0.;
  double R_bi = 0.;


  /*for the bubble on the surface*/
  scalar posy[], posx[], posxcl[]; // interface position
  double posxcl_1, posyheight_1;
  position(f, posx, {1, 0});   // x-coordinate of interface
  position(f, posy, {0, 1});   // y-coordinate of interface
  position(f, posxcl, {0, 1}); // position of contact line

  foreach (reduction(+
                     : vol))
  {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6)
    {
      posx[] = nodata;
      posy[] = nodata;
    }

    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || x > WIDTH / pow(2, LEVEL))
    {
      posxcl[] = nodata; // If vertical distance too far from surface then posycl nodata
    }

    vol += (1. - f[]) * dv();
  }
  posyheight_1 = statsf(posx).max;                                                                         // store H. position
  posxcl_1 = statsf(posxcl).max;                                                                           // store C.L. position
  
  R_bi = (pow(posxcl_1, 2) + pow(posyheight_1, 2))/2./posyheight_1;// equivalent radius

  a_c.cdi = a_c.cd + 2 * f.sigma / R_bi / Temp / Ru_cst; // interface BC at each time step
#else
  a_c.cdi = a_c.cd;
#endif
}

event acceleration(i++)
{
  face vector av = a;
  foreach_face(x)
       av.x[] -= 9.8; //gravity closed here

}

/*
event stability (i++) {
      DT = t < 0.95*tt ? 1e-8 : DT; //this helps for the stability
}*/


//1.0 define a function to calculate apparent contact angle.
double equivalent_contact_angle(double CLine1, double Height1)
{
  double x0 = 0., x1 = pi;
  while (x1 - x0 > 1e-4)
  {
    double x = (x1 + x0) / 2.;
    double f = CLine1 + CLine1*cos(x) - Height1*sin(x);
    //double f = Radius1 + Radius1 * cos(x) - Height1;
    if (f > 0.)
      x0 = x;
    else
      x1 = x;
  }
  return (x0 + x1) / 2.;
}

#if ca
//2.0 define a functiom to calculate local contact angle.
double local_contact_angle(double CLBC)
{
 Point point = locate (0., CLBC);
 //double dx = WIDTH / pow(2, LEVEL);
 //double ca = atan (1./(h.x[0,-1] - h.x[]));
 double CA_1 = atan (Delta/(h.y[-1,0] - h.y[0,0]));

return CA_1;
}
#endif


double posxcl1; // initialize position of contact line

// event log_simulation (i++)
event log_simulation(i = i +10)
// event log_simulation(t = 0; t <= TEND; t = t + 1e-6)
{
  double vol1 = 0.;
  scalar posy[], posx[], posxcl[]; // interface position
  position(f, posy, {0, 1});       // y-coordinate of interface
  position(f, posx, {1, 0});       // x-coordinate of interface
  position(f, posxcl, {0, 1});     // position of contact line

  double posxcl2 = posxcl1; // store C.L. position of last iteration
  t2 = t1;

  foreach (reduction(+
                     : vol1))
  {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6)
    {
      posx[] = nodata;
      posy[] = nodata;
    }

    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || x > WIDTH / pow(2, LEVEL))
    {
      posxcl[] = nodata; // If vertical distance too far from surface then posycl nodata
    }

    vol1 += (1. - f[]) * dv();
  }

  posxcl1 = statsf(posxcl).max; // store C.L. position
  t1 = t;

  double U_CL = (posxcl1 - posxcl2) / (t1 - t2); // calculate moving velocity of Contact line
  double ContactAngle = equivalent_contact_angle(posxcl1, statsf(posx).max)* 180. / pi;


#if dca
  fprintf(stderr, " %g %g %g %g %g %g %g %g %g %g %g  \n", t, vol1, statsf(posx).max, statsf(posy).max, posxcl1, U_CL, ucl2, Ca2, ContactAngle, theta_var, theta_local);
#else
  fprintf(stderr, "%g %g %g %g %g %g %g %g %g  \n", t, t - tt, vol1, statsf(posx).max, statsf(posy).max, posxcl1, U_CL, ContactAngle, theta0);
#endif
}

event interface(t = 0; t <= TEND; t = t +5e-7)
//event interface(i = i + 100)
{
  // output_facets (f, stderr);

  // The Following Code is important for parallel calculation!
  // Otherwise only get interface coordinate calculated by one of the cores
  char names[80];
  sprintf(names, "interface%d", pid());
  FILE *fp = fopen(names, "w");
  output_facets(f, fp);
  fclose(fp);

  char command[80];
  sprintf(command, "LC_ALL=C cat interfa* > facets/facets-%.7f", t);
  system(command);
}

event adapt(i++)
{
  double uemax = 1e-2;
  adapt_wavelet_leave_interface({a_c, u.x, u.y}, {f}, (double[]){1e-3, uemax, uemax, 1e-3}, LEVEL, basic_LEV, 1);
}

event snapshot1(t = 0; t <= TEND; t = t + 5e-7)
//event snapshot1(i = i + 100)
{
  char name[80];
  sprintf(name, "dump/dump-%.7f", t);
  dump(file = name);
}

event snapshot2(t = 0; t <= TEND; t = t + 5e-7)
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

