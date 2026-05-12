// 3D dynamci wetting drop/bubble
//#include "grid/multigrid3D.h"
//#include "grid/quadtree.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
// htg head file
#include "output_htg.h"
#include "solver/contact_line_velocity_3d_interp.h"


#define F1_ERR 1e-3 //tolerance on the interfacial cells
#define WIDTH 5.0e-3
#define MAX_LEVEL 7
#define MIN_LEVEL 5
/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential components. */

double TEND = 0.03;
//# define Ls 10e-6 Ls too small
# define Rb 1.14e-3

// Gas-air
#define rho_g 1.25  // density 
#define mu_g 1.82e-5 // viscosity 
// Liquid-H2O
#define rho_l 1000. // density 
#define mu_l 1.0e-3 // viscosity 
// surface tension
#define f_sigma_ST 0.072


vector h[];
double theta_int = 0.;
double theta0 = 90.; // equilibrium contact angle
double theta_adv =  109.5; // max adv CA
double theta_rec = 71.8; // min rec CA
double Ratio = 8.88e6; // Ratio---HT model: Length ratio: L/l_m;  
//double lambda_mkt = 1.62e-9, kappa_mkt = 6e3; // lambda_mkt---MKT: average distance;kappa_mkt---MKT: frequency;
double xi_mkt = 0.002; // xi_mkt---MKT model: xi_mkt = kB_cst*Temp/kappa/lambda^3
double Cst01 = 0.9e4; // Cst01---Pinning slope coefficient ;

scalar thetai[], ucl2[], Ca2[];

event CA_bc(i++){
  // Compute ucl2[] over the domain
  contact_line_velocity_and_normal_vector(f, u, ucl2);
  // Synchronize ucl2[] across processes
  boundary({ucl2});
 
  foreach_boundary(bottom){
    //thetai[] = theta0;
    if (f[] > F1_ERR && f[] < 1. - F1_ERR) {

  
    Ca2[]=mu_l*ucl2[]/f_sigma_ST;

    //004-Combined HT and MKT model with Pinning;
    double C_pin;
    if (Ca2[]>=0.0)
      C_pin = f_sigma_ST*(cos(theta0 * pi / 180.)-cos(theta_adv * pi / 180.));
    else
      C_pin = f_sigma_ST*(cos(theta_rec * pi / 180.)-cos(theta0 * pi / 180.));
      
    thetai[] = acos(cos(theta0*pi/180.) - xi_mkt*Ca2[]/mu_l - C_pin*tanh(Cst01*Ca2[])/f_sigma_ST )*180./pi;
    thetai[] = pow(thetai[] / 180. * pi, 3) + 9.0 * Ca2[] * log(Ratio);
    thetai[] = pow(thetai[], 1.0 / 3) * 180. / pi;

    if (thetai[]>114.)
      thetai[]=114.;

    if (thetai[]<52.)
      thetai[]=52.;

    }
    
  }
  // Synchronize thetai[] across processes if needed
  boundary({thetai});
  boundary({Ca2});
  
  h.t[bottom] = contact_angle (thetai[] *pi/180.);
  h.r[bottom] = contact_angle (thetai[] *pi/180.);

  //h.t[bottom] = contact_angle (theta0 *pi/180.);
  //h.r[bottom] = contact_angle (theta0 *pi/180.);
}


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
  G.y = -9.8;
  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */
  
  f.sigma = f_sigma_ST;

  //N = 1 << MAXLEVEL;
  //for (theta0 = 30; theta0 <= 150; theta0 += 30)
  DT=2.5e-6;

  TOLERANCE = 1e-4;
  CFL=0.25;
  NITERMAX = 5000;
  run();
}


/**
The initial drop is a quarter of a sphere. */
event init (t = 0)
{
  if (!restore (file = "restart")) {
  double Hc = 1.0*Rb * cos(theta_int / 180. * pi);
  fraction (f,  -(sq(x) + sq(y-Hc) + sq(z) - sq(Rb)));
  refine( -(sq(x) + sq(y-Hc) + sq(z) - sq(1.2*Rb)) && level < MAX_LEVEL); // under-saturation
  
  // initial spreading velocity
  foreach(){
         u.y[] = -1.0*f[];   
    }
  
  }
}

/*
event acceleration(i++)
{
  face vector av = a;
  foreach_face(y)
    //if (f[]>0.)
      av.y[] -= 9.8; //gravity closed here
      //av.x[] -= 0.0;
}
*/

/* in htg file the domain is like this
           y     top
           -------------------------
          /|                      /|
         / |                     / |
        /  |            back    /  |
       -------------------------   |
left  |    |                   |   |  right
      |    |                   |   |
      |    |  front            |   |
      |  O --------------------|---- x
      |   /                    |  /
      |  /                     | /
      | /                      |/
      --------------------------
     z           bottom
*/
/* No slip boundary condition. */
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.); 

u.n[top] = neumann(0.); // zeroGradient of normal velocity
u.n[right] = neumann(0.); // zeroGradient of normal velocity
u.n[front] = neumann(0.); // zeroGradient of normal velocity



//p[top] = dirichlet(0.);                     // zeroPressure at cell face (?)
//pf[top] = dirichlet(0.);     
//p[right] = neumann(0.); 
//pf[right] = neumann(0.);

/**
We log statistics on the maximum velocity, curvature and volume. If
the standard deviation of curvature falls below $10^{-2}$, we assume
that the steady shape is reached and we stop the calculation. */

// adaptive mesh refinement
//#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-3, 1e-2,1e-2,1e-2},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
//#endif






event interface(t = 0; t <= TEND; t = t +1e-4)
//event interface(i += 100)
{

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

event snapshot1(t = 0; t <= TEND; t = t + 1e-4)
//event snapshot1(i += 100)
{
  char name[80];
  sprintf(name, "dump/dump-%.7f", t);
  dump(file = name);
}

event snapshot2(t = 0; t <= TEND; t = t + 1e-4)
//event snapshot2(i =i+ 100)
{
  char path[] = "htg"; // no slash at the end!!
  char prefix[80];

  sprintf(prefix, "data-%.7f", t);
  output_htg((scalar *){f,p,u.x,u.y,u.z}, (vector *){uf}, path, prefix, i, t);
}


event stop_simulation(t = TEND)
{
  return 1;
}

/*
double iL = HUGE;

event stoprun (i++){
if (i < iL)
exit(1);
iL = i;
}
*/

