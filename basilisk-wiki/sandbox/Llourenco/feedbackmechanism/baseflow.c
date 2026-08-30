/**

# Baseflow at *Re=150*

This simulation establishes the baseflow definition of our system, acting as a stepping stone to launch the full simulations required to study the feedback mechanism


*/
#include "navier-stokes/centered.h"
#include "utils.h"
#include "view.h"

bid inner;

double Reynolds = 150; 
int maxlevel = 8;
double b = 1.; // jet diameter
double U0 = 1.; // for top hat
double W=10.; // edge-distance

scalar omega[];
scalar un[];
face vector muv[];
double cf=0.1; 
double du;
double Tend = 1000;
/**You can pass command-line arguments to the program to easily run parametric studies and test different configurations.*/

int main(int argc, char * argv[])
{
  if (argc > 4) {
    maxlevel = atof (argv[1]);
    Reynolds = atof (argv[2]);
    cf = atof (argv[3]);
    Tend = atof (argv[4]);
  }
  /**The interface includes interactive display controls to dynamically monitor and adjust key parameters during the simulation*/
  display_control(Reynolds,0,1000);
  display_control(cf,0,1);
  display_control(maxlevel,0,15);
  
/**The domain size is set to 32 and is vertically centered at the inlet. */   
  L0 = 32; 
  origin (0., -L0/2.);
  N = 64;
  mu = muv;
  init_grid (N);
  run();
}

event init (i = 0)
{
  if (!restore ("restart")) {
    //    fprintf (stderr, "Failed to restore\n");
    //   return 1;
  }
  /**To model the flat plate, the mesh is locally refined down to level 9 near the centerline, where a geometric mask is applied to define the solid boundary.*/
  refine(level < 9  && (fabs(y) < Delta));
  mask (fabs(y) <= L0/(1 << 9) ? inner : none);
/**Saves the current horizontal velocity field to allow for convergence tracking.*/
  foreach()
    un[] = u.x[];
}

/**We set a constant viscosity based on the jet width and the reference velocity*/
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds;
  /**Calculates the vorticity field and applies a localized mesh refinement near the nozzle exit to properly resolve the primary shear layer.*/
  vorticity (u, omega);
  refine (x < W/4. && fabs(y) < 1. && level < maxlevel);
}

/**We set the boundary conditions to form a jet. Outside of the jet, a coflow confines the flow, allowing it to encounter the flat plate properly.*/
u.n[left]  = dirichlet(fabs(y) > b/2. ? cf : 1.);
u.t[left]  = dirichlet(0.) ;

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);

u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);

/**Sets boundary conditions for the centerline and the plate: applies a free-slip/convective condition before the edge  and a strict no-slip condition on the solid plate.*/
u.n[inner] = x < W ? neumann(0.) : dirichlet(0.);
u.t[inner] = x < W ? dirichlet(1.) : dirichlet(0.);
/**
We monitor the iteration count, the time step, and the velocity change (du) to track the convergence of the simulation*/

event logfile (i++, t<=Tend)
{
  du = change (u.x, un);
  foreach() 
    un[]  = u.x[]; // allows to track convergence
  
  fprintf (stderr, "%d %g %g %g\n",
	   i, t, dt, du);
//  if (du < 5e-5)
//    {
//      dump();
//      return 1;
//    }
}

/**Dynamically adapts the mesh based on velocity errors to capture flow structures efficiently.*/
event adapt (i++) 
  adapt_wavelet ({u}, (double[]){1e-3,1e-3}, maxlevel, 4);

/**Saves the velocity components and vorticity to a .dat file and creates a simulation checkpoint*/
event dump(t=end)
{
  clear();
  view (width = 800, height = 400, fov = 15, tx = -0.3);
  squares ("u.x", linear = true);
  save ("velocity_x.png");
  
  FILE *fp;
  fp = fopen ("AeroField.dat","w");
  output_field({u.x,u.y,omega},fp,1251,box={{0,-W/2.},{1.25*W,W/2.}},linear=true);
  fclose(fp);
  dump();
}



event movie (t += 2; t <= Tend) {
  scalar m[];
  foreach()
    m[] = (fabs(y) <= L0/(1 << 9) && x > W) ? 0 : 1;
  
  output_ppm (omega, file="vorticity.mp4", n=512,
              box={{0,-16.0},{L0,16.0}},
              min=-3, max=3, linear=true, mask=m);
  output_ppm (u.x, file="ux.mp4", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-0.5, max=1.5, linear=true, mask=m);
}


/**
# Visualisations


![Animation of the vorticity field](baseflow/vorticity.mp4)(loop)

![Animation of the velocity field (u.x)](baseflow/ux.mp4)(loop)
*/

