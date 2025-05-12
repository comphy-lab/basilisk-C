/**
# Elasto-viscoplastic 2D droplet in a Couette Newtonian shear flow

We model the shear deformation of an elasto-viscoplastic (Saramito's model (2007)) droplet in a Newtonian matrix. After a time, the imposed deformation stops, and the droplet is left to retract.

## Code
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "tension.h"
#include "view.h"

#define Ca 0.6     // Capillary number
#define Re 0.3     // Reynold number
#define We (Ca*Re) // Weber number
#define MUr 1.     // ratio of outer(matrix) to inner(drop) viscosity
#define M 1.       // ratio of outer to inner density
#define Deb 0.4    // Deborah number
#define Beta 0.5   // ratio of the solvent visc. to the total viscoelastic visc.

/**
We set a maximum level of 8. For higher levels the simulation breaks.*/

int MAXLEVEL = 8;
int MINLEVEL = 6;

double tEnd = 6.;

scalar mupd[], lam[], vel[], yielded[];

/**
The top and bottom boundary conditions are those of a Couette flow. After $t = 3.0 s$ the shear deformation stops and the droplet retracts */

u.t[top] = dirichlet (t <= 3. ? y : 0);
u.t[bottom] = dirichlet (t <= 3. ? y : 0);

/** 
The domain is a 16x16 box which will be masked later to a become a 16x8 size box. */

int main() {
  L0 = 16.;
  origin (-L0/2, -L0/4.);
  periodic (right);
  DT = 1e-1;

  /**
  We set the viscosities, densities, surface tension and visco-elastic
  parameters. $\mu_p$ and $\lambda$ are defined in the properties event. */
  
  mu1 = MUr*Beta/Re;
  mu2 = 1./Re;
  rho1 = M;
  rho2 = 1.;
  f.sigma = 1./We;
  lambda = lam;
  mup = mupd;

  init_grid (1 << MAXLEVEL);
  run();
}


/**
The droplet radius is equal to 1.
*/
event init (i = 0) 
{
  mask (y > 4 ? top : none);

  fraction (f, 1. - (sq(x) + sq(y)));
  foreach()
    u.x[] = y;
  boundary ((scalar *){u});
}



//tensor shear[];
//scalar shear_norm[], eta[];
scalar tau_norm[], deviatoric[];
scalar fa[];

event properties (i++)
{
  double tau_y = 1.;  // Yield stress
  double mu_p = 1.;   // Plastic viscosity
  double epsilon = 0.01;  // Regularization parameter

  foreach()
  {
/*    shear.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

    shear_norm[] = sqrt(0.5)*sqrt(sq(shear.x.x[])+sq(shear.x.y[])+sq(shear.x.y[])+sq(shear.y.y[]));
    eta[] = (tau_y/(shear_norm[] + epsilon) + mu_p);*/

    tau_norm[] = sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]));  // Norm of the polymeric stress tensor

    if (tau_norm[] <= tau_y + epsilon)
      deviatoric[] = tau_y + epsilon;
    else
      deviatoric[] = tau_norm[];

    fa[] =  (4.*f[] + 2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) + f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    mupd[] = deviatoric[]/(deviatoric[] - tau_y) * fa[];
    lam[] = Deb * mupd[] * fa[];

    vel[] = pow(sq(u.x[])+sq(u.y[]),0.5);
  }
  boundary ({tau_norm, deviatoric, fa, mupd, lam, vel});

 foreach()
 {
    if(tau_norm[] > tau_y && f[] > 0.01)
    {
      yielded[] = 1.;
    }
    else
    {
      yielded[] = 0.; 
    }
  }
  boundary ({yielded});
}

/**
The mesh is adapted according to the errors on volume fraction,
velocity, yield surface and polymeric viscosity. */

event adapt (i++) {
  adapt_wavelet ({f, u, yielded, mupd}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);
}


/*
event uprofile (t += 0.2) 
{
  char namei1[80];
  sprintf (namei1, "interface-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f, fp1);
  fclose (fp1);

  char field2[80];
  sprintf (field2, "fields-output2-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  output_field ((scalar *){u, mupd, lam, tau_p, yielded}, fld2, n = 500, linear = true,  box = {{-2., -2.},{2.,2.}});
  fclose (fld2);
}*/


/** We generate some videos
*/

event movie_vel (t += 0.01)
{
  char movie_vel[80];
  sprintf (movie_vel, "movie_vel.mp4");
  clear(); 
  view(fov = 6., tx = 0.0, ty = 0.0);
  draw_vof ("f");
  squares("vel", linear = true);
  save(movie_vel);
}


event movie_mu (t += 0.01)
{
  char movie_vel[80];
  sprintf (movie_vel, "movie_mu.mp4");
  clear(); 
  view(fov = 6., tx = 0.0, ty = 0.0);
  draw_vof ("f");
  squares("mupd", linear = true);
  save(movie_vel);
}


event movie_tau_p (t += 0.01)
{
  scalar tau_pp[];
  foreach()
  {
      tau_pp[] = pow(0.5*(sq(tau_p.x.x[])+sq(tau_p.y.y[])+sq(tau_p.x.y[])+sq(tau_p.y.x[])),0.5);
  }
  boundary({tau_pp});
  char movie_vel[80];
  sprintf (movie_vel, "movie_tau_pp.mp4");
  clear(); 
  view(fov = 6., tx = 0.0, ty = 0.0);
  draw_vof ("f");
  squares("tau_pp", linear = true);
  save(movie_vel);
}


event movie_yield (t += 0.01)
{
  char movie_vel[80];
  sprintf (movie_vel, "movie_yield.mp4");
  clear(); 
  view(fov = 6.,tx = 0.0, ty = 0.00);
  draw_vof ("f");
  squares("yielded");
  save(movie_vel);
}

event end (t = tEnd) {}


/**
## Results

### Velocity field
![Velocity field](shear_droplet/movie_vel.mp4)

### Polymeric stress field
![Polymeric stress](shear_droplet/movie_tau_pp.mp4)

### Polymeric viscosity ($\mu_p$) field
![Pressure field](shear_droplet/movie_mu.mp4)

### Yielded region
![Yielded region](shear_droplet/movie_yield.mp4)

*/