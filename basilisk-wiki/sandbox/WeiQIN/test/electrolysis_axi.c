/**
# Start with a growing bubble in a super-saturated solution

Here we model a spherical bubble growing in a super-saturated solution. The saturation ratio $\zeta$ compares the bulk liquid concentration and the interfacial concentration (for a saturated interface):

$$
\zeta = \frac{c_{bulk}}{c_\Sigma}
$$

When $\zeta > 1$, the solution is defined super-saturated and the excess of species (with respect to the saturated concentration) is redistributed through a diffusive process across the interface; the volume of the bubble increases accordingly.

This problem was solved analytically by Scriven in [Scriven,1959](#scriven1959) and the time-evolving radius reads:

$$
R(t) = 2 \beta \sqrt{D\, t}
$$

where $\beta$ is a non-dimensional factor that depends on the gas and liquid densities, liquid bulk concentration and interfacial concentration.

Here we assume a $H_2$ - water system with an increased diffusivity of gas. The saturation ratio is set to $\zeta = 1$ for the whole domain. The volume change is turned on at $t=0.02s$ from the previous experiment, when the predicted bubble radius (according to Scriven's formula) reaches the size of the bubble initialization in the present simulation.

*/

#define F_ERR 1e-10

#include "grid/quadtree.h"
#include "axi.h"
#define PHASE_CHANGE 1
#include "ggennari/phase_change/navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "ggennari/phase_change/two-phase.h"
#include "tension.h"
#include "ggennari/phase_changediffusion.h"
#include "ggennari/phase_changephase_change_pure_species.h"
#include "ggennari/adapt_wavelet_leave_interface.h"

/**
Two-phase system properties
*/
#define diameter 1. //numerical diameter
#define WIDTH 25.*diameter // computation domain
#define D_b 0.0000127 //m diameter of the initial bubble
#define D_e 0.000127 //m diameter of the electrode
#define Mu_c 8.32e-4 //Pa.s
#define Mu_d 8.96e-6 //Pa.s
#define rho_c 996. //kg/m**3
#define rho_d 0.8   //kg/m**3
#define sigma12 75.e-7 //surface tension water/air
#define G 9.8 //m**2/s
#define M 0.020 //kg/mol
#define C_r 44 //hydrogen concentration mol/m**3
#define D_f 7.38e-9 //diffusion coefficient m**2/s
#define I 100  //current  A/m^2
#define F 96485.3 // Faraday’s constant As/mol

#define RHOR (rho_c/rho_d) //density ratio
#define MUR (Mu_c/Mu_d) //viscosity ratio

#define sat 1. //saturation ratio
#define J (I/(2*F)) //Faraday’s law


#define Sc (Mu_c/(rho_c*D_f)) //Schmidt number
#define ratio (diameter/D_b) //convert all geometrical lengths in numerical value  // the physical domain of computation is defined by WIDTH/ratio


#define SIGMA (sigma12*rho_c*D_b)/sq(Mu_c)
#define G_c (J*D_b/(D_f*C_r))


#define TEND 1000.

int LEVEL = 9;
/**
The concentration field is stored in the scalar *a_c*.
*/

scalar a_c[];
scalar * species = {a_c};
scalar m_a[];
scalar * m_list = {m_a};

int main()
{
  size (WIDTH);
  init_grid (64);

  rho1 = 1.;
  rho2 = rho1/RHOR;
  mu1 = 1.;
  mu2 = mu1/MUR;

  f.tracers = species;
  f.sigma =SIGMA;
  //f.height = h;
 
 
  a_c.inverse = false;
  a_c.D = mu1/(rho1*Sc); //Diffusion coefficient
  a_c.mol_mass = 1./RHOR; //molar mass
  a_c.cd = 1.; //uniform concentration inside the bubble
  a_c.He = 53.3; //Henry's coeff
  a_c.diff = true; //we turn on diffusion of the tracer
 
  tt = 103.6; //we switch on volume change when the bubble reaches the size of its initialization

  TOLERANCE = 1e-5;
  DT = 1e-1;
  run();
}

/**
We need outflow boundary conditions to let the liquid exit the domain as the bubble grows.
*/

//Bcs
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
a_c[left] = neumann(y < D_e*ratio/2.0 ? (G_c*f[]):0);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

/*

*/
/**
We initialize the bubble with a radius $r=0.5$ (numerical length) and set the initial concentration in the liquid domain (according to the saturation ratio $\zeta$).
*/

event init (t = 0) {
  refine (sq(x) + sq(y) < sq(0.75*diameter) && level < LEVEL);
  fraction (f, sq(x) + sq(y) - sq(0.5*diameter));

  //Init tr_c
  foreach()
    a_c[] = f[]*sat*a_c.cd/a_c.He;
}

event stability (i++) {
  DT = t < 0.05 ? 1e-4 : 1e-3; //this helps for theS stability
}

/**
We monitor the volume of the bubble.
*/

event log_simulation (i++) {
  double vol = 0.;
  foreach(reduction(+:vol))
    vol += (1. - f[])*dv();

  fprintf (stderr, "%g %g\n", t, vol*2.*pi);
  fflush (stderr);
}

event adapt (i++)
{
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({a_c,u.x,u.y}, {f}, (double[]){1e-3,uemax,uemax,1e-3}, LEVEL, 5, 1);
}

/**
We make a movie of the phase-average concentration around the interface.
*/

event viewing(t += 10.) {
    scalar c_c[]; //phase-average concentration
    foreach()
    c_c[] = f[] > F_ERR ? a_c[]/f[] : 1;
    view(fov=5, tx=0,ty=-0.131,quat ={0,0,-0.707,0.707},
     width = 1920, height = 1080);
    clear();
    draw_vof (c = "f", lw = 2);
    squares("c_c", min=0., max=0.2, linear=true);
    mirror({0,1}){
    draw_vof (c = "f", lw = 2);
    cells ();
    }
    char png_filename[100];
    sprintf (png_filename,"Z2Interface/time-%4.2f.png",timeload);
    save(png_filename);
    save ("movie.mp4");
}


event stop_simulation (t = TEND) {
  return 1;
}


event snapshots(t += 10) {
  char name[100];
 
  sprintf (name, "dump-%4.2f", t);
  dump (name);
}

/**
## References
~~~bib
@article{scriven1959,
   author = {Scriven, L. E.},
   title = {On the dynamics of phase growth},
   journal = {Chemical Engineering Science},
   volume = {10},
   number = {1},
   pages = {1-13},
   ISSN = {0009-2509},
   DOI = {https://doi.org/10.1016/0009-2509(59)80019-1},
   url = {https://www.sciencedirect.com/science/article/pii/0009250959800191},
   year = {1959},
   type = {Journal Article}
}

@article{gennari2022,
title = {A phase-change model for diffusion-driven mass transfer problems in incompressible two-phase flows},
journal = {Chemical Engineering Science},
volume = {259},
pages = {117791},
year = {2022},
issn = {0009-2509},
doi = {https://doi.org/10.1016/j.ces.2022.117791},
url = {https://www.sciencedirect.com/science/article/pii/S000925092200375X},
author = {Gabriele Gennari and Richard Jefferson-Loveday and Stephen J. Pickering}
}
~~~
*/
