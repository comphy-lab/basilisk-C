/**
# Growing bubble in a super-saturated solution

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

Here we assume a $CO_2$ - water system with an increased diffusivity of gas. The saturation ratio is set to $\zeta = 2$ and the growth factor is found to be $\beta=0.421$. The volume change is turned on at $t=0.0186$, when the predicted bubble radius (according to Scriven's formula) reaches the size of the bubble initialization in the present simulation.

The details of this example can be found in Section 4.4.1 of [Gennari et al.,2022](#gennari2022)
*/

#define F_ERR 1e-10

#include "grid/quadtree.h"
#include "axi.h"
#define PHASE_CHANGE 1
#include "../navier-stokes/centered.h"
#include "../two-phase.h"
#include "tension.h"
#include "../diffusion.h"
#include "../phase_change_pure_species.h"
#include "../../adapt_wavelet_leave_interface.h"
#include "view.h"

/** 
Two-phase system properties
*/

#define RHOR 554. //density ratio
#define MUR 72. //viscosity ratio
#define Sc 0.0526 //Schmidt number
#define sat 2. //saturation ratio

#define WIDTH 25
#define TEND 1.

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
  f.sigma = 1.;

  a_c.inverse = false;
  a_c.D = mu1/(rho1*Sc); //Diffusion coefficient
  a_c.mol_mass = 1./RHOR; //molar mass
  a_c.cd = 1.; //uniform concentration inside the bubble
  a_c.He = 5.; //Henry's coeff
  a_c.diff = true; //we turn on diffusion of the tracer
  
  tt = 0.01871; //we switch on volume change when the bubble reaches the size of its initialization

  TOLERANCE = 1e-4;
  DT = 1e-3;
  run();
}

/**
We need outflow boundary conditions to let the liquid exit the domain as the bubble grows. 
*/

//Bcs
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
a_c[right] = dirichlet(sat*a_c.cd/a_c.He);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
a_c[top] = dirichlet(sat*a_c.cd/a_c.He);

/**
We initialize the bubble with a radius $r=0.5$ and set the initial concentration in the liquid domain (according to the saturation ratio $\zeta$).
*/

event init (t = 0) {
  refine (sq(x) + sq(y) < sq(0.75) && level < LEVEL);
  fraction (f, sq(x) + sq(y) - sq(0.5));

  //Init tr_c
  foreach()
    a_c[] = f[]*sat*a_c.cd/a_c.He;
}

event stability (i++) {
  DT = t < 0.05 ? 1e-4 : 1e-3; //this helps for the stability
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

event movie (t += 0.01) {
  scalar c_c[]; //phase-average concentration
  foreach()
    c_c[] = f[] > F_ERR ? a_c[]/f[] : a_c.cd;
  view (quat = {0.002, 0.005, 0.000, 1.000},
	fov = 30, near = 0.01, far = 1000,
	tx = -0.1, ty = -0.1, tz = -1.,
	width = 600, height = 600);
  clear();
  box();
  draw_vof ("f");
  squares("c_c", min=0.2, max=0.4, linear=true);
  save ("movie.mp4");
}

event stop_simulation (t = TEND) {
  return 1;
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