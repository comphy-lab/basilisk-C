/**
# Stefan problem

Here we model the Stefan problem for a planar interface which separates the gas and liquid regions. We assume a null initial concentration of gas in the liquid domain; the diffusion of gas drives the displacement of the interface.

The analytical solution of the interface displacement is:

$$
l(t) = \frac{2}{He} \sqrt{\frac{t\, D}{\pi}}
$$

and the concentration field in the liquid domain:

$$
c(y,t) = c_\Sigma \left(1 - erf \left( \frac{y - y_\Sigma(t)}{2\sqrt{D\, t}} \right) \right)
$$

A similar example is discussed in Section 4.3 of [Gennari et al.,2022](#gennari2022).
*/

#define F_ERR 1e-10

#include "grid/multigrid.h"
#define PHASE_CHANGE 1
#include "../navier-stokes/centered.h"
#include "../two-phase.h"
#include "../diffusion.h"
#include "../phase_change_pure_species.h"
#include "view.h"

/** 
Two-phase system properties
*/

#define RHOR 1000. //density ratio
#define MUR 100. //viscosity ratio
#define Sc 10. //Schmidt number

#define WIDTH 10
#define TEND 175

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

  a_c.inverse = false;
  a_c.D = mu1/(rho1*Sc); //Diffusion coefficient
  a_c.mol_mass = 1./RHOR; //molar mass
  a_c.cd = 1.; //uniform concentration inside the bubble
  a_c.He = 1.2; //Henry's coeff
  a_c.diff = true; //we turn on diffusion of the tracer
  
  tt = 0.; //we switch on volume change at t=0

  TOLERANCE = 1e-4;
  DT = 1e-3;
  run();
}

/**
We need outflow boundary conditions to let the liquid enter the domain as the gaseous region dissolves. 
*/

//Bcs
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
a_c[top] = dirichlet(0.);

/**
We initialize the two-phase domain.
*/

event init (t = 0) {
  fraction (f, y - 4.7);
}

event stability (i++) {
  DT = t < 1. ? 1e-4 : 1e-2; //This helps the stability
}

/**
We monitor the volume of the gaseous region.
*/

event log_simulation (i++) {
  double vol = 0.;
  foreach(reduction(+:vol))
    vol += (1. - f[])*dv();

  fprintf (stderr, "%g %g\n", t, vol);
  fflush (stderr);
}

/**
We make a movie of the phase-average concentration around the interface.
*/

event movie (t++) {
  scalar c_c[]; //phase-average concentration
  foreach()
    c_c[] = f[] > F_ERR ? a_c[]/f[] : a_c.cd;
  view (quat = {0.002, 0.005, 0.000, 1.000},
	fov = 35, near = 0.01, far = 1000,
	tx = -0.5, ty = -0.5, tz = -2.,
	width = 600, height = 600);
  clear();
  box();
  draw_vof ("f");
  squares("c_c", min=0., max=0.8, linear=true);
  save ("movie.mp4");
}

event stop_simulation (t = TEND) {
  return 1;
}

/**
## References
~~~bib
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