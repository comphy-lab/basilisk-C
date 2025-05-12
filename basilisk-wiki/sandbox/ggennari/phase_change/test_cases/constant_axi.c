/**
# Axisymmetric bubble with constant mass transfer rate

We model an axisymmetric dissolving bubble with constant mass transfer rate. A similar test can be found in Section 4.2 of [Gennari et al.,2022](#gennari2022).

The analytical solution for the radius reads:
$$
R(t) = R(0) + \frac{\dot{m}}{\rho_d}
$$
*/

#define F_ERR 1e-10
#define MASS_TRANSFER_RATE (-1e-3)

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

#define RHOR 1000.
#define MUR 100.

#define WIDTH 5
#define TEND 1.

int LEVEL = 9;

/**
We add the transport of a soluble tracer. But this is not relevant in this case, since the mass transfer rate is fixed.
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
  a_c.mol_mass = 1./RHOR;
  a_c.diff = false;
  
  tt = 0.;

  TOLERANCE = 1e-4;
  DT = 1e-3;
  run();
}

/**
We need outflow boundary conditions to let the liquid enter the domain as the bubble dissolves. 
*/

//Bcs
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

event init (t = 0) {
  refine (sq(x) + sq(y) < sq(4) && level < LEVEL);
  fraction (f, sq(x) + sq(y) - sq(1.));
}

/**
We log the volume of the bubble
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
  adapt_wavelet_leave_interface ({u.x,u.y}, {f}, (double[]){uemax,uemax,1e-3}, LEVEL, 5, 1);
}

event movie (t += 0.01) {
  view (quat = {0.002, 0.005, 0.000, 1.000},
	fov = 30, near = 0.01, far = 1000,
	tx = -0.1, ty = -0.1, tz = -1.,
	width = 600, height = 600);
  clear();
  box();
  draw_vof ("f");
  squares("u.x", min=-1.5, max=1.5, linear=true);
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