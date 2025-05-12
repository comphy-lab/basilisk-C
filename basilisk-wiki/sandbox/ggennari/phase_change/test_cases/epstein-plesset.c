/**
# Dissolving bubble with steady shape

Here we model a spherical bubble with a constant radius $R$ that releases gas in the sourranding liquid. This is a pure diffusion problem:

$$
\frac{\partial c}{\partial t} = \nabla \cdot \left(D \nabla c \right)
$$

$$
c(r,0) = c_{bulk} \quad \text{for} \quad r > R
$$

$$
c(R,t) = c_\Sigma \quad \text{for} \quad t > 0
$$

$$
\lim_{r\to\infty} c(r,t) = c_{bulk} \quad \text{for} \quad t > 0
$$

where $c_\Sigma$ is the interfacial concentration (given by Henry's Law) and the bulk liquid concentration is assumed null, i.e. $c_{bulk}=0$.

This problem was solved analytically by [Epstein and Plesset,1950](#epstein1950):
$$
\frac{dR}{dt} = \frac{D\, M\, (C_{bulk} - c_\Sigma)}{\rho_d} \left[\frac{1}{R} + \frac{1}{\sqrt{\pi\, D\, t}} \right]
$$

where the evolution of the bubble radius is predicted based on the mass flux across the interface.

The concentration field in the liquid region ($r>R$) is (see [Crank,1975](#crank1975)):
$$
c(r,t) = c_{bulk} + (c_\Sigma - c_{bulk})\frac{R}{r}erfc\left(\frac{r-R}{2\sqrt{D\, t}}\right)
$$

This example is discussed in Section 4.4.2 of [Gennari et al.,2022](#gennari2022). 
*/

#define F_ERR 1e-10

#include "grid/quadtree.h"
#include "axi.h"
#define PHASE_CHANGE 1
#include "../navier-stokes/centered.h"
#include "../two-phase.h"
#include "../diffusion.h"
#include "../phase_change_pure_species.h"
#include "../../adapt_wavelet_leave_interface.h"

#define RHOR 554. //density ratio
#define MUR 72. //viscosity ratio
#define Sc 0.0526 //Schmidt number

#define WIDTH 10
#define TEND 0.01

int LEVEL = 10;
double bubble_mass = 0.;

/**
The concentration field is stored in the scalar *a_c*
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
  a_c.He = 5.; //Henry's coeff
  a_c.diff = true; //we turn on diffusion of the tracer
  
  tt = 100.; //we switch off the volume change (the bubble must have a constant radius)

  TOLERANCE = 1e-4;
  DT = 1e-5;
  run();
}

//Bcs
a_c[right] = dirichlet(0.);
a_c[top] = dirichlet(0.);

event init (t = 0) {
  refine (sq(x) + sq(y) < sq(0.75) && level < LEVEL);
  fraction (f, sq(x) + sq(y) - sq(0.5));

  //Init tracer
  foreach()
    a_c[] = 0.;
  
  bubble_mass = 2./3.*pi*pow(0.5, 3.)*rho2; //Initial bubble mass
}

/**
We monitor the mass of the bubble by computing at each time step the amount of gas released into the liquid. We also log the concentration field at three locations $r=0.6$, $r=0.8$ and $r=1$.
*/

event log_simulation (i++) {
  double vol = 0., dm = 0.;
  foreach(reduction(+:vol) reduction(+:dm)) {
    vol += (1. - f[])*dv();
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal(point, f), p;
      double alpha = plane_alpha (f[], n);
      double area = plane_area_center (n, alpha, &p);
      dm += m_a[]*area*pow(Delta, dimension - 1)*(y + p.y*Delta)*(2.*pi)*dt;
    }
  }
  bubble_mass += dm;
  coord c[3];
  c[0] = (coord){0.6*cos(M_PI/4.), 0.6*sin(M_PI/4.)};
  c[1] = (coord){0.8*cos(M_PI/4.), 0.8*sin(M_PI/4.)};
  c[2] = (coord){1.*cos(M_PI/4.), 1.*sin(M_PI/4.)};
  double V[3];
  interpolate_array ({a_c}, c, 3, V, true);
  fprintf (stderr, "%g %g %g %g %g %g\n",
	   t, vol*2.*pi, bubble_mass, 
	   V[0], V[1], V[2]);
  fflush (stderr);
}

event adapt (i++)
{
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({a_c,u.x,u.y}, {f}, (double[]){1e-3,uemax,uemax,1e-3}, LEVEL, 5, 1);
}

event stop_simulation (t = TEND) {
  return 1;
}

/**
## References

~~~bib
@article{epstein1950,
  title={On the stability of gas bubbles in liquid-gas solutions},
  author={Epstein, Paul S and Plesset, Milton S},
  journal={The Journal of Chemical Physics},
  volume={18},
  number={11},
  pages={1505--1509},
  year={1950},
  publisher={American Institute of Physics}
}

@book{crank1975,
   author = {Crank, J.},
   title = {The mathematics of diffusion / by J. Crank},
   publisher = {Clarendon Press},
   address = {Oxford [England]},
   series = {Oxford science publications.},
   ISBN = {0198533446
0198534116},
   year = {1975},
   type = {Book}
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