/**
# Energy evolution of a pressure-forced linear wave

The surface pressure contribution is proportional to the local surface
slope, as in the original forcing calculation:
$$
p_s = -\hat p\,\partial_x\eta, \qquad
\hat p = 2\nu\sqrt{gk}.
$$
The pressure gradient is added to the barotropic acceleration. We retain
the original linear wave, depth, steepness and viscosity, and use the
same grid and duration as [viscous.c](/src/test/viscous.c).

As shown below, adding the forcing term maintains a nearly constant total energy. The slight increase in energy is due to the fact that the numerical method slightly under-estimates the energy loss compared to the analytical predicition by Joseph & Wang, as clearly shown below.

~~~gnuplot Evolution of the normalised energy.
set xlabel 't/T_0'
set ylabel 'Normalised energy E/E_0'
k = 2.*pi
lambda = 1.
g = 1.
T0 = k*lambda/sqrt(g*k)
nu = 1./20000.
set xrange [0:10]
set key bottom left
plot 'log' u ($1/T0):4 every 30 w p t 'Numerical (forcing)', \
     '/home/basilisk-src/basilisk/src/test/viscous/log' u ($1/T0):4 every 30 w p t 'Numerical (no forcing)', \
     exp(-2.*nu*k**2*x*T0) w l t 'Joseph and Wang decay', \
     1 t 'constant'
~~~

The macros below add a surface-pressure gradient to the barotropic
acceleration, using $P/\rho=\hat p\,\partial_x\eta$. Following the surface-pressure forcing from Fedorov & Melville
(1998), we retain the original coefficient
$\hat p=2\nu\sqrt{gk}$. In the linear deep-water limit this compensates
the $2\nu k^2$ energy decay discussed in [viscous.c](/src/test/viscous.c).
*/

double phat;
#define p_surface(eta, i) (- phat*(eta[i+1] - eta[i-1])/(2*Delta))
#define p_baro(eta, i) (- G*eta[i] + p_surface(eta, i))
#define a_baro(eta, i) (gmetric(i)*(p_baro (eta, i) - p_baro (eta, i - 1))/Delta)

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

double k_ = 2.*pi, h_ = 2.0, g_ = 1., ak = 0.01;
double RE = 20000.;

#define T0 (k_*L0/sqrt(g_*k_))

int main()
{
  origin (-L0/2.);
  periodic (right);
  G = g_;
  nu = 1./RE;

  N = 256;
  nl = 60;
  NITERMIN = 3;
  DT = 1e-2 [0,1];

  run();
}

/**
Initialise the original linear wave on geometrically distributed layers.
The surface tension in the original calculation was zero. */

event init (i = 0)
{
  geometric_beta (1./3., true);
  foreach() {
    zb[] = -h_;
    eta[] = ak/k_*cos(k_*x);
    double H = eta[] - zb[];
    foreach_layer()
      h[] = H/nl;
  }
  vertical_remapping (h, tracers);
  foreach() {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      u.x[] = ak*sqrt(g_/k_)*cos(k_*x)*exp(k_*z);
      w[] = ak*sqrt(g_/k_)*sin(k_*x)*exp(k_*z);
      z += h[]/2.;
    }
  }
}

/**
Set the pressure coefficient and add horizontal viscosity. */

event face_fields (i++) {
  phat = 2.*nu*sqrt(g_*k_);
}

event viscous_term (i++) {
  horizontal_diffusion ({u, w}, nu, dt);
}

event profiles (t += T0/4.; t <= 10.*T0) {
  foreach (serial) {
    double H = zb[];
    foreach_layer()
      H += h[];
    printf ("%g %g\n", x, H);
  }
  printf ("\n\n");
}

/**
Log time, kinetic energy, potential energy and normalised total energy
in stderr (the sandbox's `log` file). As in decay.c, normalise using the
first timestep to skip the initial relaxation. */

event logfile (i = 1; i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  double Ec = ke/2., Ep = g_*(gpe + L0*sq(h_)/2.);
  double E = Ec + Ep;
  static double E0 = 0.;
  if (i == 1)
    E0 = E;  
  fprintf (stderr, "%g %g %g %g\n", t, Ec, Ep, E/E0);
}

/**
## References

~~~bib
@article{fedorov1998,
  title={Nonlinear gravity--capillary waves with forcing and dissipation},
  author={Fedorov, Alexey V. and Melville, W. Kendall},
  journal={Journal of Fluid Mechanics},
  volume={354},
  pages={1--42},
  year={1998},
  doi={10.1017/S0022112097007453}
}
~~~
*/