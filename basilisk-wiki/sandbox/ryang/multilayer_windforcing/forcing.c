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
set xlabel 't/T_0' font 'Times New Roman,16'
set ylabel 'Normalised energy E/E_0' font 'Times New Roman,16'
set tics font 'Times New Roman,14'
set key bottom left font 'Times New Roman,14'
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

/**
We reuse the Stokes wave initialisation routine for simplicity. The
nonlinear Stokes terms do not matter since the amplitude is small,
although this could be checked by using a pure linear wave
solution. */   

#include "test/stokes.h"

/**
The wavenumber is $k=2\pi$, the amplitude/slope $ak$ is small enough
(decreasing it has little effect on the results). The Reynolds number is 
$\text{RE} = c\lambda/\nu$ with $c = \sqrt{g/k}$ the deep-water phase speed. 
Decreasing the viscosity (or increasing it) does not
affect the results qualitatively.

Note that we set the fluid depth to h_ = 2 (i.e. twice the wavelength)
to be in the deep water limit. This is important to match the
theory. Further increasing the depth does not change the results. */

double k_ = 2.*pi [-1,0];
double h_ = 2.0 [1,0];
double g_ = 1. [1,-2];
double ak = 0.01 [0,0];
double RE = 7979. [0,0]; // Re = c*lambda/nu; this gives nu ≈ 1./20000., matching viscous.c (/src/test/viscous.c).

#define T0 (k_*L0/sqrt(g_*k_))

int main()
{
  origin (-L0/2.);
  periodic (right);
  G = g_;
  double lambda_ = 2.*pi/k_;
  double c_ = sqrt(g_/k_);
  nu = c_*lambda_/RE;

  N = 256;
  nl = 60;
  NITERMIN = 3;
  DT = 1e-2 [0,1];

  run();
}

/**
Initialise the original linear wave on geometrically distributed layers. */

event init (i = 0)
{
  geometric_beta (1./3., true);
  
  foreach() {
    zb[] = - h_;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
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
in stderr (the sandbox's `log` file). */

event logfile (i++)
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
  if (i == 0)
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