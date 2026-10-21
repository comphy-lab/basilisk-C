/**
# Viscous dissipation of a linear wave

Using the multilayer model, we try to reproduce the prediction of Lamb
(1932) that the total energy of a small amplitude (linear) water wave
should decrease like
$$
E(t) = E_0 \exp(-4 \nu k^2 t)
$$
with $\nu$ the kinematic viscosity and $k$ the wavenumber.

Joseph & Wang (2004) also showed that if a free-slip rather than a
free-surface stress-free condition at the top is applied, one gets
instead
$$
E(t) = E_0 \exp(-2 \nu k^2 t)
$$

We get the following figure, where the numerical results have been
checked for numerical convergence. The numerics agree neither with
Lamb (stronger dissipation), nor with (Joseph & Wang) (weaker
dissipation).

~~~gnuplot Evolution of the numerical and theoretical energies.
set xlabel 'Time'
set ylabel 'Normalised energy E/E_0'
k = 2.*pi
nu = 1./20000.
set xrange [0:25]
plot 'log' u 1:4 every 30 w p t 'Numerical', \
     exp(-2.*nu*k**2*x) w l t 'Joseph and Wang', \
     exp(-4.*nu*k**2*x) w l t 'Lamb'
~~~

## References

~~~bib
@article{joseph2004,
  title={The dissipation approximation and viscous potential flow},
  author={Joseph, Daniel D and Wang, Jing},
  journal={Journal of Fluid Mechanics},
  volume={505},
  pages={365--377},
  year={2004},
  publisher={Cambridge University Press}
}
~~~
*/

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
The wavenumber is $k=2\pi$, the amplitude/slope $ak$ is small. The
viscosity is $\nu=1./RE$. */

double k_ = 2.*pi, h_ = 0.5, g_ = 1., ak = 0.01;
double RE = 20000.;

/**
This is the wave period. */

#define T0  (k_*L0/sqrt(g_*k_))

int main()
{
  origin (-L0/2.);
  periodic (right);
  G = g_;
  nu = 1./RE;

  /**
  We adjust the numerical parameters for convergence. The vertical
  resolution matters, as do the convergence of the multigrid solver
  and timestep. */
  
  N = 256;
  nl = 60;
  NITERMIN = 3;
  DT = 1e-2 [0,1];
  
  run();
}

/**
We optimise layers using a geometric distribution to better resolve
the surface, and initialize the wave height and velocity
components. */

event init (i = 0)
{
  geometric_beta (1./3., true);
  
  foreach() {
    zb[] = -0.5;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
}

/**
We need to add the horizontal viscosity. This is important. To do it
correctly we copy from an [existing
example](/src/examples/breaking.c). */

event viscous_term (i++) {
  foreach()
    vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.);
  horizontal_diffusion ({u, w}, nu, dt);
}

/**
We plots some free surface profiles, as a visual check. */

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
We compute and log the kinetic, potential and total (normalised)
energies. */

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
  double E = ke/2. + g_*gpe + 0.125;
  static double E0 = 0.;
  if (i == 0)
    E0 = E;
  fprintf (stderr, "%g %g %g %g\n", t, ke/2., g_*gpe + 0.125, E/E0);
}
