/**
# Viscous dissipation of a linear wave

This is a work-in-progress version of [/src/test/viscous.c]() where we
are trying to address the [correct free-surface stress
condition](#free-surface-stress-free-condition).

As shown below we now agree much better with Lamb.

~~~gnuplot Evolution of the numerical and theoretical energies.
set xlabel 't/T_0'
set ylabel 'Normalised energy E/E_0'
k = 2.*pi
nu = 1./20000.
lambda = 1.
g = 1.
T0 = k*lambda/sqrt(g*k)
set xrange [0:10]
plot 'log' u ($1/T0):4 every 30 w p t 'Numerical', \
     exp(-2.*nu*k**2*x*T0) w l t 'Joseph and Wang', \
     exp(-4.*nu*k**2*x*T0) w l t 'Lamb'
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

## See also

* [Equivalent case with the Navier-Stokes solver](viscous-ns.c)
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
The wavenumber is $k=2\pi$, the amplitude/slope $ak$ is small enough
(decreasing it has little effect on the results). The viscosity is
$\nu=1./RE$. Decreasing the viscosity (or increasing it) does not
affect the results qualitatively.

Note that we set the fluid depth to h_ = 2 (i.e. twice the wavelength)
to be in the deep water limit. This is important to match the
theory. Further increasing the depth does not change the results. */

double k_ = 2.*pi, h_ = 2., g_ = 1., ak = 0.01;
double RE = 20000.;

/**
This is the wave period. */

#define T0  (k_*L0/sqrt(g_*k_))

int main()
{
  periodic (right);
  G = g_;
  nu = 1./RE;

  /**
  We adjust the numerical parameters for convergence. The vertical
  resolution matters, as do the convergence of the multigrid solver
  and timestep. The values below are converged i.e. increasing N, nl
  and NITERMIN or decreasing DT has little effect on the results. */
  
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
## Free-surface stress-free condition

The default (stress) condition on the free-surface is "no-slip"
i.e. $\partial_z u_x = 0$. This is not correct since the (linearised)
free-surface stress-free condition should be $\partial_z u_x +
\partial_x w = 0$. So here we redefine the stress-free condition as
$\partial_z u_x = - \partial_x w$ where we approximate the
free-surface vertical component with the vertical velocity in the last
layer (nl - 1).

At the free-surface, the boundary condition for the diffusion of the
vertical component of velocity must verify the incompressibility
condition, i.e.
$$
\partial_z w = - \partial_x u
$$
*/

vector dzu[];
scalar dzw[];

event viscous_term (i++) {
  dut = dzu;
  dwt = dzw;
  foreach() {
    dzu.x[] = - (w[1,0,nl-1] - w[-1,0,nl-1])/(2.*Delta);
    dzw[] = - (u.x[1,0,nl-1] - u.x[-1,0,nl-1])/(2.*Delta);
  }
  
  /**
  We need to add the horizontal viscosity. This is also important. */

  horizontal_diffusion ({u, w}, nu, dt);
}

/**
## Viscous pressure jump on the free-surface

To recover Lamb's dissipation, the viscous pressure jump must also be
taken into account i.e. the pressure at the (linearised) surface is
$$
p_s = g\eta + 2\nu\partial_z w
$$
to take this jump into account in the horizontal momentum, we need to add
$$
- \partial_x ( 2\nu\partial_z w )
$$
to the acceleration. */

event acceleration (i++)
{
  foreach_face() {
    double dp = - 2.*nu*(dzw[] - dzw[-1])/Delta;
    foreach_layer()
      ha.x[] += hf.x[]*dp;
  }
}

/**
We plot some free surface profiles, as a visual check. We only do 10
wave periods. This is largely sufficient to measure dissipation.

This confirms that the wave evolution matches the theoretical wave
period. The decrease in amplitude is also visible.

~~~gnuplot Evolution of the free surface profiles
set xlabel 'x'
set ylabel 'z'
set xrange [0:1]
plot 'out' w l t ''
~~~

To check whether $\partial_x w$ and $\partial_z u$ at the free-surface
agree with the theory we plot them below. While $\partial_x w$ agrees
well, $\partial_z u$ does not. We checked that when viscosity is zero,
both agree with the (inviscid) theory. So there is something to
understand here.

~~~gnuplot $\partial_x w$ and $\partial_z u$ at the free surface, and the theoretical value.
set ylabel ''
plot 'out' u 1:7 every 4 w p t 'u_z', '' u 1:8 every 4 w p t 'w_x', \
     '' u 1:9 w l t 'theory'
~~~
*/

event profiles (t = T0/4.; t += T0/4.; t <= 10.*T0) {
  double omega = sqrt(g_*k_*tanh(k_*h_));
  double a_ = ak/k_, U = a_*omega/sinh(k_*h_);
  foreach (serial) {
    double H = zb[];
    foreach_layer()
      H += h[];
    printf ("%g %g %g %g %g %g %g %g %g\n", x, H,
            u.x[0,0,nl-1], w[0,0,nl-1],
            U*cosh(k_*h_)*cos(k_*x - omega*t),
            U*sinh(k_*h_)*sin(k_*x - omega*t),            
            2.*(u.x[0,0,nl-1] - u.x[0,0,nl-2])/(h[0,0,nl-1] + h[0,0,nl-2]),
            (w[1,0,nl-1] - w[-1,0,nl-1])/(2.*Delta),
            U*k_*sinh(k_*h_)*cos(k_*x - omega*t));
  }
  printf ("\n\n");
}

/**
We compute and log the kinetic, potential and total (normalised)
energies. The initial condition seems to lead to an initial "jump"
(i.e. relaxation) in the numerical solution, so we skip the first
timestep. This may be improved. */

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
