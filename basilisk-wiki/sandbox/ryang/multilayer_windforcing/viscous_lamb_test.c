/**
# Viscous dissipation of a linear wave (By ChatGPT and only for temporary use)

The original zero top gradients and zero surface gauge pressure are
replaced by:

* $u_z=-w_x$: zero tangential viscous stress.
* $w_z=-u_x$: incompressibility, used in vertical diffusion of w.
* $(p_s-p_{atm})/\rho=2\nu w_z$: normal stress balance. The factor two
  follows from $\sigma_{zz}=-p+2\rho\nu w_z$.

Together these complete the viscous surface boundary conditions;
no empirical damping is added. Lamb's $4\nu k^2$ energy decay is the
weak-viscosity, deep-water limit.

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
/**
## Local hook for vertical-velocity diffusion

The stock NL header passes zero as the top w-gradient to vertical_diffusion.
Only while including that header, redirect its call to the wrapper below.
The wrapper invokes the original diffusion routine with w_z=-u_x. This
avoids editing/copying nh.h or adding a second vertical-diffusion step.
Horizontal-velocity diffusion, already included through hydro.h, is unchanged.
The wrapper has the same argument list as the upstream function; dst is
replaced only in this NH call. */
scalar clean_w_gradient[];

void clean_vertical_diffusion (Point point, scalar h, scalar s,
                               double dt, double D, double dst,
                               double s_b, double lambda_b)
{
  vertical_diffusion (point, h, s, dt, D,
                      clean_w_gradient[], s_b, lambda_b);
}

#define vertical_diffusion clean_vertical_diffusion
#include "layered/nh.h"
#undef vertical_diffusion
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

double k_ = 2.*pi [-1], h_ = 2. [1], g_ = 1. [1,-2], ak = 0.01;
double RE = 20000. [-2,1]; // Inverse kinematic viscosity, not dimensionless Re.
vector dutf[];
scalar surface_u[], surface_w[], surface_pressure[];

/**
This is the wave period. */

#define T0  (k_*L0/sqrt(g_*k_))

int main()
{
  periodic (right);
  G = g_;
  nu = 1./RE;
  
  N = 256;
  nl = 60;
  NITERMIN = 3;
  DT = 1e-2 [0,1];
  
  dut = dutf;
  
  run();
}

/**
We optimise layers using a geometric distribution to better resolve
the surface and initialize the wave height and velocity
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

Zero u_z is a free-slip approximation, NOT no-slip. A clean interface
instead requires u_z+w_x=0, together with the normal stress balance.

The top velocities are reconstructed from the two top layer averages.
For thicknesses ht and hb, linear extrapolation to the surface gives
u_s=((2 ht+hb) u_top-ht u_below)/(ht+hb), and similarly for w_s.
All surface gradients use centred differences on this periodic grid. */

void update_surface_stress (void)
{
  foreach() {
    double ht = h[0,0,nl-1], hb = h[0,0,nl-2];
    surface_u[] = ((2.*ht+hb)*u.x[0,0,nl-1]
                   -ht*u.x[0,0,nl-2])/(ht+hb);
    surface_w[] = ((2.*ht+hb)*w[0,0,nl-1]
                   -ht*w[0,0,nl-2])/(ht+hb);
  }
  boundary ({surface_u, surface_w});
  foreach() {
    dutf.x[] = -(surface_w[1]-surface_w[-1])/(2.*Delta);
    clean_w_gradient[] = -(surface_u[1]-surface_u[-1])/(2.*Delta);
    surface_pressure[] = 2.*nu*clean_w_gradient[];
  }
  boundary ({dutf, clean_w_gradient, surface_pressure});
}

event viscous_term (i++)
{
  // Inherited events run this update before native vertical diffusion.
  update_surface_stress();
  // Keep both horizontal velocity components' Laplacians with the SAME nu.
  horizontal_diffusion ({u, w}, nu, dt);
}

event acceleration (i++)
{
  update_surface_stress();
  // Add the prescribed surface pressure gradient before the inherited
  // implicit/NH projection. Its interior pressure correction enforces
  // incompressibility. Do not overwrite phi at the top layer centre.
  foreach_face() {
    double ap = -(surface_pressure[]-surface_pressure[-1])/Delta;
    foreach_layer()
      ha.x[] += hf.x[]*ap;
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

*/

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
energies, using the original four-column output and normalisation at t=0.
The Stokes initial condition can undergo a short viscous adjustment. */

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

/** Stop explicitly after ten periods, as in the original profile schedule. */
event stop (t = 10.*T0) { return 1; }
