/**
# Rotating convection in Europa's subsurface ocean

First-step reproduction of [Hay et al.,
2023](../REPORT.md) (JGR Planets, 128, e2022JE007648) with the
hydrostatic [multilayer solver](/src/layered/README) of Basilisk,
structured after [held-suarez.c](/src/examples/held-suarez.c).

The ocean is a 100 km deep spherical shell on a 1561 km radius body,
rotating with period 306 000 s (one Europan day). It is heated from
below and cooled at the top with a fixed temperature contrast $\Delta
T$, and experiences quadratic turbulent drag $\tau = \rho_0 c_D
|\mathbf{u}|\mathbf{u}$, $c_D = 0.002$, at *both* the seafloor and the
ice-ocean interface.

All parameters below are those of Table 1 of the paper, giving Ek =
3e-4, Pr = 1; the Rayleigh number is set through $\Delta T$ (Table 2
of the paper): $Ra = \alpha g \Delta T D^3/(\nu\kappa) = 6.85\times
10^7 \Delta T[K]$. */

#include "grid/multigrid.h"
#include "spherical.h"
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "bderembl/libs/netcdf_bas.h"
/**
## Coriolis acceleration and quadratic drag

The Coriolis parameter takes its standard (traditional approximation)
form. */

const double Omega = 2.*pi/306000.;
#define F0() (2.*Omega*sin(y*pi/180.))

/**
The quadratic drag (eq. 2 of the paper) is implemented through the
(linear) friction coefficient `K0()` of
[coriolis.h](/src/layered/coriolis.h): making it velocity-dependent,
$K_0 = c_D |\mathbf{u}|/h$, gives the quadratic stress $\tau/\rho_0 =
c_D|\mathbf{u}|\mathbf{u}$ as a deceleration of the boundary
layer. It is applied in the bottom ($l = 0$) and top ($l =$ nl - 1)
layers only, following the pattern of
[global-tides.c](/src/examples/global-tides.c) and
[ocean.h](/src/examples/ocean.h). */

const double Cd = 2e-3;
#define K0() ((point.l == 0 || point.l == nl - 1) ?		\
	      (h[] > dry ? Cd*norm(u)/h[] : HUGE) : 0.)
#include "layered/coriolis.h"

/**
## Temperature and buoyancy

The linear equation of state of the paper (density depends on
potential temperature only) is
$$
\Delta \rho = - \alpha (T - T_0)
$$
given to the [Boussinesq buoyancy module](/src/layered/dr.h). */

double T0 = 0.;
const double alphaT = 2e-4;
#define drho(T) (- alphaT*((T) - T0))
#include "layered/dr.h"

#include "layered/diffusion.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "profiling.h"

/**
## Physical and numerical constants

Values from Table 1 of the paper. */

const double DEPTH = 1.e5;            // ocean thickness (m)
const double rho0  = 1000.;           // reference density (kg/m^3)
const double kappa = 61.6;            // thermal diffusivity (m^2/s)
const double eday  = 306000.;         // rotation period (s)

/**
$\Delta T$ sets the Rayleigh number (default: Ra = 6.67e6, the middle
of the paper's range). `tend` is the run duration in Europan
days. Both can be changed from the command line (see `main()`). */

double DELTAT = 0.0973, tend = 100.;

int main (int argc, char * argv[])
{

  /**
  The global domain, in spherical coordinates with the radius of
  Europa setting the length unit (meters). */

  dimensions (2, 1);
  size (360.);
  origin (-180., -90.);
  periodic (right);
  Radius = 1561000. [1];
  G = 1.3 [1,-2];

  if (argc > 1) N = atoi (argv[1]);   else N = 256;
  if (argc > 2) DELTAT = atof (argv[2]);
  if (argc > 3) DT = atof (argv[3]);  else DT = 600. [0,1];
  if (argc > 4) tend = atof (argv[4]);

  DT = 50; //150.[0,1] * 256./N;

  nl = 10;
  nu = kappa; // vertical viscosity (Pr = 1)
  TOLERANCE = 0.1; // as in held-suarez.c

  double Ra = alphaT*G*DELTAT*cube(DEPTH)/sq(kappa);
  fprintf (stderr, "# Ra = %g, Ek = %g, N = %d, nl = %d, DT = %g, tend = %g edays\n",
	   Ra, kappa/(Omega*sq(DEPTH)), N, nl, DT, tend);

  run();
}

/**
## Initial conditions */

event init (i = 0)
{
  if (restore ("restart"))
    event ("metric");
  else {

    /**
    Flat seafloor, uniform sigma layers, linear unstable stratification
    (hot at the bottom, cold at the top) plus a random perturbation of
    order $10^{-3}\Delta T$ to seed convection, as in the paper. */

    foreach(cpu) {
      zb[] = - DEPTH;
      double z = zb[];
      foreach_layer() {
        h[] = DEPTH/nl;
        z += h[]/2.;
        double zz = (z - zb[])/DEPTH; // 0 at the bottom, 1 at the top
        T[] = T0 + DELTAT*(0.5 - zz) + 1e-3*DELTAT*noise();
        z += h[]/2.;
      }
    }
    create_nc({zb, h, eta, u, T}, "out.nc");
  }

  /**
  The vertical viscosity must be free-slip at the bottom so that the
  quadratic drag above is the *only* boundary stress (see section 2.4
  of the paper). The top is already free-slip by default in
  [diffusion.h](/src/layered/diffusion.h). A very large Navier slip
  length gives effective free-slip. */

  // foreach()
  //   foreach_dimension()
  //     lambda_b.x[] = 1.e12;
  lambda_b[] = {HUGE,HUGE,HUGE};
  /**
  The poles are closed with dry land on the non-periodic (latitude)
  boundaries, as in [ocean.h](/src/examples/ocean.h). */

  u.t[top] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  zb[top] = 100.;
  zb[bottom] = 100.;
  h[top] = 0;
  h[bottom] = 0;
}

/**
## Thermal forcing

Following the paper, the top and bottom boundary temperatures are
held fixed by relaxing the top and bottom layers towards $T_0 \mp
\Delta T/2$ on a timescale of twice the timestep ("in effect constant
temperatures"). */

event thermal_bc (i++)
{
  const double tau = 2.*DT;
  foreach() {
    foreach_layer()
      if (point.l == 0) {
        const double Tb = T0 + DELTAT/2.;
        T[] = (T[] + dt/tau*Tb)/(1. + dt/tau);
      }
      else if (point.l == nl - 1) {
        const double Tt = T0 - DELTAT/2.;
        T[] = (T[] + dt/tau*Tt)/(1. + dt/tau);
      }
  }
}

/**
## Diffusion and convective adjustment

Vertical thermal diffusion ($\kappa$ = 61.6 m^2^/s, implicit)
transports heat away from the pinned boundaries. Vertical momentum
viscosity ($\nu = \kappa$, Pr = 1) is applied by the
[diffusion.h](/src/layered/diffusion.h) module itself (with the
free-slip bottom set above). Horizontal diffusion of momentum and
temperature (also 61.6 m^2^/s) is explicit: the stability limit
$\Delta^2/4\kappa \sim 10^{10}$ s is never approached.

Since the hydrostatic layered model cannot overturn, statically
unstable columns are homogenised by thickness-weighted pairwise
mixing: this is the surrogate for the paper's resolved convection
plumes. */

event convection (i++, last)
{
  foreach()
    vertical_diffusion (point, h, T, dt, kappa, 0., 0., HUGE);

  horizontal_diffusion ({u.x, u.y, T}, kappa, dt);

	//  foreach() {
	//    double hv[nl], Tv[nl];
	//    int l = 0;
	//    foreach_layer() { hv[l] = h[], Tv[l] = T[]; l++; }
	//    for (int it = 0; it < nl; it++) {
	//      bool stable = true;
	//      for (l = 0; l < nl - 1; l++)
	// if (drho(Tv[l+1]) > drho(Tv[l])) { // dense over light
	//   double m = (hv[l]*Tv[l] + hv[l+1]*Tv[l+1])/(hv[l] + hv[l+1]);
	//   Tv[l] = Tv[l+1] = m;
	//   stable = false;
	// }
	//      if (stable) break;
	//    }
	//    l = 0;
	//    foreach_layer() { T[] = Tv[l]; l++; }
  //}
}

/**
## Diagnostics

Global kinetic energy monitors the spin-up towards a statistically
steady state (the convergence criterion used in the paper), and `Tq`
is the instantaneous axial ice-ocean torque from the quadratic stress
in the top layer,
$$
\mathcal{T}_z = \sum_\text{cells} \rho_0 c_D |\mathbf{u}_h| u_\phi
\cos\phi\, R\, dA
$$
i.e. eq. 5 of the paper without the time/zonal averaging. `HM` is the
total water mass (conservation check). */

event logfile (i++)
{
  double KE = 0., Tq = 0., HM = 0., umax = 0., Tmin = 1e30, Tmax = -1e30;
  foreach (reduction(+:KE) reduction(+:Tq) reduction(+:HM) \
	   reduction(max:umax) reduction(min:Tmin) reduction(max:Tmax)) {
    foreach_layer() {
      double u2 = sq(u.x[]) + sq(u.y[]);
      umax = max (umax, sqrt(u2));
      Tmin = min (Tmin, T[]), Tmax = max (Tmax, T[]);
      KE  += 0.5*u2*h[]*dv();
      HM  += h[]*dv();
      if (point.l == nl - 1)
	Tq += rho0*Cd*sqrt(u2)*u.x[]*(Radius*cos(y*pi/180.))*dv();
    }
  }
  fprintf (stderr, "%g %g %g %g %g %g %g %g\n",
	   t, dt, KE, Tq, HM, umax, Tmin, Tmax);
  write_nc();
}

/**
## Movies and snapshots

Maps of the top-layer zonal velocity and temperature. */

scalar ut[], Tt[];

event movies (t += eday/2.)
{
  foreach() {
    point.l = nl - 1;
    ut[] = u.x[];
    Tt[] = T[];
  }
  output_ppm (ut, file = "ut.mp4", n = 512, linear = true, map = cool_warm);
  output_ppm (Tt, file = "Tt.mp4", n = 512, linear = true, map = cool_warm);
}

event snapshot (t += 10.*eday) {
  dump();
}

event ending (t = tend*eday);

/**
## Main setup

Command-line parameters: `N` (default 256), `DELTAT` (K, default
0.0973 i.e. Ra = 6.67e6), `DT` (s, default 600), `tend` (Europan
days, default 100). */


