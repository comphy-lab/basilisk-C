/**
# Bagnold Piston

Active pressurization of a channel of length $l$, filled with liquid on the left
and with an ideal gas ullage on the right. Liquid is injected from the left
boundary with mass flux $\dot{m}_{in}$ and pushes the interface like a piston,
with velocity $u_\Gamma = \dot{m}_{in}/\rho_l$. Without heat conduction, the
ullage is compressed isentropically:
$$
  P_0(t) = \left(\dfrac{l - x_{\Gamma,0}}{l - x_{\Gamma,0} - u_\Gamma t}
    \right)^\gamma P_0(0)
$$
$$
  T_g(t) = \left(\dfrac{P_0(0)}{P_0(t)}\right)^{\frac{1-\gamma}{\gamma}} T_0
$$
Compared to [pressurization.c](pressurization.c), the pressurization model of
[low-mach.h](../src/navier-stokes/low-mach.h) is driven by an inflow boundary
instead of a heat flux.

![Evolution of the interface and the gas phase
temperature](piston/movie.mp4)(width="100%")
*/

/**
## Simulation Setup

The channel is a multigrid with aspect ratio $l/h = 10$. The phase change is
switched off using the fixed flux model with a null vaporization rate. */

#include "grid/multigrid.h"
#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Model Data

A hydrogen-like ideal gas. The initial interface position $x_{\Gamma,0} = 0.3\,l$
reproduces the reference temperature rise (20.4 K to 24.3 K in 100 s). */

#define RGAS    5049.     // Specific gas constant [J/kg/K]
#define GAMMA   1.78      // Heat capacity ratio of the gas [-]
#define PIN     103.e3    // Initial thermodynamic pressure [Pa]
#define TIN     20.4      // Initial and inflow temperature [K]
#define MDOTIN  0.1       // Inflow mass flux [kg/m2/s]
#define RHORAT  70.       // Liquid to gas density ratio [-]
#define LENGTH  1.        // Length of the channel [m]
#define HEIGHT  0.1       // Height of the channel [m]
#define XGAMMA0 0.3       // Initial position of the interface [m]

int maxlevel;
double rhog0, uin;

/**
The gas density follows the ideal gas law. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = p;
  return ts->P/(RGAS*ts->T);
}

/**
### Boundary Conditions

Liquid enters from the left; the other boundaries are adiabatic slip walls
(default). */

u.n[left] = dirichlet (uin);
u.t[left] = dirichlet (0.);
f[left] = dirichlet (1.);
T[left] = dirichlet (TIN);

int main (void) {

  NGS = 1, NLS = 1;

  /**
  Null conductivity and diffusivity give an adiabatic compression. Viscosities
  and liquid heat capacity are those of hydrogen at 20 K and do not affect the
  plug flow. */

  P0 = PIN;
  TG0 = TIN, TL0 = TIN;

  rhog0 = P0/(RGAS*TIN);
  rho1 = RHORAT*rhog0, rho2 = rhog0;
  mu1 = 1.3e-5, mu2 = 1.1e-6;
  Dmix1 = 0., Dmix2 = 0.;
  lambda1 = 0., lambda2 = 0.;
  cp1 = 9.7e3, cp2 = GAMMA*RGAS/(GAMMA - 1.);
  dhev = 0.;

  uin = MDOTIN/rho1;

  /**
  The domain volume is fixed: the inflow rate is passed to the closed-system
  projection as a negative net outflow, which compresses the ullage. */

  closed = true;
  closed_flux = -uin*HEIGHT;

  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  /**
  The divergence in the ullage is $\sim u_\Gamma/l \approx 10^{-3}$ s$^{-1}$,
  hence the reduced tolerance. */

  TOLERANCE = 1.e-6;
  DT = 0.1;

  L0 = LENGTH;
  int nx = round (LENGTH/HEIGHT);
  dimensions (nx = nx, ny = 1);

  for (maxlevel = 3; maxlevel <= 5; maxlevel++) {
    init_grid (nx*(1 << maxlevel));
    run();
  }
}

/**
### Initial Conditions

The interface is flat at $x_{\Gamma,0}$. The velocity is initialized with the
exact solution: starting from rest, the Poisson problem is not compatible with
`closed_flux` during the first time steps. */

event init (i = 0) {
  fraction (f, XGAMMA0 - x);

  foreach()
    u.x[] = (x < XGAMMA0) ? uin : uin*(LENGTH - x)/(LENGTH - XGAMMA0);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = P0, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = P0, tsg.x = (double[]){1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){R_GAS*1.e3/RGAS});
  phase_set_properties (gas, MWs = (double[]){R_GAS*1.e3/RGAS});

  tp2.rhov = gasprop_density_idealgas;
  tp2.betaT = gasprop_thermal_expansion;
  tp2.chiT = gasprop_isothermal_compressibility;
}

/**
### Thermodynamic Pressure

$P_0$ is integrated in time with the rate computed by the projection. */

event pressurization (i++) {
  P0 += dt*dP0dt;
}

/**
## Post-Processing

We write $P_0$, the average ullage temperature, the interface position and the
mass of gas, together with the analytic values. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  /**
  `TG` is a tracer (multiplied by $1 - f$). The gas mass uses
  $\rho_g = P_0/(R_g T_g)$ rather than the density field, which lags one step
  behind the advected volume fraction. */

  scalar TG = gas->T;
  double vliq = 0., vgas = 0., tgas = 0., mgas = 0.;
  foreach (reduction(+:vliq) reduction(+:vgas) reduction(+:tgas)
      reduction(+:mgas)) {
    vliq += f[]*dv();
    vgas += (1. - f[])*dv();
    tgas += TG[]*dv();
    if (1. - f[] > F_ERR && TG[] > 0.)
      mgas += P0/RGAS*sq(1. - f[])/TG[]*dv();
  }
  double Tgas = tgas/vgas;
  double xgamma = vliq/HEIGHT;

  static double mgas0 = 0.;
  if (i == 0)
    mgas0 = mgas;

  double P0exact = pow ((LENGTH - XGAMMA0)/(LENGTH - XGAMMA0 - uin*t),
      GAMMA)*PIN;
  double Tgexact = pow (PIN/P0exact, (1. - GAMMA)/GAMMA)*TIN;

  fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, P0, P0exact, Tgas, Tgexact,
      xgamma, XGAMMA0 + uin*t, (mgas - mgas0)/mgas0),
    fflush (fp);
}

/**
We log $P_0$ and the ullage temperature (for testing), and we write the movie. */

event logger (t += 10.) {
  double Tgas = interpolate (gas->T, 0.9*LENGTH, 0.5*HEIGHT);
  fprintf (stderr, "%d %g %.6g %.6g\n", maxlevel, t, P0, Tgas);
}

event movie (t += 1; t <= 100) {
  if (maxlevel == 5) {
    clear();
    view (fov = 4.2, tx = -0.5, ty = -0.5*HEIGHT, width = 4000, height = 800);
    box (notics = true);
    draw_vof ("f", lw = 4.);
    squares ("T", min = TIN, max = 24.5, linear = true);
    save ("movie.mp4");
  }
}

/**
## Results

The pressure, the ullage temperature and the interface position follow the
analytic solution. Mass conservation improves with refinement.

~~~gnuplot Evolution of the thermodynamic pressure
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "P_0 [kPa]"
set size square

plot "OutputData-4" u 1:($3/1e3) every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:($2/1e3) w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:($2/1e3) w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:($2/1e3) w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

~~~gnuplot Evolution of the ullage temperature
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "T_g [K]"
set size square

plot "OutputData-4" u 1:5 every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:4 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:4 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:4 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

~~~gnuplot Position of the interface
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "x_Γ [m]"
set size square

plot "OutputData-4" u 1:7 every 50 w p ps 0.9 pt 6 lc rgb "black" \
       t "Analytic", \
     "OutputData-3" u 1:6 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:6 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:6 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~

~~~gnuplot Relative variation of the gas mass
reset
set grid
set key top left
set xlabel "t [s]"
set ylabel "(m_g - m_g^0)/m_g^0 [-]"
set size square

plot "OutputData-3" u 1:8 w l lw 2 lc 1 t "LEVEL 3", \
     "OutputData-4" u 1:8 w l lw 2 lc 2 t "LEVEL 4", \
     "OutputData-5" u 1:8 w l lw 2 dt 2 lc 3 t "LEVEL 5"
~~~
*/
