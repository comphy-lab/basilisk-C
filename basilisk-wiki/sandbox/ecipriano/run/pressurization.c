/**
# Pressurization of a Closed Tank

When the domain is closed, the volume of the system cannot change, and the
expansion of the fluid must be balanced by an increase of the thermodynamic
pressure $P_0$, which is uniform in space for a low Mach number system. These
test cases verify the pressurization model implemented in
[low-mach.h](../src/navier-stokes/low-mach.h) in configurations with an
analytic solution.

A closed tank is filled with a liquid layer at the bottom and with an ideal gas
ullage at the top. The phase change is switched off, the liquid is
incompressible, and the volume $V_g$ of the ullage does not change. The energy
supplied to the ullage gives:
$$
  \dfrac{dP_0}{dt} = \left(\gamma - 1\right)\dfrac{Q_g}{V_g}
$$
where $Q_g$ is the heat supplied to the gas phase per unit time. The heat is
supplied in two ways:

* **Volumetric source** (default): a uniform heat source $\dot{q}$ is applied
  to the ullage, and the thermal conductivity is null in both phases, therefore
  the ullage remains uniform, $Q_g/V_g = \dot{q}$ and
  $T_g(t) = T_g^0 + \gamma\dot{q}t/(\rho_g^0 c_{p,g})$. The velocity must
  remain zero everywhere.
* **Wall heat flux** (`-DHEATFLUX=1`): a uniform heat flux $q_w$ enters the
  tank through the walls, and it is transported by conduction. This
  configuration reproduces the passive pressurization of cryogenic tanks.

Each configuration is computed in two tanks: a square box of side $L_0$, and a
circle of radius $R$ inscribed in the box, obtained with an embedded boundary
(`-DSOLID=1`). With the interface at a height $H_l$ from the bottom of the box
($H_g = L_0 - H_l$), or at a distance $h$ above the center of the circle
($\theta = \arccos (h/R)$), the wall heat flux gives:
$$
  \dfrac{Q_g}{V_g} = \dfrac{q_w\left(L_0 + 2H_g\right)}{L_0 H_g}
  \qquad
  \dfrac{Q_g}{V_g} = \dfrac{2 q_w R\theta}{R^2\theta - h\sqrt{R^2 - h^2}}
$$

![Volumetric source (square box)](pressurization/movie.mp4)

![Volumetric source (embedded boundary)](pressurization-embed/movie.mp4)

![Wall heat flux (square box)](pressurization-heatflux/movie.mp4)

![Wall heat flux (embedded boundary)](pressurization-heatflux-embed/movie.mp4)
*/

/**
## Simulation Setup

We use the low Mach number Navier--Stokes solver, together with the variable
properties module. The phase change is switched off using the fixed flux model
with a null vaporization rate.
*/

#include "grid/multigrid.h"
#if SOLID
# include "embed.h"
#endif
#include "navier-stokes/low-mach.h"
#define P_ERR 1.e-6
#include "two-phase-varprop.h"
#include "basilisk-properties.h"
#include "two-phase.h"
#include "phasechange.h"
#include "fixedflux.h"
#include "view.h"

/**
### Case Parameters

The heating mechanism sets the duration of the simulation and the output
intervals. The finest grid can be increased (e.g. `-DMAXLEVEL=8`) for a
convergence study. */

#ifndef MAXLEVEL
# define MAXLEVEL 6
#endif

#if HEATFLUX
# define TEND 1.
# define LOG_EVERY 0.1
# define MOVIE_EVERY 0.01
# define DTMOVIE 100.
#else
# define TEND 0.1
# define LOG_EVERY 0.02
# define MOVIE_EVERY 0.002
# define DTMOVIE 15.
#endif

/**
### Boundary Conditions

Without heat flux, the default symmetry conditions make the box closed and
adiabatic. With the wall heat flux, the walls are no-slip and impermeable:
without the condition on the face velocity `uf`, the Poisson equation of `pf`
does not converge during the first time steps. The heat flux is imposed in the
`init` event, where the temperature fields of the phases are available. */

#if HEATFLUX
# if EMBED
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
# else
u.t[right] = dirichlet (0.);
u.t[left] = dirichlet (0.);
u.t[top] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);

uf.n[right] = dirichlet (0.);
uf.n[left] = dirichlet (0.);
uf.n[top] = dirichlet (0.);
uf.n[bottom] = dirichlet (0.);
# endif
#endif

/**
### Model Data

The volumetric heat source, the wall heat flux, the radius of the circular
tank, and the initial liquid level. */

int maxlevel;
double q0 = 1.e5, qw = 100., R0, level0 = 0.5;
double gamma0, mg0;

#define MWG 29.

/**
The density of the gas phase follows the ideal gas law. */

double gasprop_density_idealgas (void * p) {
  ThermoState * ts = p;
  return ts->P*MWG/(R_GAS*1.e3*ts->T);
}

int main (void) {
  NGS = 1, NLS = 1;

  /**
  We set the material properties of water and air. The thermal conductivity is
  null with the volumetric source, to keep the ullage uniform. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.8e-5;
  Dmix1 = 0., Dmix2 = 0.;
#if HEATFLUX
  lambda1 = 0.6, lambda2 = 0.026;
#else
  lambda1 = 0., lambda2 = 0.;
#endif
  cp1 = 4184., cp2 = 1004.5;
  dhev = 0.;

  P0 = 101325.;
  TG0 = 300., TL0 = 300.;

  closed = true;

  pcm.isomassfrac = true;
  pcm.divergence = true;
  mEvapVal = 0.;

  L0 = 0.01;
  R0 = 0.45*L0;
  DT = 1.e-3;

  gamma0 = cp2/(cp2 - R_GAS*1.e3/MWG);

  /**
  With the wall heat flux, the default tolerance of the Poisson solver gives
  velocity oscillations during the transient. The tolerance of the temperature
  equation is a residual, with the units of $\rho c_p T/\Delta t$. */

#if HEATFLUX
  TOLERANCE = 1e-6 [*];
  TOLERANCE_TEMPERATURE = 1e-2;
#endif

  for (maxlevel = 4; maxlevel <= MAXLEVEL; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
We initialize a flat interface which lies on a grid line on all the grids. */

#define circle(x,y,R) (sq(R) - sq(x - 0.5*L0) - sq(y - 0.5*L0))

event init (i = 0) {
#if EMBED
  solid (cs, fs, circle (x, y, R0));
  fractions_cleanup (cs, fs);
#endif
  fraction (f, level0*L0 - y);

  ThermoState tsl, tsg;
  tsl.T = TL0, tsl.P = P0, tsl.x = (double[]){1.};
  tsg.T = TG0, tsg.P = P0, tsg.x = (double[]){1.};

  phase_set_thermo_state (liq, &tsl);
  phase_set_thermo_state (gas, &tsg);

  phase_set_properties (liq, MWs = (double[]){18.});
  phase_set_properties (gas, MWs = (double[]){MWG});

  tp2.rhov = gasprop_density_idealgas;
  tp2.betaT = gasprop_thermal_expansion;
  tp2.chiT = gasprop_isothermal_compressibility;

  /**
  The heat flux is a Neumann condition on the temperature of both phases, on
  the boundaries of the domain or on the embedded boundary. */

#if HEATFLUX
  scalar TL = liq->T, TG = gas->T;
# if EMBED
  TL[embed] = neumann (qw/lambda1);
  TG[embed] = neumann (qw/lambda2);
# else
  for (int bid = 0; bid < nboundary; bid++) {
    TL[bid] = neumann (qw/lambda1);
    TG[bid] = neumann (qw/lambda2);
  }
# endif
#endif

  /**
  We store the initial mass of the ullage, computed with the virtual volumes
  of the cells (see [low-mach.h](../src/navier-stokes/low-mach.h)), which are
  consistent with the conservation properties of the solver. */

  mg0 = 0.;
  foreach (reduction(+:mg0))
    mg0 += tp2.rhov (&tsg)*(1. - f[])*dvv();
}

/**
### Volumetric Heat Source

The source is added in the `chemistry` event, after the source terms are reset
and before the `divergence` event, so that it enters both the velocity
divergence and the temperature equation. It is weighted on the volume fraction
of the gas phase. */

#if !HEATFLUX
event chemistry (i++) {
  scalar STexpG = gas->STexp;
  foreach()
    STexpG[] += q0*cmv()*(1. - f[]);
}
#endif

/**
### Pressurization

The thermodynamic pressure is integrated in time using the pressurization rate
computed by the projection step. */

event pressurization (i++) {
  P0 += dt*dP0dt;
}

/**
## Post-Processing

### Output Files

We write the time, the numerical and analytic thermodynamic pressure and
pressurization rate, the relative variation of the mass of the ullage, the
maximum velocity, and the volume-averaged temperature of the ullage. The
temperature of the gas phase is stored as a tracer. The mass and the volume of
the ullage are computed with the virtual metrics, consistently with the vof
transport in interfacial cells (ignoring the solid). */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

#if HEATFLUX
# if EMBED
  double h = (level0 - 0.5)*L0, theta = acos (h/R0);
  double dP0dtexact = (gamma0 - 1.)*qw*2.*R0*theta /
    (sq(R0)*theta - h*sqrt (sq(R0) - sq(h)));
# else
  double Hg = (1. - level0)*L0;
  double dP0dtexact = (gamma0 - 1.)*qw*(L0 + 2.*Hg)/(L0*Hg);
# endif
#else
  double dP0dtexact = (gamma0 - 1.)*q0;
#endif

  scalar TG = gas->T;
  double mg = 0., umax = 0., vg = 0., tg = 0.;
  foreach (reduction(+:mg) reduction(max:umax) reduction(+:vg)
      reduction(+:tg)) {
    double fg = 1. - f[];
    if (fg > F_ERR) {
      double rhog = P0*MWG/(R_GAS*1.e3*TG[]/fg);
      mg += rhog*fg*dvv();
      vg += fg*dvv();
      tg += TG[]*dvv();
    }
    umax = max (umax, norm (u));
  }

  fprintf (fp, "%g %g %g %g %g %g %g %g\n", t, P0, 101325. + dP0dtexact*t,
      dP0dt, dP0dtexact, (mg - mg0)/mg0, umax, tg/vg),
    fflush (fp);
}

/**
### Logger

We output the thermodynamic pressure of the tank in time (for testing). */

event logger (t += LOG_EVERY) {
  fprintf (stderr, "%d %.2f %.6g\n", maxlevel, t, P0);
}

event stop (t = TEND);

/**
### Movie

We write the animation with the evolution of the temperature field. */

event movie (t += MOVIE_EVERY) {
  if (maxlevel == 6) {
    clear();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f", lw = 2.);
#if EMBED
    draw_vof ("cs", "fs", filled = -1, fc = {1.,1.,1.});
    draw_vof ("cs", "fs", lw = 4.);
#else
    box (notics = true, lw = 4.);
#endif
    squares ("T", min = TG0, max = TG0 + DTMOVIE);
    save ("movie.mp4");
  }
}

/**
## Results

In all the figures, the left panel refers to the square box and the right panel
to the circular tank.

### Volumetric Heat Source

~~~gnuplot Evolution of the thermodynamic pressure
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key top left
set xlabel "t [s]"
set ylabel "P_0 [Pa]"
set size square

set title "Square box"
plot "../pressurization/OutputData-5" u 1:3 every 5 w p pt 6 lc rgb "black" \
       t "Analytic", \
     for [l=4:6] "../pressurization/OutputData-".l u 1:2 w l lw 2 lc l-3 \
       t "LEVEL ".l

set title "Embedded boundary"
plot "../pressurization-embed/OutputData-5" u 1:3 every 5 w p pt 6 \
       lc rgb "black" t "Analytic", \
     for [l=4:6] "../pressurization-embed/OutputData-".l u 1:2 w l lw 2 \
       lc l-3 t "LEVEL ".l
unset multiplot
~~~

Since the mass and the volume of the ullage do not change, the gas is heated
along the isochore $P_0 = \rho_g^0 R T_g$ (with $R$ the specific gas
constant).

~~~gnuplot P-T diagram of the ullage (volumetric source)
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "T_g [K]"
set ylabel "P_0 [Pa]"
set size square
isochore(T) = 101325.*T/300.

set title "Square box"
plot "../pressurization/OutputData-6" u 8:(isochore($8)) \
       w l lw 2 lc rgb "black" dt 2 \
       t "Ideal gas {/Symbol r} = {/Symbol r}_g^0", \
     for [l=4:6] "../pressurization/OutputData-".l u 8:2 \
       every 8 w p pt 2*l-6 lc l-3 t "LEVEL ".l

set title "Embedded boundary"
plot "../pressurization-embed/OutputData-6" u 8:(isochore($8)) \
       w l lw 2 lc rgb "black" dt 2 \
       t "Ideal gas {/Symbol r} = {/Symbol r}_g^0", \
     for [l=4:6] "../pressurization-embed/OutputData-".l u 8:2 \
       every 8 w p pt 2*l-6 lc l-3 t "LEVEL ".l
unset multiplot
~~~

The spurious velocity grows from round-off and saturates without affecting the
pressure of the ullage.

~~~gnuplot Spurious velocity
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "t [s]"
set ylabel "max |u| [m/s]"
set logscale y
set format y "10^{%T}"
set yrange [1e-20:1e-1]
set size square

set title "Square box"
plot for [l=4:6] "../pressurization/OutputData-".l u 1:7 w l lw 2 lc l-3 \
       t "LEVEL ".l

set title "Embedded boundary"
plot for [l=4:6] "../pressurization-embed/OutputData-".l u 1:7 w l lw 2 \
       lc l-3 t "LEVEL ".l
unset multiplot
~~~

### Wall Heat Flux

~~~gnuplot Pressurization rate
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "t [s]"
set ylabel "dP_0/dt [Pa/s]"
set size square

set title "Square box"
set yrange [15500:16500]
plot "../pressurization-heatflux/OutputData-5" u 1:5 w l lw 2 \
       lc rgb "black" t "Analytic", \
     for [l=4:6] "../pressurization-heatflux/OutputData-".l u 1:4 \
       w l lw 2 lc l-3 t "LEVEL ".l

set title "Embedded boundary"
set yrange [13000:18500]
plot "../pressurization-heatflux-embed/OutputData-5" u 1:5 w l lw 2 \
       lc rgb "black" t "Analytic", \
     for [l=4:6] "../pressurization-heatflux-embed/OutputData-".l u 1:4 \
       w l lw 2 lc l-3 t "LEVEL ".l
unset multiplot
~~~

The error of the circular tank is due to the polygonal approximation of the
tank given by the embedded boundary.

With the wall heat flux, the temperature of the ullage is not uniform and the
volume-averaged temperature lies slightly to the right of the isochore (the
pressure is about $10^{-4}$ below it at the end of the simulation, on all the
grids). The mass-averaged temperature, which for an ideal gas at uniform
pressure is the harmonic mean over the volume, lies exactly on the isochore as
long as the mass of the ullage is conserved.

~~~gnuplot P-T diagram of the ullage (wall heat flux)
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key bottom right
set xlabel "T_g [K]"
set ylabel "P_0 [Pa]"
set size square
isochore(T) = 101325.*T/300.

set title "Square box"
plot "../pressurization-heatflux/OutputData-6" u 8:(isochore($8)) \
       w l lw 2 lc rgb "black" dt 2 \
       t "Ideal gas {/Symbol r} = {/Symbol r}_g^0", \
     for [l=4:6] "../pressurization-heatflux/OutputData-".l u 8:2 \
       every 50 w p pt 2*l-6 lc l-3 t "LEVEL ".l

set title "Embedded boundary"
plot "../pressurization-heatflux-embed/OutputData-6" u 8:(isochore($8)) \
       w l lw 2 lc rgb "black" dt 2 \
       t "Ideal gas {/Symbol r} = {/Symbol r}_g^0", \
     for [l=4:6] "../pressurization-heatflux-embed/OutputData-".l u 8:2 \
       every 50 w p pt 2*l-6 lc l-3 t "LEVEL ".l
unset multiplot
~~~

~~~gnuplot Convergence of the pressurization rate (embedded boundary)
reset
set print "errors"
do for [l=4:6] {
  stats "<tail -n 1 ../pressurization-heatflux-embed/OutputData-".l \
    u (abs ($4 - $5)/$5) nooutput name "EMBED"
  print sprintf ("%d %.12f", 2**l, EMBED_mean)
}
unset print

set xlabel "Resolution"
set ylabel "Relative Error"
set xr[2**3:2**7]
set yr[1e-2:1]
set key bottom left
set size square
set grid
set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "errors" u 1:2 via a,b
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", a, b)

plot "errors" u 1:2 pt 8 lc 2 t "Embedded boundary", \
     f(x) w l lc 2 t ftitle(a, b)
~~~

~~~gnuplot Relative variation of the mass of the ullage
reset
set term @SVG size 900,450
set multiplot layout 1,2
set grid
set key top right
set xlabel "t [s]"
set ylabel "(m_g - m_g^0)/m_g^0"
set size square

set title "Square box"
plot for [l=4:6] "../pressurization-heatflux/OutputData-".l u 1:6 \
       w l lw 2 lc l-3 t "LEVEL ".l

set title "Embedded boundary"
plot for [l=4:6] "../pressurization-heatflux-embed/OutputData-".l u 1:6 \
       w l lw 2 lc l-3 t "LEVEL ".l
unset multiplot
~~~
*/
