/**
# Growth of an Electrolytic Bubble

A hydrogen bubble grows on a planar electrode, fed by the dissolved hydrogen
produced at the electrode surface. The setup reproduces the single-bubble
simulations of [Gennari et al., 2022](#gennari2022), Sec. 5, based on the
experiments of [Glas and Westwater, 1964](#glas1964). The bubble is pure
hydrogen and it is the *dispersed* phase, while the liquid is the *continuous*
phase, made of dissolved hydrogen and water, the inert species (see
[evaporation.h](../src/evaporation.h)). Hydrogen moves from the continuous to
the dispersed phase.

The electrode is the flat end of a wire of diameter $D_e$, on the left
boundary. The molar flux of hydrogen produced at the electrode follows from
Faraday's law, $J = I/(2F)$, with $I$ the current density. The liquid is
initially saturated, and the interface mass fraction follows from Henry's law
$\hat{c}_l = \hat{c}_g/He$, i.e. $\hat\omega = \rho_g/(He\rho_l)$. The bubble
nucleates with diameter $D_e/10$ and contact angle $\theta = 90^\circ$, so it
is a half sphere sitting on the electrode, and it is kept frozen for the
nucleation time $t_n$, while the concentration boundary layer develops. Then
the bubble grows following Scriven's law $R(t) = 2\beta\sqrt{\mathcal{D}t}$,
where the growth constant $\beta$ increases with the current density.

![Dissolved hydrogen around the growing bubble](electrolysis/movie.mp4)

To keep this test short, the domain is half the one of the paper ($12.5
D_b$ instead of $25 D_b$), the grid is coarser (level 6, $\Delta \approx
2.5~\mu$m, vs level 10), and the simulation stops at $t = 0.2$ s. The growth
constant changes by less than 1% compared to a level 8 simulation on the
full domain up to $t = 0.3$ s ($\beta = 0.327$). */

#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "view.h"

/**
### Problem data

The current density [A/m²] can be changed at compilation time. */

#ifndef CURRENT
# define CURRENT 300.
#endif

#ifndef MAXLEVEL
# define MAXLEVEL 6
#endif

#ifndef TMAX
# define TMAX 0.2
#endif

#define TNUC 0.02

int maxlevel = MAXLEVEL, minlevel = 4;
double De = 0.127e-3, R0, He = 53.3, Dl = 7.38e-9, MWH2;
double wint, gradw;

/**
### Boundary conditions

The axis of symmetry is the bottom boundary, the electrode is on the left
boundary, and the right boundary is an outflow. The left and top walls are
slip walls (the default symmetry conditions), which is consistent with the
$90^\circ$ contact angle. The species boundary conditions are set in *init*,
since the species fields are created in *defaults*. */

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

int main (void) {

  /**
  The gas phase contains only hydrogen, while the liquid contains dissolved
  hydrogen and water. */

  NGS = 1, NLS = 2;
  gas_species = (char *[]){"H2"};
  liq_species = (char *[]){"H2", "H2O"};

  /**
  Material properties from Table 8 of [Gennari et al., 2022](#gennari2022). As
  in the paper, the gas density is increased by a factor of 10 to speed up the
  pressure solver, and the molar mass is increased by the same factor, so that
  the volume change $\propto M/\rho_g$ is not altered. The surface tension is
  reduced by a factor $10^4$, which is enough to keep the bubble spherical
  without gravity. */

  rho1 = 996., rho2 = 0.8;
  mu1 = 8.32e-4, mu2 = 8.96e-6;
  Dmix1 = Dl, Dmix2 = Dl;
  MWH2 = 2.016e-2;

  f.sigma = 0.075e-4;

  /**
  The molecular weights are equal, so that mole and mass fractions coincide,
  and the Henry constant becomes the interface mass fraction *YIntVal*. */

  MW1 = 1., MW2 = 1.;
  wint = rho2/(He*rho1);
  YIntVal = wint;

  /**
  The mass flux of hydrogen produced at the electrode $J M$ gives the normal
  gradient of the dissolved hydrogen mass fraction. */

  double Faraday = 96485.3;
  gradw = CURRENT/(2.*Faraday)*MWH2/(rho1*Dl);

  /**
  Single velocity field, with the phase tracers advected consistently with
  the volume fraction. The system is isothermal. */

  nv = 1;
  pcm.consistent = true;
  pcm.isothermal = true;

  R0 = 0.5*De/10.;
  L0 = 12.5*2.*R0;
  TOLERANCE = 1e-6;

  init_grid (1 << maxlevel);
  run();
}

/**
The default settings of [evaporation.h](../src/evaporation.h) are for
droplets. For a bubble, the expansion term is shifted to the gas, so that the
liquid carries the Stefan flow, and the interface regression uses the liquid
density. This event runs after all the *init* events. During the nucleation
stage the bubble volume is frozen: the interface is not moved by the phase
change and the expansion term is not added, but the bubble still consumes the
dissolved gas through the species source terms. */

event end_init (i = 0) {
  pcm.shifting = SHIFT_TO_GAS;
  pcm.byrhogas = false;
  pcm.advection = NO_ADVECTION;
  pcm.expansion = false;
}

event nucleation (t = TNUC) {
  pcm.advection = ADVECTION_VELOCITY;
  pcm.expansion = true;
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, -circle (x, y, R0));

  scalar YLH2 = liq->YList[0], YLH2O = liq->YList[1];
  scalar YGH2 = gas->YList[0];
  foreach() {
    YLH2[] = f[]*wint;
    YLH2O[] = f[]*(1. - wint);
    YGH2[] = 1. - f[];
  }

  /**
  On the electrode, the outward normal gradient is set where the wall is wet:
  the hydrogen enters the liquid, and the water mass fraction compensates. */

  YLH2[left] = neumann (y < 0.5*De ? gradw*f[] : 0.);
  YLH2O[left] = neumann (y < 0.5*De ? -gradw*f[] : 0.);
}

#if TREE
event adapt (i++) {
  scalar YLH2 = liq->YList[0];
  adapt_wavelet_leave_interface ({YLH2, u.x, u.y}, {f},
      (double[]){5e-2*wint, 1e-6, 1e-6}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-processing

The bubble is a half sphere on the electrode, so its radius follows from the
volume in the axisymmetric domain. */

double effective_radius() {
  double vol = 0.;
  foreach (reduction(+:vol))
    vol += (1. - f[])*dv();
  return cbrt (3.*vol);
}

event output_data (t += 0.002) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  /**
  We also write the supersaturation on the electrode, outside the bubble. */

  scalar YLH2 = liq->YList[0];
  double xp = 0.5*L0/(1 << maxlevel), yp = 0.375*De;
  double ff = interpolate (f, xp, yp);
  fprintf (fp, "%g %g %g %g\n", t, sqrt (t), effective_radius(),
      (ff > 1e-6) ? interpolate (YLH2, xp, yp)/ff/wint : 0.);
  fflush (fp);
}

event logger (t += 0.02) {
  fprintf (stderr, "%d %.2f %.3g\n", i, t, effective_radius());
}

#if MOVIE
event movie (t += 0.004; t <= TMAX) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("YLH2", min = 0., max = 3.*wint, linear = true);
  save ("movie.mp4");
}
#endif

event stop (t = TMAX);

/**
## Results

The growth constant is the slope of $R$ vs $\sqrt{t}$ divided by
$2\sqrt{\mathcal{D}}$, fitted once the bubble grows steadily.

~~~gnuplot Bubble radius
reset
set xlabel "t^{1/2} [s^{1/2}]"
set ylabel "R [μm]"
set yrange [0:*]
set key top left
set size square
set grid
D = 7.38e-9
f(x) = a*x + b
fit [0.3:] f(x) "OutputData-6" u 2:3 via a, b
x1 = 0.28; x2 = 0.37; dy = 2.
y1 = f(x1)*1e6 - dy; y2 = f(x2)*1e6 - dy
set arrow from x1,y1 to x2,y1 nohead lw 1.5 lc -1
set arrow from x2,y1 to x2,y2 nohead lw 1.5 lc -1
set arrow from x1,y1 to x2,y2 nohead lw 1.5 lc -1
set label sprintf ("β = %.3f", a/(2.*sqrt(D))) \
  at 0.5*(x1 + x2), y1 - 1. center
plot "OutputData-6" every 2 u 2:($3*1e6) w p pt 6 t "LEVEL 6", \
     f(x)*1e6 w l lw 2 lc -1 t "linear fit"
~~~

## References

~~~bib
@article{gennari2022,
  title={A phase-change model for diffusion-driven mass transfer problems in
    incompressible two-phase flows},
  author={Gennari, Gabriele and Jefferson-Loveday, Richard and Pickering,
    Stephen J},
  journal={Chemical Engineering Science},
  volume={259},
  pages={117791},
  year={2022},
  publisher={Elsevier}
}

@article{glas1964,
  title={Measurements of the growth of electrolytic bubbles},
  author={Glas, JP and Westwater, JW},
  journal={International Journal of Heat and Mass Transfer},
  volume={7},
  number={12},
  pages={1427--1443},
  year={1964},
  publisher={Elsevier}
}
~~~
*/
