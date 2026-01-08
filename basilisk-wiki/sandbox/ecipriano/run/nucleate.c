/**
# Nucleate Boiling with Conjugate Heat Transfer

Simulation setup for nucleate boiling with conjugate heat transfer. This
setup is very similar to [bubblcecontact.c](../run/bubblcecontact.c) but
including the solution of the energy equation within the solid heater and the
energy transfer at the fluid--solid interface.

This setup aims to be the same described by [Torres et al., 2024](#torres2024)
which consists of the expansion of a FC-72 bubble expanding over a solid surface
in microgravity conditions. Experiments for this configuration were performed
in the context of the [RUBI project](https://www.eoportal.org/satellite-missions/iss-rubi)
which allowed boiling experiments to be carried out on the International Space
Station. The absence of gravity allows the bubble growth process to be observed
without the influence of buoyancy-driven flows which would promote fast
detachment of the bubble from the solid surface.

![Nucleate boiling experiment in microgravity conditions, showing a bubble growing in contact line regime over a solid surface, surrounded by several thermocouples.](https://www.eoportal.org/api/cms/documents/d/eoportal/rubi_auto4-jpeg){width="70%"}

The simulation shows that the Joule effect, localized on the top of the solid heater, increases
the solid temperature and it is responsible for the heating of the fluid environment. After a few
seconds, a small bubble is nucleated and it starts growing on the solid surface with the specified
contact angle. The physical consistency of the oscillations observed around the contact line should
be further investigated.

![Evolution of the interface and temperature field from the numerical simulation](nucleate/movie.mp4)(width="70%")
*/

#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "boiling.h"
#include "conjugate.h"
#include "view.h"

/**
## Boundary Conditions

Top and right sides are open boundaries, while no-slip conditions are imposed on
the left boundary which represents the solid wall. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

/**
A constant contact angle model is used. */

double theta0 = 27.;
vector h[];
h.t[left] = contact_angle (theta0*pi/180.);

/**
## Simulation setup

We initialize the maximum level of refinement, and the radius of the nucleus:
since we cannot simulate the nucleation process from the molecular scale, we
just create a very small bubble at a specific time instant. */

int maxlevel = 7;
double R0 = 0.25e-3;

int main (void) {

  /**
  We set properties of the technical fluid FC-72, in liquid (1) and in gas (2)
  phase. The solid (3) is sapphire. The solid heats up the fluid system through
  Joule effect, whose thermal power is set and included in the solution of the
  solid temperature equation. */

  rho1 = 1653.9, rho2 = 7.02, rho3 = 4890.;
  mu1 = 5.49e-4, mu2 = 1.228e-5;
  cp1 = 1069.2, cp2 = 2034., cp3 = 410.;
  lambda1 = 5.44e-2, lambda2 = 0.024, lambda3 = 10.9;
  dhev = 89.77e+3, qjoule = 1e4;

  /**
  The interface is set to the saturation temperature (considering P = 500 mbar),
  All the other initial temperatures are set to 0.5 K above the saturation. See
  [Torres et al., 2024](#torres2024) for a comprehensive discussion about it. */

  TIntVal = 310.85;
  TL0 = TIntVal + 0.5, TG0 = TL0, TS0 = TL0;

  /**
  We add the possibility to include an interfacial heat transfer resistance, computed
  using Hertz-Knudsen theory, and with a pre-factor `acc` which is a function of the
  accommodation coefficient. It is not the accommodation coefficient itself because
  different works have different definitions of such quantity. */

  RInt = sqrt (2.*pi*8.314*pow (TIntVal, 3.)/0.338)/rho2/sq (dhev);
  acc = 0;

  f.height = h;

  nv = 1;
  pcm.consistent = true;

  DT = 0.01;
  size (8e-3);
  for (maxlevel = 7; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  /**
  We initialize the solid heater geometry similarly to what we do with embed. */

  solid (cw, fw, -(x - X0) + 5e-3 + 1e-5);
  foreach()
    f[] = 1.;

  /**
  The boundary conditions for temperature are set here because the phase model
  initializes those in the `defaults` event. Both the solid and the fluid
  temperature boundary conditions are set to the temperature of the solid
  boundary and the temperatures of the fluid boundaries, respectively. Those
  values are computed in [conjugate.h](../src/conjugate.h). */

  scalar TL = liq->T, TG = gas->T;
  TS[left] = dirichlet (TSB[]);
  TL[left] = dirichlet (TLB[]);
  TG[left] = dirichlet (TGB[]);

  TL[right] = dirichlet (TL0);
}

/**
## Nucleation

After 5 seconds of simulation we create a small bubble. The surface tension is
set here in order to avoid to reduce the time step when the interface is not
present. */

event nucleation (t = 5) {
#if TREE
  refine (circle (x, y, 2.*R0)  > 0. && level < maxlevel);
#endif
  fraction (f, -circle (x - R0*cos(theta0*pi/180.), y, R0));
  f.sigma = 9.91e-3;

  scalar TL = liq->T, TG = gas->T;
  foreach() {
    TL[] *= f[];
    TG[] = (1. - f[])*TG0;
    T[] = TL[] + TG[];
  }
}

/**
## Adaptive grid

We adapt the simulation according to the interface, fluid temperature, and the
velocity field. Do we need to include also the solid temperature? I think that
it might be unnecessary since the fluid temperature will automatically refine
the region around the solid-fluid interface (the left boundary in this case).
This is convenient because we avoid introducing additional fluid cells due to
the solution of the solid phase on the same grid. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e-2,1e-2,1e-2,1e-2}, maxlevel, 5, 1);
}
#endif

/**
## Post-processing

We output the bubble equivalent diameter in time. */

event logger (t += 0.01) {
  if (t >= 5) {
    double volume = 0.;
    foreach (reduction(+:volume))
      volume += (1. - f[])*dv();
    double deq = cbrt (12.*volume);
    fprintf (stderr, "level %d %.4g %.4g %.4g\n", maxlevel, t, volume, deq),
            fflush (stderr);
  }
}

/**
We write a video with the evolution of the interface and the temperature
fields in the fluid and solid phases. */

#include "custom-cmaps.h"

trace
bool draw_line (coord o, coord d, float lc[3] = {0}, float lw = 1.) {
  bview * view = draw();
  draw_lines (view, lc, lw) {
    glBegin (GL_LINES);
    glvertex2d (view, o.x, o.y);
    glvertex2d (view, d.x, d.y);
    glEnd ();
  }
  return true;
}

event movie (t += 0.05; t <= 8) {
  clear();
  view (theta = 0., phi = 0., psi = -pi/2., ty = -0.3);
  squares ("T", spread = -1, map = RdBu);
  draw_vof ("f", lw = 2);
  draw_line ({0.,0.}, {0.,L0}, lw = 2.);
  mirror ({0,1}) {
    squares ("T", spread = -1, map = RdBu);
    draw_vof ("f", lw = 2);
    draw_line ({0.,0.}, {0.,L0}, lw = 2.);
  }
  mirror ({1,0}) {
    draw_vof ("cw", "fw", filled = -1, fc = {1.,1.,1.});
    draw_vof ("cw", "fw", lw = 2);
    squares ("TS", spread = -1, map = RdBu);
    mirror ({0,1}) {
      draw_vof ("cw", "fw", filled = -1, fc = {1.,1.,1.});
      draw_vof ("cw", "fw", lw = 2);
      squares ("TS", spread = -1, map = RdBu);
    }
  }
  save ("movie.mp4");
}

/**
## References

~~~bib
@article{torres2024,
  title={On the coupling between direct numerical simulation of nucleate boiling and a micro-region model at the contact line},
  author={Torres, Loric and Urbano, Annafederica and Colin, Catherine and Tanguy, S{\'e}bastien},
  journal={Journal of Computational Physics},
  volume={497},
  pages={112602},
  year={2024},
  publisher={Elsevier}
}
~~~
*/
