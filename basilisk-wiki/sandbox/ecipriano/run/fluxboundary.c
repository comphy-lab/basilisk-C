/**
In this test case, a non-diffusive tracer, representing a mass-fraction field,
is advected by the velocity field. The goal is to compute the mass balance:
$$
  \dfrac{dm}{dt} = \dot{m}_{in} - \dot{m}_{out}
$$
by computing the accumulation term, and by integrating the flux of tracer
through the boundaries of the domain.
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tracer.h"
#include "view.h"

/**
We define the tracer field, and the list of tracers
to be advected using by the [advection.h](/src/advection.h) functions. */

scalar tr[];
scalar * tracers = {tr};

/**
We define the maximum level of refinemet, the inlet
velocity, and the total simulation time. */

int sim;
int maxlevel = 7;
double uin = 1.e-3;
double tEnd = 2.;
double mass_in_domain, mass_in_domain_old;

/**
The tracer is set to 0 everywhere at the beginning
of the simulation, while it is injected from the
bottom of the domain. */

f[left]  = dirichlet ( (y <= 0.5e-3) ? 1. : 0. );
tr[left] = dirichlet ( (y <= 0.5e-3) ? 1. : 0. );
u.n[left] = dirichlet (uin); u.t[left] = dirichlet (0.);
u.n[right] = neumann (0.); u.t[right] = neumann (0.);
p[left] = neumann (0.); p[right] = dirichlet (0.);

int main (void) {

  /**
  We define the physical properties for the simulation. */

  rho1 = 1000., rho2 = 1000.;
  mu1 = 1.e-3, mu2 = 1.e-3;

  f.sigma = 0.;
  tr.gradient = minmod2;

  /**
  If tracers is empty, use f.tracers,
  while if f.tracers is empty, tracers must be used. */

  tracers = NULL;
  f.tracers = {tr};

  L0 = 1.e-3;
  origin (0., 0.);
  size (L0);

  double dtlist[] = {1.e-2, 1.e-3, 1.e-4};
  for (sim = 0; sim < 3; sim++) {
    DT = dtlist[sim];
    init_grid (1 << maxlevel);
    run();
  }
}

event init (i = 0) {

  /**
  We set the initial conditions for the vof-field,
  and for the velocity, in order to avoid a sudden
  jump in the velocity field at time = 0. */

  foreach() {
    u.x[] = uin;
    u.y[] = 0.;
    f[] = (y <= 0.5e-3) ? 1. : 0.;
  }
  boundary({f});
  mass_in_domain = rho1 * statsf(tr).sum;
}

event compute_fluxes (i++) {
  mass_in_domain_old = mass_in_domain;
}

/**
# The output_data event computes the mass balance:
* *mass_inlet* is inlet mass flowrate of tracer $\dot{m}_{in}$ [kg/s]
* *mass_outlet* is the outlet mass flowrate of tracer $\dot{m}_{out}$ [kg/s]
* *delta_mass_over_dt* is the accumulation term, computed as the finite
difference between the total mass in the domain at time $t+1$ and at time $t$
*/

event output_data (i++,last) {

  /**
  We compute the total mass in the domain. */

  //double mass_inlet = rho1*uin*0.5*L0;
  mass_in_domain = rho1 * statsf(tr).sum;
  double delta_mass_over_dt = (mass_in_domain - mass_in_domain_old) / dt;

  /**
  We exploit the tracer_fluxes function, in [bcg.h](/src/bcg.h), even when
  the tracer is advected using the vof fluxes, because it easier, but it
  can be improved. */

  // Outlet mass flowrate
  face vector trfluxes[];
  tracer_fluxes (tr, uf, trfluxes, dt, zeroc);

  double mass_outlet = 0.;
  foreach_boundary (right, reduction(+:mass_outlet)) {
    mass_outlet += rho1*trfluxes.x[1]*Delta;
  }
  double mass_inlet = 0.;
  foreach_boundary (left, reduction(+:mass_inlet)) {
    mass_inlet += rho1*trfluxes.x[0]*Delta;
  }
  double rhs = delta_mass_over_dt + mass_outlet;

  char name[80];
  sprintf (name, "OutputData-%d-%d", maxlevel, sim);
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%f %f %f %f %f\n", t, mass_inlet, mass_outlet, delta_mass_over_dt, rhs);
  fflush (fp);
}

event movie (t += 0.1) {
  clear ();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("tr", min = 0., max = 1., linear = false);
  save ("t.mp4");
}

event stop (t = tEnd) {
}

/**
# Results
![Evolution of the tracer field](fluxboundary/t.mp4)(width="800" height="600")
~~~gnuplot Plot of the mass balance at three different maximum time steps.
set xlabel "time [s]"
set ylabel "mass flowrate [kg/s]"

p "OutputData-7-0" every 20:20 u 1:2 w p title "inlet mass flowrate", \
  "OutputData-7-0" u 1:4 w l title "dmdt", \
  "OutputData-7-0" u 1:3 w l title "outlet mass flowrate", \
  "OutputData-7-0" u 1:5 w l title "dmdt + outlet mass flowrate at DT = 10^{-2}", \
  "OutputData-7-1" u 1:5 w l title "dmdt + outlet mass flowrate at DT = 10^{-3}", \
  "OutputData-7-2" u 1:5 w l title "dmdt + outlet mass flowrate at DT = 10^{-4}"
~~~
There is a small discrepancy between the accumulation term and the outlet mass flowrate,
which becomes smaller by decreasing the time step.
*/