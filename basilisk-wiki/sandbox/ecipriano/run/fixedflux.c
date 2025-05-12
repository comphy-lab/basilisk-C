/**
# Fixed Flux Droplet Evaporation

Evaporation of a liquid droplet with a fixed vaporization
flowrate. In these conditions, the droplet consumption
is obtained from a mass balance on the liquid droplet and
the analytic solution can be easily derived:
$$
  \dfrac{dR}{dt} = \dfrac{\dot{m}}{\rho_l}
$$

We want to test and compare the droplet consumption with
the exact solution, considering the presence of the Stefan
flow. The simulation setup used here was adapted from
[Malan et al, 2021](#malan2021geometric).
The animation shows the consuption of the liquid droplet,
that maintains a perfectly spherical shape throughtout the
entire lifetime, as well as the radial velocity field
caused by the phase change.
![Evolution of the interface position and the velocity field](fixedflux/movie.mp4)(width="400" height="400")
*/

/**
## Phase Change Setup

We divide the mass balance by the density of the gas phase,
according with the non-dimensional treatment reported in
the [paper](#malan2021geometric). There is no need to shift
the volume expansion term because we do not transport the
temperature or mass fractions in this test case.
*/

#define BYRHOGAS
#define NOSHIFTING

/**
## Simulation Setup

We use the centered Navier-Stokes equations solver with
the evaporation source term in the projection step. The
double pressure velocity coupling approach is adopted to
obtain a divergence-free velocity which can be used for
the VOF advection equation. The evporation model is
combined with the fixed flux mechanism, that imposes
a constant vaporization flowrate.
*/

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "fixedflux.h"
#include "view.h"

/**
### Boundary Conditions

Outflow boundary conditions are imposed on every
domain boundary since the droplet is initially placed
at the center of the domain. Because the
centered-doubled method is used, the boundary conditions
must be imposed also for the *extended* fields.
*/

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
uext.n[bottom] = neumann (0.);
uext.t[bottom] = neumann (0.);
pext[bottom] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

/**
### Model Data

We set the initial radius of the droplet, and
the value of the vaporization rate.
*/

int maxlevel;
double R0 = 0.23;
double mEvapVal = -0.05;

int main(void) {
  L0 = 1.;

  /**
  We set the density and viscosity values. */

  rho1 = 2., rho2 = 1.;
  mu1 = 1., mu2 = 0.1;

  /**
  The surface tension is set to zero in order to assess
  the sphericity of the droplet without being influenced
  by the surface tension. */

  f.sigma = 0.;

  /**
  We make lists with the levels of refinement and the
  corresponding maximum time steps. We need to impose the
  time steps because the surface tension is set to zero. */

  double mllist[] = {4, 5, 6, 7, 8};
  double dtlist[] = {0.005, 0.002, 0.001, 0.0005, 0.0001};

  /**
  We modify the Poisson equation solver tolerance and
  we create a folder where the facet outputs will be
  written. */

  system ("mkdir facets");
  TOLERANCE = 1.e-4;

  /**
  We perform the simulation at different levels of
  refinement. */

  for (int sim = 0; sim < 3; sim++) {
    maxlevel = mllist[sim];
    DT = dtlist[sim];
    TOLERANCE = 1.e-4;
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x - 0.5*L0) - sq(y - 0.5*L0))

/**
We initialize the volume fraction field.
*/

event init (i = 0) {
  fraction (f, circle (x, y, R0));
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the analytic
solution for the droplet volume.
*/

double exact (double time) {
  return pi*sq (R0 + mEvapVal*time);
}

/**
### Output Files

We write the droplet volume from the numerical simulation
as well as the exact solution.
*/

event output (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  double dropvol = statsf(f).sum;
  double relerr = fabs (exact(t) - dropvol) / exact(t);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%.12f %.12f %.12f %.12f\n", t, dropvol, exact(t), relerr);
  fflush (fp);
}

/**
### Facets

We write the vof facets in such a way that it works also
in parallel: every processor writes its own facets in
a different file and the files are then gathered in a
single facets output file, ready to be plotted.
*/

event output_facets (t += 1) {
  char names[80];
  sprintf (names, "interfaces%d", pid());
  FILE * fpf = fopen (names,"w");
  output_facets (f, fpf);
  fclose (fpf); fflush (fpf);
  char command[80];
  sprintf(command, "LC_ALL=C cat interfa* > facets/facets-%d-%.1f", maxlevel, t);
  system(command);
}

/**
### Movie

We write the animation with the volume fraction field
and the velocity vectors.
*/

event movie (t += 0.1; t <= 4) {
  if (maxlevel == 6) {
    clear();
    view (tx = -0.5, ty = -0.5);
    draw_vof ("f", lw = 1.5);
    vectors ("u", scale = 7.e-2);
    box();
    save ("movie.mp4");
  }
}

/**
## Results

The droplet is consumed by the evaporation is a smooth way,
and the comparison between the numerical and the analytic
droplet show a very good agreement and the convergence to
the exact solution.

~~~gnuplot Droplet Volume
reset
set xlabel "t [s]"
set ylabel "Droplet Volume [m^3]"
set size square
set grid

plot "OutputData-6" every 100 u 1:3 w p ps 2 t "Analytic", \
     "OutputData-4" u 1:2 w l lw 2 t "LEVEL 4", \
     "OutputData-5" u 1:2 w l lw 2 t "LEVEL 5", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-4" using 4 nooutput name "LEVEL4"
stats "OutputData-5" using 4 nooutput name "LEVEL5"
stats "OutputData-6" using 4 nooutput name "LEVEL6"

set print "Errors.csv"

print sprintf ("%d %.12f", 2**4, LEVEL4_mean)
print sprintf ("%d %.12f", 2**5, LEVEL5_mean)
print sprintf ("%d %.12f", 2**6, LEVEL6_mean)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[8:128]
set size square
set grid

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
  10*x**(-1) lw 2 title "1^{st} order", \
  30*x**(-2) lw 2 title "2^{nd} order"
~~~

~~~gnuplot Droplet Sphericity
reset
set size square
set xrange[0.5:0.75]
set yrange[0.5:0.75]
set xlabel "radius [m]"
set ylabel "radius [m]"
set grid

array r[3]
r[1] = 0.18
r[2] = 0.13
r[3] = 0.08

set style fill transparent solid 0.2 noborder

set object 1 circle front at 0.5,0.5 size r[1] fillcolor rgb "black" lw 1
set object 2 circle front at 0.5,0.5 size r[2] fillcolor rgb "black" lw 1
set object 3 circle front at 0.5,0.5 size r[3] fillcolor rgb "black" lw 1

p \
  "facets/facets-4-1.0" w l lw 2 lc 1 t "LEVEL 4", \
  "facets/facets-4-2.0" w l lw 2 lc 1 notitle, \
  "facets/facets-4-3.0" w l lw 2 lc 1 notitle, \
  "facets/facets-5-1.0" w l lw 2 lc 2 t "LEVEL 5", \
  "facets/facets-5-2.0" w l lw 2 lc 2 notitle, \
  "facets/facets-5-3.0" w l lw 2 lc 2 notitle, \
  "facets/facets-6-1.0" w l lw 2 lc 3 t "LEVEL 6", \
  "facets/facets-6-2.0" w l lw 2 lc 3 notitle, \
  "facets/facets-6-3.0" w l lw 2 lc 3 notitle
~~~

## References

~~~bib
@article{malan2021geometric,
  title={A geometric VOF method for interface resolved phase change and conservative thermal energy advection},
  author={Malan, LC and Malan, Arnaud G and Zaleski, St{\'e}phane and Rousseau, PG},
  journal={Journal of Computational Physics},
  volume={426},
  pages={109920},
  year={2021},
  publisher={Elsevier}
}
~~~
*/
