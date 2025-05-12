/**
# Stefan Problem

A thin vapor layer is initialized close to a superheated
wall. The high wall temperature heats up the vapor layer,
which is in contact with a liquid phase that always remains
at saturation temperature leading to the phase change. The
simulation setup used here was adapted from [Malan et al, 2021](#malan2021geometric).
The analytic solution to the problem describes the evolution
of the vapor layer thickness:
$$
  \delta(t) = 2\lambda\sqrt{\alpha_g t}
$$
where $\alpha_g$ is the thermal diffusivity, defined as
$\alpha_g = \lambda_g/\rho_g/cp_g$, while $\lambda$ is
a value, specific to the physical properties under
investigation, which can be found from the solution
of the transcendental equation:
$$
  \lambda\exp(\lambda^2) \text{erf}(\lambda) =
  \dfrac{cp_g(T_{wall}-T_{sat})}{\Delta h_{ev}\sqrt{\pi}}
$$
The animation shows that, at the beginning of the simulation,
the superheated wall heats up the vapor layer which becomes
hotter than the liquid phase, leading to the evaporation
which establishes the interface velocity jump.
![Evolution of the interface and the temperature field](stefanproblem/movie.mp4)(width="400" height="400")
*/

/**
## Phase Change Setup

We move the interface using the velocity *uf*, with the
expansion term shifted toward the gas-phase. In this way
*uf* is divergence-free at the interface. The double
pressure velocity couping is used to obtain an extended
velocity, used to transport the gas phase tracers. */

#define INT_USE_UF
#define CONSISTENTPHASE2
#define SHIFT_TO_GAS

/**
## Simulation Setup

We use the centered solver with the divergence source term,
and the extended velocity is obtained using the centered-doubled
procedure. The evaporation model is used in combination
with the temperature-gradient mechanism, which manages
the solution of the temperature field.
*/

#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "two-phase.h"
#include "tension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "view.h"

/**
We declare the variables required by the
temperature-gradient model. */

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Twall;

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are
imposed on the left wall, while the temperature on the
superheated wall is set to the superheated value.
*/

u.n[left] = neumann (0);
u.t[left] = neumann (0);
p[left] = dirichlet (0);
uext.n[left] = neumann (0);
uext.t[left] = neumann (0);
pext[left] = dirichlet (0);

TG[left]  = dirichlet (Tsat);
TG[right] = dirichlet (Twall);
TL[left]  = dirichlet (Tsat);
TL[right] = dirichlet (Twall);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, and the initial thickness of
the vapor layer.
*/

int maxlevel, minlevel = 3;
double lambdaval = 0.06779249298045148;
double  delta0 = 322.5e-6;
double tshift, teff;


int main (void) {
  /**
  We set the material properties of the
  fluids. */

  rho1 = 958., rho2 = 0.6;
  mu1 = 2.82e-4, mu2 = 1.23e-5;
  lambda1 = 0.68, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2080.;
  dhev = 2.256e6,

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tsat = 373., Twall = 383.;
  TL0 = Tsat, TG0 = Tsat; TIntVal = Tsat;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 10e-3;
  f.sigma = 0.059;

  /**
  We run the simulation for different maximum
  levels of refinement. */

  for (maxlevel = 5; maxlevel <= 6; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
We initialize the volume fraction field and the temperature
in the gas and in liquid phase. */

event init (i = 0) {
  fraction (f, -(x - 0.0096775));

  foreach() {
    TL[] = f[]*TL0;
    TG[] = (1. - f[])*TG0;
    T[] = f[]*TL0 + (1. - f[])*TG0;
  }

  /**
  At simulation time equal to zero, the thickness of the
  vapor layer is not zero. Therefore, we compute a time
  *shifting factor* (just for post-processing purposes). */

  double effective_height;
  effective_height = (sq(L0) - statsf(f).sum)/L0;
  tshift = sq(effective_height/2./lambdaval)*rho2*cp2/lambda2;
}

/**
We refine the interface and the region where the
temperature field changes. */

#if TREE
event adapt (i++) {
  adapt_wavelet_leave_interface ({T}, {f},
      (double[]){1.e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the exact solution to the
thickness of the vapor layer, and the analytic temperature profile.
*/

double exact (double time) {
  return 2.*lambdaval*sqrt(lambda2/rho2/cp2*time);
}

double tempsol (double time, double x) {
  return Twall + ((Tsat - Twall)/erf(lambdaval))*
    erf(x/2./sqrt(lambda2/rho2/cp2*time));
}

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

event output (t += 0.1) {
  double effective_height = 0.;
  foreach(reduction(+:effective_height))
    effective_height += (1. - f[])*dv();
  effective_height /= L0;

  double relerr = fabs (exact(t+tshift) - effective_height) / exact(t+tshift);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%g %g %g %g\n", t+tshift, effective_height, exact (t+tshift), relerr);
  fflush (fp);
}

/**
### Temperature Profile

We write on a file the temperature profile at the
final time step. */

event profiles (t = end) {
  char name[80];
  sprintf (name, "Temperature-%d", maxlevel);

  FILE * fpp = fopen (name, "w");
  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    double R = exact (t+tshift);
    double tempexact = (x >= L0-R) ? tempsol (t+tshift, L0-x) : Tsat;
    fprintf (fpp, "%g %g %g\n", x, interpolate (T, x, 0.), tempexact);
  }
  fflush (fpp);
  fclose (fpp);
}

/**
### Movie

We write the animation with the evolution of the
temperature field and the gas-liquid interface.
*/

event movie (t += 0.1; t <= 10) {
  if (maxlevel == 5) {
    clear();
    view (tx = -0.5, ty = -0.5);
    box();
    draw_vof ("f", lw = 1.5);
    squares ("T", min = Tsat, max = Twall, linear = true);
    vectors ("u", scale = 1.e-1, lc = {1.,1.,1.});
    save ("movie.mp4");
  }
}

/**
## Results

~~~gnuplot Evolution of the vapor layer thickness
set xlabel "t[s]"
set ylabel "Vapor Layer Thickness [m]"
set key left top
set size square
set grid

plot "OutputData-5" every 10 u 1:3 w p ps 2 t "Analytic", \
     "OutputData-5" u 1:2 w l lw 2 t "LEVEL 5", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-4" using (last4=$4) nooutput
stats "OutputData-5" using (last5=$4) nooutput
stats "OutputData-6" using (last6=$4) nooutput
#stats "OutputData-7" using (last7=$4) nooutput

set print "Errors.csv"

print sprintf ("%d %.12f", 2**5, last5)
print sprintf ("%d %.12f", 2**6, last6)
#print sprintf ("%d %.12f", 2**7, last7)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**4:2**7]
set yr[1e-6:1e-2]

set size square
set grid

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
  0.05*x**(-1) lw 2 title "1^{st} order", \
  0.5*x**(-2)  lw 2 title "2^{nd} order"
~~~

~~~gnuplot Temperature Profile
reset
set xlabel "Length [m]"
set ylabel "Temperature [K]"

set xr[0.0075:0.01]

set key bottom right
set size square
set grid

plot "Temperature-5" u 1:3 w p ps 2 t "Analytic", \
     "Temperature-5" u 1:2 w l lw 2 t "LEVEL 5", \
     "Temperature-6" u 1:2 w l lw 2 t "LEVEL 6"
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
