/**
# Scriven Problem

Expansion of a bubble in a superheated environment.
A spherical bubble is initialized on the bottom of
an AXI domain. The bubble is at saturation temperature
while the liquid environment is superheated. The
phase change pheomena, promoted by the temperature
gradient between the surface of the bubble and the liquid
environment leads to the expansion of the bubble.
The simulation setup used here was inspired
by [Tanguy et al., 2014](#tanguy2014benchmarks), and
is characterized by $Ja=3$. The analytic solution to
this problem was proposed by [Scriven, 1959](#scriven1959dynamics),
and it describes the dynamic of the bubble radius in time:
$$
  R(t) = 2\beta\sqrt{\alpha_l t}
$$
where $\alpha_g$ is the thermal diffusivity, defined as
$\alpha_l = \lambda_l/\rho_l/cp_l$, while $\beta$ is the
growth constant, which is obtained as a function of the
physical properties of the simulation:
$$
  \dfrac{\rho_l cp_l (T_{bulk}-T_{sat})}{\rho_g
  (\Delta h_{ev} + (cp_l - cp_g)(T_{bulk} - T_{sat}))} =
  2\beta^2 \int_{0}^1
  \exp \left(-\beta^2((1-x)^{-2}
  -2(1-\dfrac{\rho_l}{\rho_g})x - 1)\right)dx
$$

The animation shows the map of the temperature field, which
is initialized with the analytic solution at the beginning
of the simulation. We let the bubble expand until reaching
a radius which is twice the initial radius. The Stefan flow
contribution must be considered because it provides an
additional transport of the temperature field in liquid
phase. The analytic solution takes into account this
phenomena.

![Evolution of the interface and temperature field](scrivenproblem/movie.mp4)(width="400" height="400")
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
#define INIT_TEMP

/**
## Simulation Setup

We use the centered solver with the divergence source term,
and the extended velocity is obtained using the centered-doubled
procedure. The evaporation model is used in combination
with the temperature-gradient mechanism, which manages
the solution of the temperature field. */

#include "axi.h"
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
double TL0, TG0, TIntVal, Tsat, Tbulk;

/**
### Boundary conditions

Outflow boundary conditions for velocity and pressure are
imposed on the top, left, and right walls. The temperature
on these boundaries is imposed to the bulk value. */

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

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

TL[top] = dirichlet (Tbulk);
TL[right] = dirichlet (Tbulk);
TL[left] = dirichlet (Tbulk);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, the initial radius of the droplet,
the growth constant, and additional post-processing
variables. */

int maxlevel, minlevel = 5;
double R0 = 1.e-3;
double XC, YC;
double V0, tshift;
double betaGrowth = 3.32615013;

/**
### Initial Temperature

We use gsl to compute the integral required to
obtain the analytic temperature profile:
$$
  T(r,t) = T_{bulk} - 2\beta^2
  \left(
  \dfrac{\rho_g (\Delta h_{ev} +
  (cp_l-cp_g)(T_{bulk} - T_{sat}))}{\rho_l cp_l}
  \right)
  \int_{1-R/r}^1
  \exp \left(-\beta^2((1-x)^{-2}
  -2(1-\dfrac{\rho_l}{\rho_g})x - 1)\right)dx
$$
*/

#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas

double intfun (double x, void * params) {
  double beta = *(double *) params;
  return exp(-sq(beta)*(pow(1. - x, -2.) - 2.*(1. - rho2/rho1)*x - 1 ));
}

double tempsol (double r, double R) {
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error;
  double beta = betaGrowth;
  gsl_function F;
  F.function = &intfun;
  F.params = &beta;
  gsl_integration_qags (&F, 1.-R/r, 1., 1.e-9, 1.e-5, 1000,
                        w, &result, &error);
  gsl_integration_workspace_free (w);
  return Tbulk - 2.*sq(beta)*(rho2*(dhev + (cp1 - cp2)*(Tbulk - Tsat))/rho1/cp1)*result;
}

int main (void) {
  /**
  We set the material properties of the fluids. */

  rho1 = 958.; rho2 = 0.59;
  mu1 = 2.82e-4; mu2 = 1.23e-6;
  lambda1 = 0.6, lambda2 = 0.026;
  cp1 = 4216., cp2 = 2034.;
  dhev = 2.257e+6;

  /**
  The initial bubble temperature and the interface
  temperature are set to the saturation value. */

  Tbulk = 373.989, Tsat = 373., TIntVal = 373.;
  TL0 = Tbulk, TG0 = TIntVal;

  /**
  We change the dimension of the domain,
  the surface tension coefficient, and the coordinates
  of the center of the bubble. */

  L0 = 12.e-3;
  XC = 0.5*L0, YC = 0.;
  f.sigma = 0.001;

  /**
  We reduce the tolerance of the Poisson
  equation solver. */

  TOLERANCE = 1.e-6;

  /**
  We compute the time shifting factor (post-processing),
  since the bubble Radius at simulation time t=0 is not
  zero. */

  double alpha = lambda1/rho1/cp1;
  tshift = sq(R0/2./betaGrowth)/alpha;

  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */

  double dtlist[] = {1.e-4, 5.e-5, 5.e-5, 5.e-5};
  int mllist[] = {6, 7, 8, 9};

  /**
  We run the simulation for different
  levels of refinement. */

  for (int sim=0; sim<4; sim++) {
    DT = dtlist[sim];
    maxlevel = mllist[sim];
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x, y, R) (sq(R) - sq(x - XC) - sq(y - YC))

event init (i = 0) {
  fraction (f, circle(x,y,R0));

  foreach()
    f[] = 1. - f[];

  /**
  We initialize the temperature field. The liquid
  phase temperature is set to the analytic value. */

  foreach() {
    double r = sqrt( sq(x - XC) + sq(y - YC) );
#ifdef INIT_TEMP
      TL[] = f[]*tempsol (r, R0);
#else
    TL[] = f[]*TL0;
#endif
    TG[] = (1. - f[])*TG0;
    T[]  = TL[] + TG[];
  }
}

/**
We refine the domain according to the interface and the
temperature field. */

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({T,u}, {f},
      (double[]){1e-3,uemax,uemax,uemax,1e-3}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the exact
solution to the thickness of the vapor layer.
*/

double exact (double time) {
  return 2.*betaGrowth*sqrt(lambda1/rho1/cp1*time);
}

/**
### Output Files

We write the bubble radius and the analytic solution on 
a file. */

event output (t += 0.005) {
  scalar fg[];
  foreach()
    fg[] = 1. - f[];

  double effective_radius = pow(3./2.*statsf(fg).sum, 1./3.);

  double rsol = exact (t+tshift);
  double relerr = (rsol > 0.) ? fabs (rsol - effective_radius) / rsol : 0.;

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g %g\n",
      t+tshift, effective_radius, exact (t+tshift), relerr);
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
  for (double x = 0.5*L0; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    double r = x - 0.5*L0;
    double R = exact (t+tshift);
    double tempexact = (r >= R) ? tempsol (r, R) : Tsat;
    fprintf (fpp, "%g %g %g\n", x, interpolate (T, x, 0.), tempexact);
  }
  fflush (fpp);
  fclose (fpp);
}

/**
### Movie

We write the animation with the evolution of the
temperature field and the gas-liquid interface. */

event movie (t += 0.01; t <= 0.5) {
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Tbulk);
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Evolution of the bubble radius
reset
set xlabel "t [s]"
set ylabel "Bubble Radius [m]"
set key top left
set size square
set grid

plot "OutputData-6" every 8 u 1:3 w p ps 2 t "Analytic", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "OutputData-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "OutputData-8" u 1:2 w l lw 2 t "LEVEL 8", \
     "OutputData-9" u 1:2 w l lw 2 t "LEVEL 9"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-6" using (last6=$4) nooutput
stats "OutputData-7" using (last7=$4) nooutput
stats "OutputData-8" using (last8=$4) nooutput
stats "OutputData-9" using (last9=$4) nooutput

#stats "OutputData-6" using 4 nooutput name "LEVEL6"
#stats "OutputData-7" using 4 nooutput name "LEVEL7"
#stats "OutputData-8" using 4 nooutput name "LEVEL8"
#stats "OutputData-9" using 4 nooutput name "LEVEL9"

set print "Errors.csv"

#print sprintf ("%d %.12f", 2**6, LEVEL6_mean)
#print sprintf ("%d %.12f", 2**7, LEVEL7_mean)
#print sprintf ("%d %.12f", 2**8, LEVEL8_mean)
#print sprintf ("%d %.12f", 2**9, LEVEL9_mean)

print sprintf ("%d %.12f", 2**6, last6)
print sprintf ("%d %.12f", 2**7, last7)
print sprintf ("%d %.12f", 2**8, last8)
print sprintf ("%d %.12f", 2**9, last9)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**5:2**10]
set yr[1e-4:10]

set size square
set grid

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
  50*x**(-1) lw 2 title "1^{st} order", \
  1000*x**(-2) lw 2 title "2^{nd} order"
~~~

~~~gnuplot Temperature Profile
reset
set xlabel "Radius [m]"
set ylabel "Temperature [K]"

set yr[372.9:374.1]

set key bottom right
set size square
set grid

plot "Temperature-6" u 1:3 w p ps 2 t "Analytic", \
     "Temperature-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "Temperature-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "Temperature-8" u 1:2 w l lw 2 t "LEVEL 8", \
     "Temperature-9" u 1:2 w l lw 2 t "LEVEL 9"
~~~

## References

~~~bib
@article{tanguy2014benchmarks,
  title={Benchmarks and numerical methods for the simulation of boiling flows},
  author={Tanguy, S{\'e}bastien and Sagan, Micha{\"e}l and Lalanne, Benjamin and Couderc, Fr{\'e}d{\'e}ric and Colin, Catherine},
  journal={Journal of Computational Physics},
  volume={264},
  pages={1--22},
  year={2014},
  publisher={Elsevier}
}

@article{scriven1959dynamics,
  title={On the dynamics of phase growth},
  author={Scriven, LE},
  journal={Chemical engineering science},
  volume={10},
  number={1-2},
  pages={1--13},
  year={1959},
  publisher={Elsevier}
}
~~~
*/
