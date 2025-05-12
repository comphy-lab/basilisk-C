/**
# Sucking Interface Problem

Planar boiling configuration. A thin gas layer is initialized on the
left side of the 2D domain. The gas is at saturation temperature and
it is in contact with a superheated liquid phase. The temperature gradient
between the interface, at saturation temperature, and the surrounding
superheated environment leads to the phase change. This boiling
problem is similar to the [Scriven test case](/sandbox/ecipriano/run/scrivenproblem.c),
but with a planar geometry. This test case setup was inspired by
[Welch and Wilson, 2000](#welch2000volume), and [Zhao et al., 2022](#zhao2022boiling).
The analytical solution describing the gas layer thickness in time
reads:

$$
  \delta(t) = 2\beta\sqrt{\alpha_g t}
$$

where $\alpha_g$ is the thermal diffusivity, defined as
$\alpha_g = \lambda_g/\rho_g/cp_g$, while $\beta$ is the
growth constant, which is obtained as a function of the
physical properties of the simulation, as reported in
[Boyd et al., 2023](#boyd2023consistent):
$$
\exp (\beta^2) \text{erf} (\beta)
\left(
  \beta - \dfrac{\left(T_{bulk} - T_{sat}\right)c_{p,g}\lambda_l
  \sqrt{\alpha_g} \exp\left(-\beta^2
  \dfrac{\rho_g^2\alpha_g}{\rho_l^2\alpha_l}\right)}{\Delta h_{ev}
  \lambda_g \sqrt{\pi \alpha_l}\text{erfc}\left(
  \beta\dfrac{\rho_g\sqrt{\alpha_g}}{\rho_l\sqrt{\alpha_l}} \right)}
\right)
= 0
$$

The animation shows the map of the temperature field, which
is initialized with the analytic solution at the beginning
of the simulation. The velocity field due to the expansion term in
the continuity equation points toward the liquid phase. Therefore,
this test case considers the additional transport of temperature from
the Stefan flow in liquid phase.

![Evolution of the gas layer thickness](suckingproblem/movie.mp4)(width="400" height="400")
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

Outflow boundary conditions for velocity and pressure are imposed on
the right wall, while the temperature is set to the bulk value on
the wall in contact with the liquid phase, and to the saturation
value on the wall adjacent to the gaseous layer. */

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

TL[right] = dirichlet (Tbulk);
TG[left] = dirichlet (Tsat);

/**
### Problem Data

We declare the maximum and minimum levels of refinement,
the $\lambda$ parameter, the initial thickness of the vapor layer. */

int maxlevel, minlevel = 1;
double tshift, betaGrowth, effective_height0;

/**
### Growth Constant and Analytical Temperature

We use gsl, and the [fsolve.h](/sandbox/ecipriano/src/fsolve.h)
module, to find the value of the growth constant, which is then
used to compute the analytical temperature profile:

$$
  T(x,t) = T_{bulk} -
  \left(\dfrac{T_{bulk}-T_{sat}}{\text{erfc}\left( \beta
  \dfrac{\rho_g\sqrt{\alpha_g}}{\rho_l\sqrt{\alpha_l}} \right)} \right)
  \text{erfc}\left( \dfrac{x}{2\sqrt{\alpha_l t}} +
  \beta\dfrac{\rho_g -
  \rho_l}{\rho_l}\sqrt{\dfrac{\alpha_g}{\alpha_k}} \right)
$$
*/

#define USE_GSL
#include "fsolve-gsl.h"

/**
We define the function to zero for the growth constant. */

int betafun (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  double beta = xdata[0];
  double alpha1 = lambda1/rho1/cp1;
  double alpha2 = lambda2/rho2/cp2;

  fdata[0] = exp(sq(beta))*erf(beta)*(beta -
          ( (Tbulk - Tsat)*cp2*lambda1*sqrt (alpha2)*exp
          (-sq(beta)*sq(rho2)*alpha2/sq(rho1)/alpha1) ) /
          (dhev*lambda2*sqrt(pi*alpha1)*
          erfc(beta*rho2*sqrt(alpha2)/rho1/sqrt(alpha1))));

  return GSL_SUCCESS;
}

/**
We define the function for the exact temperature field. */

double tempexact (double x, double beta, double t) {
  double alpha1 = lambda1/rho1/cp1;
  double alpha2 = lambda2/rho2/cp2;

  return Tbulk - ((Tbulk - Tsat)/erfc(beta*rho2*sqrt(alpha2)/rho1/sqrt(alpha1)))
      * erfc (x/2./sqrt(alpha1*t) + beta*(rho2 - rho1)/rho1*sqrt(alpha2/alpha1));
}

int main (void) {
  /**
  We set the material properties of the fluids. */

  rho1 = 958.4, rho2 = 0.597;
  mu1 = 2.80e-4, mu2 = 1.26e-5;
  lambda1 = 0.679, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2030.;
  dhev = 2.26e+6;

  /**
  We set the bulk and the saturation temperature values. The
  interface must remain at saturation temperature, while the values
  `TL0` an `TG0` are used to set the initial and boundary conditions.
  */

  Tbulk = 378.15, Tsat = 373.15, TIntVal = Tsat;
  TL0 = Tbulk, TG0 = TIntVal;

  /**
  We change the dimension of the domain, and the surface tension
  coefficient. */

  L0 = 10.e-3;
  f.sigma = 0.0059;

  /**
  We run the simulation at four different levels of refinement. */

  for (maxlevel = 6; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

event init (i = 0) {

  /**
  The gas phase layer thickness is initially set to 0.476 mm. */

  fraction (f, (x - 0.476e-3));
  effective_height0 = L0 - 1./L0*statsf(f).sum;

  /**
  We use the root finding procedure implemented in gsl_multiroots
  for solving the transcendental equation to find the growth
  constant. */

  Array * arrUnk = array_new();
  {
    double betafg = 0.9;
    array_append (arrUnk, &betafg, sizeof(double));
    double * unks = (double *)arrUnk->p;
    fsolve (betafun, arrUnk, NULL);
    betaGrowth = unks[0];
  }
  array_free (arrUnk);

  /**
  We find the time of the analytical solution corresponding to the
  value of the vapor layer thickness that we just initialized. */

  tshift = rho2*cp2/lambda2*sq (effective_height0/2./betaGrowth);

  /**
  We initialize the temperature field, setting the liquid phase
  temperature to the analytical value. We also set an initial value
  for the velocity field to speed up convergence of the Poisson
  equation at the first iteration. */

  foreach() {
    TL[] = f[]*tempexact (x, betaGrowth, t+tshift);
    TG[] = (1. - f[])*Tsat;
    T[]  = TL[] + TG[];

    u.x[] = 11.e-3*f[];
    //u.x[] = 0.;
    u.y[] = 0.;
  }
}

/**
We refine the domain according to the interface, velocity, and the
temperature field. */

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e-3,uemax,uemax}, maxlevel, minlevel, 1);
}
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the exact solution for the
thickness of the gas phase layer in time. */

double exact (double time) {
  return 2.*betaGrowth*sqrt(lambda2/rho2/cp2*time);
}

/**
### Output Files

We write the thickness of the gas layer from this simulation, the
exact solution, and the relative error. */

event output (t += 0.01) {
  double effective_height = L0 - 1./L0*statsf(f).sum;

  double relerr = fabs (exact(t+tshift) - effective_height) / exact(t+tshift);

  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  // Average velocity
  double uavg = 0., favg = 0.;
  foreach(reduction(+:uavg) reduction(+:favg)) {
    double ux = 0.5*(uf.x[] + uf.x[1]);
    uavg += f[]*ux;
    favg += f[];
  }
  uavg /= favg;

  fprintf (fp, "%g %g %g %g %g\n", t + tshift, effective_height, exact (t+tshift), relerr, uavg);
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
  for (double x = 0; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    double tempanal = tempexact (x, betaGrowth, t+tshift);
    double radius = exact (t+tshift);
    tempanal = (x <= radius) ? Tsat : tempanal;
    fprintf (fpp, "%g %g %g\n", x, interpolate (T, x, 0.5*L0), tempanal);
  }
  fflush (fpp);
  fclose (fpp);
}

/**
### Movie

We write the animation with the evolution of the temperature field,
and the gas-liquid interface. We reduce the total simulation time
in order to speed up the simulation on the basilisk server. */

//event movie (t += 0.01; t <= 1) {
event movie (t += 0.01; t <= 0.2) {
  clear();
  view (tx = -0.5, ty = -0.5);
  box();
  draw_vof ("f", lw = 1.5);
  squares ("T", min = Tsat, max = Tbulk, linear = true);
  save ("movie.mp4");
}

/**
## Results

~~~gnuplot Evolution of the vapor layer thickness
reset
set xlabel "t [s]"
set ylabel "Vapor Layer Thickness [m]"
set key top left
set size square
set grid

plot "OutputData-8" every 1 u 1:3 w p ps 2 t "Analytic", \
     "OutputData-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "OutputData-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "OutputData-8" u 1:2 w l lw 2 t "LEVEL 8"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-6" using (last6=$4) nooutput
stats "OutputData-7" using (last7=$4) nooutput
stats "OutputData-8" using (last8=$4) nooutput

set print "Errors.csv"

print sprintf ("%d %.12f", 2**6, last6)
print sprintf ("%d %.12f", 2**7, last7)
print sprintf ("%d %.12f", 2**8, last8)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[2**5:2**9]
set yr[1e-3:1]

set size square
set grid

plot "Errors.csv" w p pt 8 ps 2 title "Results", \
  30*x**(-1) lw 2 title "1^{st} order", \
  600*x**(-2) lw 2 title "2^{nd} order"
~~~

~~~gnuplot Temperature Profile
reset
set xlabel "Radius [m]"
set ylabel "Temperature [K]"

set key bottom right
set size square
set grid
set yr[373:379]

plot "Temperature-8" every 5 u 1:3 w p ps 2 t "Analytic", \
     "Temperature-6" u 1:2 w l lw 2 t "LEVEL 6", \
     "Temperature-7" u 1:2 w l lw 2 t "LEVEL 7", \
     "Temperature-8" u 1:2 w l lw 2 t "LEVEL 8"
~~~

## References

~~~bib
@article{welch2000volume,
  title={A volume of fluid based method for fluid flows with phase change},
  author={Welch, Samuel WJ and Wilson, John},
  journal={Journal of computational physics},
  volume={160},
  number={2},
  pages={662--682},
  year={2000},
  publisher={Elsevier}
}

@article{zhao2022boiling,
  title={Boiling and evaporation model for liquid-gas flows: A sharp and conservative method based on the geometrical VOF approach},
  author={Zhao, Shuo and Zhang, Jie and Ni, Ming-Jiu},
  journal={Journal of Computational Physics},
  volume={452},
  pages={110908},
  year={2022},
  publisher={Elsevier}
}

@article{boyd2023consistent,
  title={A consistent volume-of-fluid approach for direct numerical simulation of the aerodynamic breakup of a vaporizing drop},
  author={Boyd, Bradley and Ling, Yue},
  journal={Computers \& Fluids},
  volume={254},
  pages={105807},
  year={2023},
  publisher={Elsevier}
}
~~~
*/
