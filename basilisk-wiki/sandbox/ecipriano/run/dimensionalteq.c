/**
# Dimensional Analysis for a Temperature Equation

This is a test to familiarize with the dimensional analysis, applied to the solution of a simple temperature diffusion equation:
$$
\rho c_p \dfrac{\partial T}{\partial t} =
\nabla\cdot\left(\lambda\nabla T\right)
$$
In principle, the solver should only implement the diffusion step,
getting the dimensions of $\rho$, $c_p$, and $\lambda$, from the user.
[This file](dimensionalteq/dimensionalteq.dims) contains the output from the dimensional analysis.

![Evolution of the temperature](dimensionalteq/movie.mp4)
*/

#include "diffusion.h"
#include "run.h"
#include "view.h"

/**
For convience, we define a convention on the principal dimensions,
and other macros that define common dimensions.
*/

// Base dimensions
#define SI_Length  [1]
#define SI_Time    [0,1]
#define SI_Mass    [0,0,1]
#define SI_Temp    [0,0,0,1]
#define SI_Moles   [0,0,0,0,1]
#define SI_Current [0,0,0,0,0,1]
#define SI_Lumin   [0,0,0,0,0,0,1]
// Common dimensions
#define SI_Density [-3,0,1]
#define SI_Joule [2,-2,1]
#define SI_Watt [2,-3,1]
#define SI_SpecificHeat [2,-2,0,-1]
#define SI_ThermalCond [1,-3,1,-1]
#define SI_MassFrac [0]
#define SI_MoleFrac [0]
#define SI_Diffusivity [2,-1]

/**
We define the inputs necessary for solving the diffusion equation.
Problably it is not necessary to set the units for all these
quantities, but it seems easier to me to just set the units of
every input value. */

double T0 = 273. SI_Temp;
double TL = 373. SI_Temp;
double TR = 173. SI_Temp;
double rho = 1.  SI_Density;
double cp = 1.   SI_SpecificHeat;
double k = 1.e-2 SI_ThermalCond;

/**
We create the temperature field and the thermal diffusivity,
and we set the boundary conditions for temperature. */

scalar T[];
face vector lambda[];

T[left] = dirichlet (TL);
T[right] = dirichlet (TR);

int main (void) {
  DT = 0.1 SI_Time;
  init_grid (1 << 7);
  run();
}

/**
We set the initial value of the temperature field, and and the thermal diffusivity. */

event init (i = 0) {
  foreach()
    T[] = T0;
  
  foreach_face()
    lambda.x[] = k/rho/cp;
}

/**
The timestep is imposed to be equal to `DT`. */

event stability (i++) {
  dt = dtnext (DT);
}

/**
We solve the diffusion equation. The unit of the temperature should be automatically deduced by qcc. */

event tracer_diffusion (i++) {
  diffusion (T, dt, D=lambda);
}

/**
We print the dimensions of some constants used in the code and of the
temperature field. Outside the `foreach` loop the dimensions of `T[]` are
not printed unless a `Point` is defined. */

event logfile (t = end) {
  show_dimension (rho);
  show_dimension (lambda);
  show_dimension (k);

  foreach()
    show_dimension (T[]);
}

/**
Movie with the evolution of the temperature field.
*/

event movie (t += 0.1; t <= 10) {
  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("T", min = TR, max = TL,
      linear = true);
  save ("movie.mp4");
}

