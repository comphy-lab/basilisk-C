/**
# Fick Corrected Approach

The solution of diffusion equations for chemical species mass
fractions or concentrations usually exploit Fick's law. Such law
is valid for a binary mixture, where the diffusivity values are
the same for both chemical species. In that case, the diffusion of a
chemical species is balanced by the counter-diffusion of the
other species. Therefore, the closure of the diffusive fluxes is
respected:

$$
  \sum_{i=1}^{NS} j_i^F = 0
$$

When mixtures with more than two chemical species are solved, Fick's
law does not guarantee the closure of the diffusive fluxes anymore.
This happens because, usually, the diffusivity value is not the same
for every chemical species. To overcome this problem, different
approaches where proposed. The Maxwell-Stefan model is most elegant
approach for the solution of diffusive fluxes in multicomponent
mixtures. However, this approach is computationally expensive, and
not suitable for the segregated solution of the transport equations.
A simpler approach that we want to test here is the *Fick Corrected*
model, which defines the diffusive flux for each chemical species $i$
as:

$$
  j_i = j_i^F - \omega_i \sum_{j=1}^{NS} j_j^F
$$

where $j_i^F = -\rho\mathcal{D}_i\nabla\omega_i$ is the classic Fick's law,
$\omega_j$ is the mass fraction of the chemical species $i$, and the
RHS is the error on the diffusive fluxes, re-distributed on every
species on the basis of the mass fraction. */

#include "grid/multigrid.h"
#include "diffusion.h"
#include "run.h"

/**
We define an iterator to simplify loops over lists of chemical
species. */

#define foreach_elem(list, index) \
for (int index=0; index<list_len (list); index++)

/**
The diffusivity value is constant and uniform but it changes
depending on the species. */

attribute {
  double Dmix;
}

/**
We create 3 chemical species mass fraction fields (`Y1`, `Y2`, `Y3`)
and a scalar field that collects these species. In the same way, we
define the diffusive flux for each species, the total diffusive flux
and a list that collects the different fluxes. */

scalar Y1[], Y2[], Y3[];
scalar * YList = {Y1, Y2, Y3};

scalar J1[], J2[], J3[], Jtot[];
scalar * JList = {J1, J2, J3};

/**
## Problem setup

We set dirichlet boundary conditions everywhere on the left and right
boundaries. */

Y1[left] = dirichlet (1.);
Y2[left] = dirichlet (0.);
Y3[left] = dirichlet (0.);

Y1[right] = dirichlet (0.);
Y2[right] = dirichlet (1.);
Y3[right] = dirichlet (0.);

/**
The boolean value `fick_corrected` allows different `run()` with or
without the fick corrected approach, in order to see the difference.
*/

bool fick_corrected = true;

int main (void) {

  /**
  We set the diffusivity of each chemical species choosing three
  different values at three different orders of magnitude. */

  Y1.Dmix = 1.e-2;
  Y2.Dmix = 1.e-3;
  Y3.Dmix = 1.e-4;

  /**
  We set the CFL value, create the grid, and run the simulation. */

  CFL = 0.5;
  init_grid (1 << 7);
  run();
  fick_corrected = false;
  run();
}

/**
We set the initial conditions on the mass fractions. */

event init (i = 0) {
  foreach() {
    Y1[] = 0.;
    Y2[] = 0.;
    Y3[] = 1.;
  }
}

/**
## Stability Conditions

The correction is added as an explicit source term, therefore a
stability condition on the diffusion process is computed. */

event stability (i++) {
  double Dmix = -HUGE;
  for (scalar Y in YList)
    Dmix = max (Dmix, Y.Dmix);

  double dmin = HUGE;
  foreach_face (reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (Delta < dmin) dmin = Delta;
    }
  double dtmax = sq (dmin)/4./Dmix;
  dt = dtnext (CFL*dtmax);
}

/**
## Diffusion Equation

We solve the diffusion equation including the correction to the
diffusive fluxes. */

event tracer_diffusion (i++) {

  /**
  We compute the diffusive flux for each chemical species. */

  foreach_elem (YList, j) {
    scalar Y = YList[j];
    scalar J = JList[j];

    face vector JJ[];
    foreach_face()
      JJ.x[] = Y.Dmix*face_gradient_x (Y, 0);

    foreach() {
      J[] = 0.;
      foreach_dimension()
        J[] += (JJ.x[1] - JJ.x[])/Delta*cm[];
    }
  }

  /**
  We compute the total diffusive flux. */

  foreach() {
    Jtot[] = 0.;
    if (fick_corrected)
      for (scalar J in JList)
        Jtot[] -= J[];
  }

  /**
  We use the values in the vector `JList` to store the diffusion
  correction. */

  foreach() {
    foreach_elem (YList, j) {
      scalar Y = YList[j];
      scalar J = JList[j];

      J[] = Jtot[]*Y[];
    }
  }

  /**
  We call the diffusion equations solver. */

  foreach_elem (YList, j) {
    scalar Y = YList[j];
    scalar J = JList[j];

    face vector Dmixf[];
    foreach_face()
      Dmixf.x[] = Y.Dmix;
    diffusion (Y, dt, Dmixf, r=J);
  }
}

/**
## Post-Processing

We plot the chemical species mass fraction profiles along the length
of the domain for each chemical species. */

event profiles (t = 2.) {
  scalar Ytot[];
  foreach() {
    Ytot[] = 0.;
    for (scalar Y in YList)
      Ytot[] += Y[];
  }

  for (double x=0; x<L0; x += 0.5*L0/(1 << grid->maxdepth))
    fprintf (stderr, "%g %g %g %g %g\n", x,
        interpolate (Y1, x, 0.5*L0), interpolate (Y2, x, 0.5*L0),
        interpolate (Y3, x, 0.5*L0), interpolate (Ytot, x, 0.5*L0));
  fprintf (stderr, "\n\n");
}

/**
## Results

We can see that the sum of the mass fractions without the Fick
corrected approach give unphysical values grater than 1. The
correction modifies the diffusive fluxes in order to enforce this
constraint.

~~~gnuplot Evolution of mass fraction with classic Fick's law
reset
set xlabel "length [m]"
set ylabel "Mass Fractions [-]"
set yrange[-0.2:2]
set grid

plot "log" index 1 u 1:2 w l lw 2 t "Species 1", \
     "log" index 1 u 1:3 w l lw 2 t "Species 2", \
     "log" index 1 u 1:4 w l lw 2 t "Species 3", \
     "log" index 1 u 1:5 w l lw 2 t "Total"
~~~

~~~gnuplot Evolution of mass fraction with corrected Fick's law
reset
set xlabel "length [m]"
set ylabel "Mass Fractions [-]"
set yrange[-0.2:1.4]
set grid

plot "log" index 0 u 1:2 w l lw 2 t "Species 1", \
     "log" index 0 u 1:3 w l lw 2 t "Species 2", \
     "log" index 0 u 1:4 w l lw 2 t "Species 3", \
     "log" index 0 u 1:5 w l lw 2 t "Total"
~~~
*/

