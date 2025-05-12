/**
# Interface Regression Velocity

This module test the capability of the interface regression
model implemented in [evaporation.h](../src/evaporation.h).
The interface gradients are computed at cell centeres, and
so is the vaporization rate. Therefore, interpolations are
needed to transform the vaporization rate in an interface
regression velocity, defined on cell faces, and that allows
the interface to be transported by a modified velocity.
This velocity is not conservative, but the "amount" of non
conservation reflects the amount of reference phase to be
added or removed by the phase change phenomena.

**PROS**: This approach works well also when the velocity
field is zero everywhere, and the only contribution to the
interface motion is the phase change (i.e., isolated droplet
with zero surface tension). While these conditions are ideal,
they can be exploited to decouple the effect of phase change
from other physical phenomena. Apart from these numerical
conditions, this method limits the possibility of under-
and over-shoots, that easily arise when applying the phase
change term in the VOF transport equation as an explicit
source. In that case, the use of a redistribution algorithm
is crucial.

**CONS**: The interface regression velocity is distributed
in a manner that guarantees the correct consumption of the
reference phase (according to the material balance). This
does not mean that the divergence of the velocity used is
equal to the vaporization rate, localized at the interface.
Therefore, it cannot be used for the transport of tracers
associated with the vof field, using an advection equation
in conservative form (e.g., conserved tracers in the all-mach
solver). The transport of tracers in non-conservative form,
as currently implemented in [vof.h](/src/vof.h), works fine.

## Redistribution Procedure

The vaporization rate, in each interfacial cell, is
distributed among the faces in contact with that interfacial
cell. There are two possible cases:

1. If the face connects an interfacial cell with a pure cell,
the interface regression velocity on that face is taken from
the interfacial cell, weighted by the corresponding normal
component.

2. If the face connects two interfacial cells, the interface
regression velocity is computed from a linear interpolation
between the two consecutive vaporization rates, weighted by
the interface normal.

![Interface Regression Velocity Distribution](interfaceregression.png){ width="80%"}
*/

/**
## Phase Change Setup

We suppress the expansion term in the continuity equation,
it is not relevant for this test. Therefore, we do not need
a method to compute the *extended* velocity. For this
reason, after the declaration of the field $\mathbf{u}_f$
we set the extended velocity to be equal to this one-field
face velocity. */

#define DIFFUSIVE

#include "grid/multigrid.h"
#include "navier-stokes/centered-evaporation.h"
#define ufext uf
#include "two-phase.h"
#include "evaporation.h"
#include "fixedflux.h"
#include "view.h"

/**
## Simulation Setup

We set the value of the vaporization rate per unit of
interface surface. We also declare the index of the
simulation case (we run 3 different cases). */

double mEvapVal = -0.02;
int sim = 0;

int main (void) {
  origin (-0.5, -0.5);
  DT = 1.e-2;
  for (sim=0; sim<3; sim++) {
    init_grid (1 << 6);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))
#define iplane(x,y)(x-y+1.e-5)
#define ellipse(x,y,R)(1 - sq(x/(1.2*R)) - sq(y/(0.8*R)))

event init (i = 0) {
  switch (sim) {
    case 0: fraction (f, circle(x,y,0.23)); break;
    case 1: fraction (f, iplane(x,y)); break;
    case 2: fraction (f, ellipse(x,y,0.3)); break;
  }
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Compute Divergence

We compute the divergence of the modified velocity field,
used by the vof advection. The interface regression velocity
is computed in the vof event. Therefore, we compute the
divergence in the tracer_advection event, which is executed
right after vof. */

scalar divu[];

event tracer_advection (i++) {
  foreach() {
    divu[] = 0.;
    foreach_dimension()
      divu[] += ((uf.x[1] - vpc.x[1]) - (uf.x[] - vpc.x[]));
    divu[] /= Delta;
  }
}

/**
### Time-Derivative of the Volume Fraction

We want to check that, even if we distribute the vaporization
rate, we obtain a variation of the volume fraction which is
coherent with the cell-centered vaporization rate:

$$
  \dfrac{\partial f}{\partial t} =
  \dfrac{f^{t+1} - f^t}{\Delta t} =
  \dfrac{\dot{m}}{\rho_l}
$$

To compare the two quantities, we compute the volume
integrals of the vaporization rate and of the vof fraction
time derivative.

First, we store the old time before the vof advection. */

scalar fold[];

event vof (i++) {
  foreach()
    fold[] = f[];
}

/**
Then, we compute the volume integrals. */

scalar mEvapEff[], mEvapVol[];

event tracer_advection (i++) {
  foreach() {
    mEvapEff[] = (f[] - fold[])/dt;
    mEvapEff[] = fabs (mEvapEff[]) > F_ERR ? mEvapEff[] : 0.;
    vofrecon vr = vof_reconstruction (point, f);
    mEvapVol[] = mEvapTot[]*vr.dirac;
  }
  double totmEvapEff = rho1*statsf(mEvapEff).sum;
  double totmEvapVol = rho1*statsf(mEvapVol).sum;

  fprintf (stderr, "%g %g %g\n", t, totmEvapEff, totmEvapVol);
}

event closefile (t = end,last) {
  fprintf (stderr, "\n\n");
}

/**
### Movie

We write the animation with the divergence of the velocity
field and the interface facets. */

void write_movie (char * name) {
  clear();
  draw_vof ("f", lw = 3);
  squares ("divu", spread=-1);
  cells();
  save (name);
}

event movie (t += 0.1; t <= 10) {
  switch (sim) {
    case 0: write_movie ("case1.mp4"); break;
    case 1: write_movie ("case2.mp4"); break;
    case 2: write_movie ("case3.mp4"); break;
  }
}

/**
## Results

The animations show the divergence of the velocity that
transports the volume fraction.

![Sphere](interfaceregression/case1.mp4)

![Plane](interfaceregression/case2.mp4)

![Ellipse](interfaceregression/case3.mp4)

The following plot shows the comparison between the variation
of liquid volume according to the vaporization rate, and the
actual variation from the interpolation presented in this
module.

~~~gnuplot Comparison between time derivative and vaporization rate
reset
set xlabel "time [s]"
set ylabel "Vaporization Rate [kg/s]"
set size square
set key bottom right
set grid

plot "log" index 0 u 1:2 w l lw 2 t "Time Derivative Sphere", \
     "log" index 0 u 1:3 w l lw 2 t "Vaporization Rate Sphere", \
     "log" index 1 u 1:2 w l lw 2 t "Time Derivative Plane", \
     "log" index 1 u 1:3 w l lw 2 t "Vaporization Rate Plane", \
     "log" index 2 u 1:2 w l lw 2 t "Time Derivative Ellipse", \
     "log" index 2 u 1:3 w l lw 2 t "Vaporization Rate Ellipse"
~~~
*/
