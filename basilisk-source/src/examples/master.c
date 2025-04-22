/**
# Coupling two solvers

In this example we couple different versions of the Navier-Stokes
solver: one with embedded boundaries and the other without, but the
technique is applicable to any combination of solvers, grids,
dimensions etc.

The goal of this example is to use the output of one solver (the
"[slave](slave.c)"), which generates a well-known [von Karman vortex
street](karman.c), as inflow boundary condition for another solver (the
"[master](master.c)") which simply solves the flow in a straight channel.

To do so, we run the two simulations simultaneously and we sample the
slave fields (which vary in space and time) along a vertical line in
the middle of the domain and use these values as inflow boundary
conditions for the left boundary of the master.

We generate movies (and other outputs) in the usual manner for the
master and the slave.

The result is illustrated in the two movies below.

![Animation of the vorticity field for the slave.](master/vort-slave.mp4)(autoplay)

![Animation of the vorticity field for the master.](master/vort-master.mp4)(autoplay)

Note that only aspects specific to the one-way coupling are documented
here. The master and slave simulations are very similar to the [von
Karman vortex street]() example which should be consulted for
complementary explanations.

Note also that coupling the two solvers requires some simple
compilation and linking tricks. See the `slave.o` target in the
[Makefile]() for details. */

#include "grid/multigrid.h"

/**
We must include the master driver *before* the solver. */

#include "master.h"

/**
We can autolink `slave.o` or link it manually. */

#pragma autolink slave.o

/**
Otherwise the simulation follows the setup described in [karman.c](). */

#include "navier-stokes/centered.h"

double Reynolds = 160.;
face vector muv[];

int main()
{
  L0 = 8. [1];
  origin (-L0/2., -L0/2.);
  N = 128;
  mu = muv;

  /**
  The slave simulation is advanced to $t=15$ before starting the
  master simulation, so that the von Karman vortex street is already
  developed. */
  
  slave_start = 15.;
  
  run(); 
}

double D = 0.5, U0 = 1.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Reynolds;
}

/**
## Boundary conditions

This is where the coupling happens: we set the left boundary condition
for the velocity on the master equal to the velocity components
interpolated along a vertical line located at $x =$ `-0.5 + L0/2.`
which corresponds to the middle of the slave simulation. */

u.n[left]  = dirichlet (slave_interpolate ("u.x", - 0.5 + L0/2., y, 0, true));
u.t[left]  = dirichlet (slave_interpolate ("u.y", - 0.5 + L0/2., y, 0, true));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

/**
The other boundary condition is just outflow. */

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
## Initialization, logfiles and outputs 

Nothing special here. */

event init (t = 0)
{
  foreach()
    u.x[] = U0;
}

event logfile (i++)
  fprintf (stderr, "master %d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);

event movies (t += 0.05; t <= 30.)
{
  scalar omega[];
  foreach()
    omega[] = (u.y[1] - u.y[-1] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
  output_ppm (omega, file = "vort-master.mp4", min = -5, max = 5, linear = true, n = 256);
}
