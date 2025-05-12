/**
# Lid-driven cavity at Re=1000

We use the multigrid implementation (rather than the default tree
implementation) and either the MAC or the centered Navier--Stokes
solver. */

#include "grid/cartesian.h"


/**
Here we define the domain geometry: a square box of size unity
centered on (0,0). We also set the viscosity and some parameters 
controlling the numerical scheme. */

int main()
{ 
  // coordinates of lower-left corner
  origin (-0.5, -0.5, -0.5);
  // number of grid points
  init_grid (64);
  // viscosity

  const face vector muc[] = {1e-3,1e-3,1e-3};
  mu = muc;

  // maximum timestep
  DT = 0.1;
  // CFL number
  CFL = 0.8;

  /**
  We then call the `run()` method of the Navier--Stokes solver. */

  run();
}

/**
The default boundary conditions are symmetry (i.e. slip walls). We
need no-slip on three boundaries and $u=1$ on the top
boundary i.e. */

u.t[top] = dirichlet(1);
u.r[top] = dirichlet(0);



u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
u.t[front] = dirichlet(0);
u.t[back]   = dirichlet(0);
/**
For the colocated solver, imposing boundary conditions for the normal
components of the (face-centered) advection velocity improves the
results. Ideally, this should be done automatically by the solver. */



/**
We define an auxilliary function which computes the total kinetic
energy. The function works both for face and centered
discretisations of the velocity. */



/**
We add an option to restore the simulation from a previous dump. */



/**
We want the simulation to stop when we are close to steady state. To
do this we store the `u.x` field of the previous timestep in an
auxilliary variable `un`. */

scalar un[];

/**
Every 0.1 time units we check whether $u$ has changed by more than
10^-5^. If it has not, the event returns 1 which stops the
simulation. We also output the evolution of the kinetic energy on
standard error. */

event logfile (t += 0.1; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-5)
    return 1; /* stop */
  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

