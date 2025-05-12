/**
# An advection solver

We wish to solve the advection equations
$$
\partial_tf_i+\mathbf{u}\cdot\nabla f_i=0
$$
where $\mathbf{u}$ is the velocity field and $f_i$ are a list of
passive tracers.  This can be done with a flux-based advection scheme
such as the 2nd-order, unsplit, upwind scheme of [Bell-Collela-Glaz,
1989](references.bib#bell89).

The main time loop is defined in [run.h](). A stable timestep needs to
respect the [CFL
condition](http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition). */

#include "run.h"
#include "timestep.h"

/**
We allocate the (face) velocity field. For compatibility with the
other solvers, we allocate it as `uf` and define an alias. The
`gradient` function is used to set the type of slope-limiting
required. The default is to not use any limiting (i.e. a purely
centered slope estimation). */

face vector uf[];
vector u[];
vector * u_ebit = {u};

double (* gradient) (double, double, double) = NULL;

/**
Here we set the gradient functions for each tracer (as defined in the
user-provided `tracers` list). */

extern scalar * tracers;

event defaults (i = 0) {
  for (scalar f in tracers)
    f.gradient = gradient;
}

/**
User initialisation happens here. */

event init (i = 0);

/**
The velocity event does nothing, since we use a constant timestep for
kinematice test cases.
The tracer.h is removed since it change the vof event order, resulting
in some issues when the EBIT is coupled with AMR. (patch 8dae9c096 on Nov. 9, 2024)*/

event velocity (i++, last) {}
