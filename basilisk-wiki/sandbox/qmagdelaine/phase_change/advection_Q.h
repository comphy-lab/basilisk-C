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
#define u uf

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
We apply boundary conditions after user initialisation. */

event init (i = 0) {
  boundary ((scalar *){u});
  boundary (tracers);
}

/**
The timestep is set using the velocity field and the CFL criterion. The
integration itself is performed in the events of [tracer.h]() and [vof.h]()
if it is included. 
An alias for stability is defined for backward compability: */

#define velocity stability

event stability (i++, last) {
  dt = dtnext (timestep (u, DT));
}

/**
We define empty events in order to fixe the order of the events: 
* 'stability()' where the velocity and the timestep are set,
* 'vof()' which advects a VOF tracer if [vof.h]() is included,
* 'tracer_advection()' and 'tracer_diffusion()' to advect (with [tracer.h]())
and make diffuse tracers,
* 'adapt()' to adapt the mesh. */

event vof (i++, last);
event tracer_advection (i++, last);
event tracer_diffusion (i++, last);

#if TREE
event adapt (i++,last);
#endif
#include "tracer.h"
