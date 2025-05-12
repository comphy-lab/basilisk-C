/**

Adapted from the Hele-Shaw solver in src. */

#include "advection_reaction.h"
#include "poisson.h"

/**
We allocate the pressure $p$ and divergence field $\zeta$. The
$\beta$ coefficients need to be defined at face locations to
match the locations of the face pressure gradients (and the
face velocity components). These two sets of coefficients are
stored in a vector field. We also allocate space to store the
statistics of the Poisson solver. */

scalar p[], zeta[], m[], opening[], ones[];
face vector beta[];
mgstats mgp;


/**
At every timestep, but after all the other events for this timestep
have been processed (the '`last`' keyword), we update the pressure
field $p$ by solving the Poisson equation with variable coefficient
$\beta$. */

event pressure (i++, last)
{
  mgp = poisson (p, zeta, beta);
  
  /**
  We then update the velocity field by computing the face pressure
  gradients. */
  
  foreach_face()
    u.x[] = beta.x[]*(p[] - p[-1])/Delta;
  

}
