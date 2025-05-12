#include "vof.h"
#include "tracer.h"

scalar d[], f[], * interfaces = {f}, * tracers = {d};

#include "two-phase-generic.h"

/**
The initial volume fraction is computed from the initial distance
field, which must be initialised by the user. */

event init (i = 0)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
  fractions (phi, f);
}

/**
The distance function is reinitialised at each timestep. */

#include "height_distance.h"

event properties (i++)
{
  vector hf[];
  heights(f,hf);

  height_distance (hf, d, 1.); // Only using distance calculated from HF
}
