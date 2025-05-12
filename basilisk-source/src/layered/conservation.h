/**
# Conservation of the volume of each layer

This module enforces conservation of the volume of each layer using a
spatially-constant, but varying in time, "upwelling" velocity within
each layer which compensates any volume loss. This can be interpreted
as a vertical volume/mass transfer between layers. 

This upwelling velocity is only applied when the thickness of a layer
is larger than `hmin/10`. The value of `hmin` is typically that set in
[entrainment.h](). */

extern double * hmin;

/**
The domain considered needs to be simply-connected in order to avoid
the transfer of mass between basins which are not connected. To
enforce this we use the functions in [/src/tag.h](). */

#include "tag.h"

/**
This stores the initial volume of each layer. */

static struct {
  double * sum;
} Conservation = {0};

event init (i = 0)
{

  /**
  At initialisation we make sure that the domain is simply connected
  by retaining only the largest "oceanic" basin. We first define an
  indicator function which is zero on dry land and one in wet
  areas. */
  
  scalar d[];
  foreach() {
    double H = 0.;
    foreach_layer()
      H += h[];
    d[] = H > dry;
  }

  /**
  We remove cells which are only connected diagonally with their
  neighbors. */
  
  foreach()
    if ((d[] && !d[0,-1]) &&
	((d[-1,-1] && !d[-1]) || (d[1,-1] && !d[1])))
      d[] = 0.;

  /**
  We only keep the largest basin. */

  remove_droplets (d, -1);

  /**
  We "dry out" all the other basins. */
  
  foreach()
    if (!d[]) {
      foreach_layer()
	h[] = 0.;
      zb[] = HUGE;
    }

  /**
  And finally we store the initial volume of each layer. */

  Conservation.sum = malloc (nl*sizeof (double));
  foreach_layer()
    Conservation.sum[_layer] = statsf (h).sum;
}

/**
At each timestep... */

event viscous_term (i++)
{
  foreach_layer() {

    /**
    ... we compute the current volume of each layer and the area of
    cells which are "thick enough". */
    
    double area = 0., volume = 0.;
    foreach (reduction(+:area) reduction(+:volume)) {
      if (h[] > hmin[_layer]/10.)
	area += dv();
      volume += h[]*dv();
    }

    /**
    The ratio of these two quantities gives the vertical displacement,
    which is then applied in each cell which is "thick enough". */

    assert (area > 0.);
    double dh = (Conservation.sum[_layer] - volume)/area;
    assert (dh > - hmin[_layer]/10.);
    foreach()
      if (h[] > hmin[_layer]/10.)
	h[] += dh;
  }
}

/**
We free memory at the end of the simulation to avoid memory leaks. */

event cleanup (t = end)
{
  free (Conservation.sum);
}
