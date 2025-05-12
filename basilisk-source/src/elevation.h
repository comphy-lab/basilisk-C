/**
# Conservation of water surface elevation 

When using the default adaptive reconstruction of variables, the
[Saint-Venant solver](saint-venant.h) or the [layered solver](hydro.h)
will conserve the water depth when cells are refined or
coarsened. However, this will not necessarily ensure that the
"lake-at-rest" condition (i.e. a constant water surface elevation) is
also preserved. In what follows, we redefine the *prolongation()* and
*restriction()* methods of the water depth $h$ so that the water
surface elevation $\eta$ is conserved.

We start with the reconstruction of fine "wet" cells: */

#if TREE
static double default_sea_level = 0.;

static void refine_elevation (Point point, scalar h)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  if (h[] >= dry) {
    double eta = zb[] + h[];   // water surface elevation  
    coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
	g.x = gradient (zb[-1] + h[-1], eta, zb[1] + h[1])/4.;
    else
      foreach_dimension()
	g.x = (zb[1] - zb[-1])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta;
      foreach_dimension()
	etac += g.x*child.x;
      h[] = max(0., etac - zb[]);
    }
  }
  else {

    /**
    The "dry" case is a bit more complicated. We look in a 3x3
    neighborhood of the coarse parent cell and compute a depth-weighted
    average of the "wet" surface elevation $\eta$. We need to do this
    because we cannot assume a priori that the surrounding wet cells are
    necessarily close to e.g. $\eta = 0$. */

    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood
    foreach_neighbor(1)
      if (h[] >= dry) {
	eta += h[]*(zb[] + h[]);
	v += h[];
      }
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else

      /**
      If none of the surrounding cells is wet, we set a default sealevel. */

      eta = default_sea_level;

    /**
    We then reconstruct the water depth in each child using $\eta$ (of the
    parent cell i.e. a first-order interpolation in contrast to the wet
    case above) and $z_b$ of the child cells. */
    
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0., eta - zb[]);
  }
}

/**
Cell restriction is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */

static void restriction_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }

  /**
  ... and use this in combination with $z_b$ (of the coarse cell) to
  compute the water depth $h$.  */
    
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

/**
We also need to define a consistent prolongation function. For cells
which are entirely surrounded by wet cells, we can use the standard
linear refinement function, otherwise we use straight injection from
the parent cell. */

static void prolongation_elevation (Point point, scalar h)
{
  bool wet = true;
  foreach_neighbor(1)
    if (h[] <= dry) {
      wet = false;
      break;
    }
  if (wet)
    refine_linear (point, h);
  else {
    double hc = h[], zc = zb[];
    foreach_child() {
      h[] = hc;
      zb[] = zc;
    }
  }
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h.refine  = refine_elevation;
  h.prolongation = prolongation_elevation;
  h.restriction = restriction_elevation;
  h.dirty = true;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif
