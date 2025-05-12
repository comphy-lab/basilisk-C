/**
# Conservation of $h - z_b$ when reconstructing $a$ 

When using the default adaptive reconstruction of variables, the
[blood flow solveer](bloodflow-hr.h) will conserve the cross-sectional
area $a$ when cells are refined or coarsened. However, this will not
necessarily ensure that the "lake-at-rest" condition (i.e. $\eta = h -
z_b = \mathrm{cst}$) is also preserved. In what follows, we redefine
the *prolongation()* and *restriction()* methods of the
cross-sectional area $a$ so that $\eta$ is conserved.

We start with the reconstruction of fine "wet" cells: */

#if TREE
static void refine_elevation (Point point, scalar a)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  double eta = k[]*sqrt (a[]) - zb[]; // water surface elevation  
  coord g; // gradient of eta
  if (gradient)
    foreach_dimension()
      g.x = gradient (k[-1]*sqrt (a[-1]) - zb[-1], eta, k[1]*sqrt (a[1]) - zb[1])/4.;
  else
    foreach_dimension()
      g.x = (k[1]*sqrt (a[1]) - zb[1] - k[-1]*sqrt (a[-1]) + zb[-1])/(2.*Delta);
  // reconstruct water depth h from eta and zb

  foreach_child() {
    double etac = eta;
    foreach_dimension()
      etac += g.x*child.x;
    a[] = sq (1./k[]*max (0., etac + zb[]));
  }
}

/**
Cell restriction is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */

static void restriction_elevation (Point point, scalar a)
{
  double eta = 0., v = 0.;
  foreach_child() {
      eta += a[]*(k[]*sqrt (a[]) - zb[]);
      v += a[];
  }

  /**
  ... and use this in combination with $z_b$ (of the coarse cell) to
  compute the water depth $h$.  */
    
  a[] = sq (1./k[]*max (0., eta/v + zb[]));
}

/**
We also need to define a consistent prolongation function. For cells
which are entirely surrounded by wet cells, we can use the standard
linear refinement function, otherwise we use straight injection from
the parent cell. */

static void prolongation_elevation (Point point, scalar a)
{
  bool wet = true;
  foreach_neighbor(1)
    if (a[] <= dry) {
      wet = false;
      break;
    }
  if (wet)
    refine_linear (point, a);
  else {
    double ac = a[], zc = zb[];
    foreach_child() {
      a[] = ac;
      zb[] = zc;
    }
  }
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  a.refine  = refine_elevation;
  a.prolongation = prolongation_elevation;
  a.restriction = restriction_elevation;
  a.dirty = true;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif
