double s_clean = 1.e-10;

scalar cs2[];
face vector fs2[];

event metric (i = 0)
{
  foreach()
    cs2[] = 1.;
  foreach_face()
    fs2.x[] = 1.;
#if TREE
  cs2.refine = embed_fraction_refine;

  /**
  For prolongation we cannot use the same function since the surface
  fraction field *fs2* is not necessarily defined for prolongation
  cells. So we switch back to the default fraction refinement (which
  is less accurate but only relies on *cs2*). */

  cs2.prolongation = fraction_refine;
  foreach_dimension()
    fs2.x.prolongation = embed_face_fraction_refine_x;
  
  /**
  Note that we do not need to change the `refine` method since the
  default `refine` method calls the prolongation method for each
  component. */
  
#endif
  boundary ({cs2, fs2});
  restriction ({cs2, fs2});
  // cm2 = cs2; // needs to be added in common.h & tree-common.h
  // fm2 = fs2;
  
}



static inline void refine_embed_linear2 (Point point, scalar s)
{
  foreach_child() {
    if (!cs2[])
      s[] = 0.;
    else {
      assert (coarse(cs2));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fs2.x,i) && coarse(fs2.y,0,j) &&
      (coarse(cs2) == 1. || coarse(cs2,child.x) == 1. ||
       coarse(cs2,0,child.y) == 1. || coarse(cs2,child.x,child.y) == 1.)) {
    assert (coarse(cs2,child.x) && coarse(cs2,0,child.y));
    if (coarse(fs2.x,i,child.y) && coarse(fs2.y,child.x,j)) {
      // bilinear interpolation
      assert (coarse(cs2,child.x,child.y));
      s[] = (9.*coarse(s) + 
         3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
         coarse(s,child.x,child.y))/16.;
    }
    else
      // triangular interpolation     
      s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(cs2,child.x,child.y) &&
           ((coarse(fs2.x,i) && coarse(fs2.y,child.x,j)) ||
        (coarse(fs2.y,0,j) && coarse(fs2.x,i,child.y)))) {
    // diagonal interpolation
    s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fs2.x,i) > 0.25 && coarse(fs2.y,0,j) > 0.25 &&
      coarse(fs2.z,0,0,k) > 0.25 &&
      (coarse(cs2) == 1. || coarse(cs2,child.x) == 1. ||
       coarse(cs2,0,child.y) == 1. || coarse(cs2,child.x,child.y) == 1. ||
       coarse(cs2,0,0,child.z) == 1. || coarse(cs2,child.x,0,child.z) == 1. ||
       coarse(cs2,0,child.y,child.z) == 1. ||
       coarse(cs2,child.x,child.y,child.z) == 1.)) {
    assert (coarse(cs2,child.x) && coarse(cs2,0,child.y) &&
        coarse(cs2,0,0,child.z));
    if (coarse(fs2.x,i,child.y) && coarse(fs2.y,child.x,j) &&
        coarse(fs2.z,child.x,child.y,k) &&
        coarse(fs2.z,child.x,0,k) && coarse(fs2.z,0,child.y,k)) {
      assert (coarse(cs2,child.x,child.y) && coarse(cs2,child.x,0,child.z) &&
          coarse(cs2,0,child.y,child.z) &&
          coarse(cs2,child.x,child.y,child.z));
      // bilinear interpolation
      s[] = (27.*coarse(s) + 
         9.*(coarse(s,child.x) + coarse(s,0,child.y) +
             coarse(s,0,0,child.z)) + 
         3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
             coarse(s,0,child.y,child.z)) + 
         coarse(s,child.x,child.y,child.z))/64.;
    }
    else
      // tetrahedral interpolation
      s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
         coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(cs2,child.x,child.y,child.z) &&
           ((coarse(fs2.z,child.x,child.y,k) &&
         ((coarse(fs2.x,i) && coarse(fs2.y,child.x,j)) ||
          (coarse(fs2.y,0,j) && coarse(fs2.x,i,child.y))))
        ||
        (coarse(fs2.z,0,0,k) &&
         ((coarse(fs2.x,i,0,child.z) && coarse(fs2.y,child.x,j,child.z)) ||
          (coarse(fs2.y,0,j,child.z) && coarse(fs2.x,i,child.y,child.z))))
        ||
        (coarse(fs2.z,child.x,0,k) &&
         coarse(fs2.x,i) && coarse(fs2.y,child.x,j,child.z))
        ||
        (coarse(fs2.z,0,child.y,k) &&
         coarse(fs2.y,0,j) && coarse(fs2.x,i,child.y,child.z))
        ))
    // diagonal interpolation
    s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
    // Pathological cases, use 1D gradients.
    s[] = coarse(s);
    foreach_dimension() {
      if (coarse(fs2.x,(child.x + 1)/2) && coarse(cs2,child.x))
        s[] += (coarse(s,child.x) - coarse(s))/4.;
      else if (coarse(fs2.x,(- child.x + 1)/2) && coarse(cs2,- child.x))
        s[] -= (coarse(s,- child.x) - coarse(s))/4.;
    }
      }
    }
  }
}

/**
## Surface force and vorticity

We first define a function which computes
$\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ while taking the boundary
conditions on the embedded surface into account. */

attribute {
  void (* embed_gradient2) (Point, scalar, coord *);
}

static inline
coord embed_gradient2 (Point point, vector u, coord p, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs2, n, p, vb, &val);
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

/**
## Restriction of cell-centered fields

We now define restriction and prolongation functions for cell-centered
fields. The goal is to define second-order operators which do not use
any values from cells entirely contained within the embedded boundary
(for which *cs = 0*). 

When restricting it is unfortunately not always possible to obtain a
second-order interpolation. This happens when the parent cell does not
contain enough child cells not entirely contained within the embedded
boundary. In these cases, some external information (i.e. a boundary
gradient condition) is required to be able to maintain second-order
accuracy. This information can be passed by defining the
*embed_gradient()* function of the field being restricted. */

static inline void restriction_embed_linear2 (Point point, scalar s)
{  
  // 0 children
  if (!cs2[]) {
    s[] = 0.;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero embedded fractions), we return the
  corresponding value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if (fine(cs2,0,i,j) && fine(cs2,1,!i,!j))
  val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are defined
  (there is at least one). */
  
  coord p = {0.,0.,0.};
  foreach_child()
    if (cs2[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;

  /**
  If the gradient is defined and if the variable is not using
  homogeneous boundary conditions, we improve the interpolation using
  this information. */
  
  if (s.embed_gradient2 && s.boundary[0] != s.boundary_homogeneous[0]) {
    coord o = {x,y,z}, g;
    s.embed_gradient2 (point, s, &g);
    foreach_dimension()
      s[] += (o.x - p.x/nv)*g.x;
  }
}