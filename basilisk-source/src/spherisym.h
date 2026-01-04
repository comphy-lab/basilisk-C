/**
# Spherically-symmetric coordinates

This file defines the metric coefficients for a (one-dimensional)
[spherically-symmetric](https://en.wikipedia.org/wiki/Circular_symmetry#Spherical_symmetry)
coordinate system.

The radial coordinate $r$ is *x*.  Note that *x* (and so *X0*) cannot
be negative.

We first define a macro which will be used in some geometry-specific
code (e.g. [viscous stress tensor](viscosity.h)). */

#define SPHERISYM 1

/**
On trees we need refinement functions. */

#if TREE
static void refine_cm_spherisym (Point point, scalar cm)
{
  fine(cm,0) = sq (x - Delta/4.);
  fine(cm,1) = sq (x + Delta/4.);
}

static void refine_face_x_spherisym (Point point, scalar fm)
{
  if (!is_refined(neighbor(-1)))
    fine(fm,0) = sq (x - Delta/2.);
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors)
    fine(fm,2) = sq (x + Delta/2.);
  fine(fm,1) = sq(x);
}
#endif // TREE

event metric (i = 0) {

  /**
  By default *cm* is a constant scalar field. To make it variable, we
  need to allocate a new field. We also move it at the begining of the
  list of variables: this is important to ensure the metric is defined
  before other fields. */

  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }

  /**
  The volume/area of a cell is proportional to $r^2$ (i.e. $x^2$). We need
  to set boundary conditions at the top and bottom so that *cm* is
  interpolated properly when refining/coarsening the mesh. */

  scalar cmv = cm;
  foreach()
    cmv[] = x*x;
  cm[left] = dirichlet(x*x);
  cm[right] = dirichlet(x*x);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the center of spherical symmetry is zero ($x=r=0$ in the
  center). To avoid division by zero we set it to epsilon (note that
  mathematically the limit is well posed). */

  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face()
    fmv.x[] = max(x*x, 1e-20);
  
  /**
  We set our refinement/prolongation functions on trees. */

#if TREE
  cm.refine = cm.prolongation = refine_cm_spherisym;
  fm.x.prolongation = refine_face_x_spherisym;
#endif
}
