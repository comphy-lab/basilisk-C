/**
# Axisymmetric coordinates

For problems with a symmetry of revolution around the $z$-axis of a
[cylindrical coordinate
system](http://en.wikipedia.org/wiki/Cylindrical_coordinate_system). The
longitudinal coordinate ($z$-axis) is *x* and the radial coordinate
($\rho$- or $r$-axis) is *y*. Note that *y* (and so *Y0*) cannot be
negative.

We first define a macro which will be used in some geometry-specific
code (e.g. [curvature computation](curvature.h)). */

#define AXI 1

/**
On trees we need refinement functions. */

#if TREE
static void refine_cm_axi (Point point, scalar cm)
{
#if !EMBED
  fine(cm,0,0) = fine(cm,1,0) = y - Delta/4.;
  fine(cm,0,1) = fine(cm,1,1) = y + Delta/4.;
#else // EMBED
  if (cs[] > 0. && cs[] < 1.) {
    coord n = interface_normal (point, cs);
    // Better? involve fs (troubles w prolongation)
    // coord n = facet_normal (point, cs, fs);
    foreach_child() {
      if (cs[] > 0. && cs[] < 1.) {
	coord p;
    	double alpha = plane_alpha (cs[], n);
	plane_center (n, alpha, cs[], &p);
	cm[] = (y + Delta*p.y)*cs[];
      }
      else
	cm[] = y*cs[];
    }
  }
  else
    foreach_child()
      cm[] = y*cs[];
#endif // EMBED
}

static void refine_face_x_axi (Point point, scalar fm)
{
#if !EMBED
  if (!is_refined(neighbor(-1))) {
    fine(fm,0,0) = y - Delta/4.;
    fine(fm,0,1) = y + Delta/4.;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors) {
    fine(fm,2,0) = y - Delta/4.;
    fine(fm,2,1) = y + Delta/4.;
  }
  fine(fm,1,0) = y - Delta/4.;
  fine(fm,1,1) = y + Delta/4.;
#else // EMBED
  double sig = 0., ff = 0.;
  if (cs[] > 0. && cs[] < 1.) {
    coord n = facet_normal (point, cs, fs);
    sig = sign(n.y)*Delta/4.;
  }
  if (!is_refined(neighbor(-1))) {
    ff = fine(fs.x,0,0);
    fine(fm,0,0) = (y - Delta/4. - sig*(1. - ff))*ff;
    ff = fine(fs.x,0,1);
    fine(fm,0,1) = (y + Delta/4. - sig*(1. - ff))*ff;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors) {
    ff = fine(fs.x,2,0);
    fine(fm,2,0) = (y - Delta/4. - sig*(1. - ff))*ff;
    ff = fine(fs.x,2,1);
    fine(fm,2,1) = (y + Delta/4. - sig*(1. - ff))*ff;
  }
  ff = fine(fs.x,1,0);
  fine(fm,1,0) = (y - Delta/4. - sig*(1. - ff))*ff;
  ff = fine(fs.x,1,1);
  fine(fm,1,1) = (y + Delta/4. - sig*(1. - ff))*ff;
#endif // EMBED
}

static void refine_face_y_axi (Point point, scalar fm)
{
#if !EMBED
  if (!is_refined(neighbor(0,-1)))
    fine(fm,0,0) = fine(fm,1,0) = max(y - Delta/2., 1e-20);
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors)
    fine(fm,0,2) = fine(fm,1,2) = y + Delta/2.;
  fine(fm,0,1) = fine(fm,1,1) = y;
#else // EMBED
  if (!is_refined(neighbor(0,-1))) {
    fine(fm,0,0) = (max(y - Delta/2., 1e-20))*fine(fs.y,0,0) ;
    fine(fm,1,0) = (max(y - Delta/2., 1e-20))*fine(fs.y,1,0);
  }
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors) {
    fine(fm,0,2) = (y + Delta/2.)*fine(fs.y,0,2);
    fine(fm,1,2) = (y + Delta/2.)*fine(fs.y,1,2);
  }
  fine(fm,0,1) = y*fine(fs.y,0,1);
  fine(fm,1,1) = y*fine(fs.y,1,1);
#endif // EMBED
}
#endif

/**
If embedded solids are presents, *cm*, *fm* and the fluxes need to be
updated consistently with the axisymmetric cylindrical coordinates and
the solid fractions. */

#if EMBED
double axi_factor (Point point, coord p) {
  return (y + p.y*Delta);
}

void cm_update (scalar cm, scalar cs, face vector fs)
{
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {
      coord p, n = facet_normal (point, cs, fs);
      double alpha = plane_alpha (cs[], n);
      plane_center (n, alpha, cs[], &p);
      cm[] = (y + Delta*p.y)*cs[];
    }
    else
      cm[] = y*cs[];
  }
  cm[top] = dirichlet(y*cs[]);
  cm[bottom] = dirichlet(y*cs[]);
}

void fm_update (face vector fm, scalar cs, face vector fs)
{
  foreach_face(x) {
    double sig = 0.;
    if (cs[] > 0. && cs[] < 1.) {
      coord n = facet_normal (point, cs, fs);
      sig = sign(n.y)*Delta/2.;
    }
    fm.x[] = (y - sig*(1. - fs.x[]))*fs.x[];
  }
  foreach_face(y)
    fm.y[] = max(y, 1e-20)*fs.y[];
  fm.t[top] = dirichlet(y*fs.t[]);
  fm.t[bottom] = dirichlet(y*fs.t[]);
}
#endif // EMBED

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
  Metric factors must be taken into account for fluxes on embedded
  boundaries. */

#if EMBED
  metric_embed_factor = axi_factor;
#endif  
  
  /**
  The volume/area of a cell is proportional to $r$ (i.e. $y$). We need
  to set boundary conditions at the top and bottom so that *cm* is
  interpolated properly when refining/coarsening the mesh. */

  scalar cmv = cm;
  foreach()
    cmv[] = y;
  cm[top] = dirichlet(y);
  cm[bottom] = dirichlet(y);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the axis of revolution is zero ($y=r=0$ on the axis). To avoid
  division by zero we set it to epsilon (note that mathematically the
  limit is well posed). */

  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face()
    fmv.x[] = max(y, 1./HUGE);
  fm.t[top] = dirichlet(y);
  fm.t[bottom] = dirichlet(y);
  
  /**
  We set our refinement/prolongation functions on trees. */

#if TREE
  cm.refine = cm.prolongation = refine_cm_axi;
  fm.x.prolongation = refine_face_x_axi;
  fm.y.prolongation = refine_face_y_axi;
#endif
}

/**
## See also

* [Axisymmetric streamfunction](axistream.h)
*/
