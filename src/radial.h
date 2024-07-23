/**
# Radial/cylindrical coordinates

This implements the radial coordinate mapping illustrated below.

![Radial coordinate mapping](figures/radial.svg)

It works in 1D, 2D and 3D. The 3D version corresponds to cylindrical
coordinates since the $z$-coordinate is unchanged. 

The only parameter is $d\theta$, the total angle of the sector. */

double dtheta = pi/3.;

/**
For convenience we add definitions for the radial and angular
coordinates $(r, \theta)$. */

map()
{
  double r = x, theta = y*dtheta/L0;
  NOT_UNUSED(r); NOT_UNUSED(theta);
}

event metric (i = 0) {

  /**
  We initialise the scale factors, taking care to first allocate the
  fields if they are still constant. */

  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }
  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }

  /**
  The area (in 2D) of a mapped element is the area of an annulus
  defined by the two radii $r-\Delta/2$ and $r+\Delta/2$, divided by
  the total number of sectors $N=2\pi L0/(d\theta\Delta)$, this gives
  $$
  \frac{\pi [(r + \Delta / 2)^2 - (r - \Delta / 2)^2]}{2 \pi L 0 / (d \theta
  \Delta)} = \frac{2 \pi r \Delta}{2 \pi L 0 / (d \theta \Delta)} = \frac{rd
  \theta}{L 0} \Delta^2
  $$
  By definition, the (area) metric factor `cm` is the mapped area divided
  by the unmapped area $\Delta^2$. */

  scalar cmv = cm;
  foreach()
    cmv[] = r*dtheta/L0;

  /**
  It is important to set proper boundary conditions, in particular
  when refining the grid. */
  
  cm[left] = dirichlet (r*dtheta/L0);
  cm[right] = dirichlet (r*dtheta/L0);
  
  /**
  The (length) metric factor `fm` is the ratio of the mapped length of
  a face to its unmapped length $\Delta$. In the present case, it is
  unity for all dimensions except for the $x$ coordinates for which it
  is the ratio of the arclength by the unmapped length $\Delta$. We
  also set a small minimal value to avoid division by zero, in the
  case of a vanishing inner radius. */
  
  face vector fmv = fm;
  foreach_face()
    fmv.x[] = 1.;
  foreach_face(x)
    fmv.x[] = max(r*dtheta/L0, 1e-20);
}
