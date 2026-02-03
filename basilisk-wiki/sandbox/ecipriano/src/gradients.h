/**
## Dirichlet boundary condition

This function returns the gradient of scalar *s*, normal to the
embedded boundary defined by *cs*, of unit normal vector *n*
(normalised using the Euclidean norm, not the box norm) and of
centroid *p*. The Dirichlet boundary condition *bc* is imposed on the
embedded boundary.

The calculation follows [Johansen and Colella, 1998](#johansen1998)
and is summarised in the figure below (see also Figure 4 of Johansen
and Colella and Figure 2 of [Schwartz et al, 2006](#schwartz2006) for
the 3D implementation).

![Third-order normal gradient scheme](figures/dirichlet_gradient.svg) 

For degenerate cases, a non-zero value of *coef* is returned and
`coef*s[]` must be added to the value returned to obtain the gradient. */

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s,
             scalar cs, face vector fs,
             coord n, coord p, double bc,
             double * coef, bool vof = false)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
          cs[i,j-1] && cs[i,j] && cs[i,j+1])
        v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
        if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
            !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
            !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
          defined = false;
      if (defined)
        // bi-quadratic interpolation
        v[l] =
          quadratic (z,
               quadratic (y1,
              (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
               quadratic (y1,
              (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
               quadratic (y1,
              (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
        break;
    }
  if (v[0] == nodata) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
  
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata && vof) // vof-biased scheme
    return (cs[]*(bc - v[0])/d[0] + (1. - cs[])*(bc - v[1])/d[1])/Delta;
  else if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  else
    return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient (Point point, scalar s, scalar cs, face vector fs,
         coord n, coord p, double bc, double * coef, bool vof = false)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, fs, n, p, bc, coef, vof);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, fs, n, p, bc, coef, vof);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, fs, n, p, bc, coef, vof);
  return dirichlet_gradient_z (point, s, cs, fs, n, p, bc, coef, vof);
#endif // dimension == 3
  return nodata;
}

/**
The interface gradients for VOF-tracers can be computed using the following
function, provided the volume and surface fractions `cs` and `fs` which define
the region where the tracers `s` is defined, and the value of the field at the
interface (`bc`). */

double plic_gradient (Point point, scalar s, scalar cs, face vector fs,
    double bc, bool vof = true)
{
  coord m = interface_normal (point, cs), p;
  double alpha = plane_alpha (cs[], m);
  plane_area_center (m, alpha, &p);
  normalize (&m);
  double coef = 0.;
  double grad = dirichlet_gradient (point, s, cs, fs, m, p, bc, &coef, vof);
  return (grad == nodata || coef != 0.) ? 0. : grad;
}

