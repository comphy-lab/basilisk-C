/**
This is a copy from [dirichlet-boundary-condition](/src/embed.h#dirichlet-boundary-condition). This function is used to compute the gradient, normal to the gas-liquid interface, of a tracer concentration field confined to one side of the interface. The scalar *cs* is the vof field of the phase where the tracer is confined and *bc* is the (known) value of concentration at the interface centroid *p*. The boolean variable *third* is used to decide whether the gradient is estimated with the third-order scheme or the vof-average one.
*/

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double concentration_gradient_x (Point point, scalar s, scalar cs, face vector fs,
					   coord n, coord p, double bc,
					   bool third)
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
    This is a degenerate case, we set the gradient to zero. */
	
    return 0.;
  }

  /**
  For non-degenerate cases, the gradient can be  obtained using
  second-order, third-order or vof-averaged estimates. */
 
  if (v[1] != nodata && third) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  else if (v[1] != nodata)
    return (cs[]*(bc - v[0])/d[0] + (1. - cs[])*(bc - v[1])/d[1])/Delta; //vof-avg
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double concentration_gradient (Point point, scalar s, scalar cs, face vector fs,
				coord n, coord p, double bc, bool third)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return concentration_gradient_y (point, s, cs, fs, n, p, bc, third);
  return concentration_gradient_z (point, s, cs, fs, n, p, bc, third);
#endif // dimension == 3
  return nodata;
}