/*******************************************************************************
  embed_evaluate function
******************************************************************************/

#define quadratic(x,a1,a2,a3)						\
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))


bool emerged = true; // emerged = true -> updated cell
scalar csm1[]; // cs at previous iteration

foreach_dimension()
static inline void embed_evaluate_x (Point point, scalar s, scalar cs,
				     coord n, coord b,
				     double * d0, double * v0,
				     double * d1, double * v1)
{
  assert ((cs[]) > 0. && (cs[]) < 1.);
  
  /**
  We use the normal pointing from solid to fluid. */
  
  foreach_dimension()
    n.x = -n.x;

  /**
  We compute the distance from the cell center *d* and the
  interpolated values *v*. */
  
  double d[2], v[2] = {nodata,nodata};
  
  /**
  We then assess if neighboring cells in the fluid domain are
  accessible (fs > 0). */
  
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;

  /**
  If neighbors are available, we perform a quadric (2D) or biquadratic
  (3D) interpolations. If not all neighbors are available, we do not
  compute a value to avoid using values of cells of the grid that are
  not topologically connected, which could prevent convergence of the
  multigrid solver. */
  
  if (defined)
    for (int l = 0; l <= 1; l++) {
      //distance from p (point on the solid boundary)
      int i = (l + 1)*sign(n.x);
      d[l] = (i - b.x)/(n.x);
      // projection in the y-direction
      double y1 = (b.y) + (d[l])*(n.y);
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
      d[l] = (i - b.x)/(n.x)-fabs(b.x/n.x); // distance from the cell center
#if dimension == 2
      // quadratic interpolation
    if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] && 
    cs[i,j-1] && cs[i,j] && cs[i,j+1] && 
    (emerged || (csm1[i,j-1] && csm1[i,j] && csm1[i,j+1])))  
	    v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      // projection in the z-direction
      double z1 = (b.z) + (d[l])*(n.z);
      int k = z1 > 0.5 ? 1 : z1 < -0.5 ? -1 : 0;
      z1 -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1] ||
	    (!emerged && (!csm1[i,j+m,k-1] || !csm1[i,j+m,k] || !csm1[i,j+m,k+1])))
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z1,
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
  *d0 = d[0]; *v0 = v[0];
  *d1 = d[1]; *v1 = v[1];
}

void embed_evaluate (Point point, scalar s, scalar cs,
		     coord n, coord b,
		     double * d0, double * v0,
		     double * d1, double * v1)
{
#if dimension == 2
  if (fabs(n.x) >= fabs(n.y))
    embed_evaluate_x (point, s, cs, n, b, d0, v0, d1, v1);
  else
    embed_evaluate_y (point, s, cs, n, b, d0, v0, d1, v1);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      embed_evaluate_x (point, s, cs, n, b, d0, v0, d1, v1);
    else
      embed_evaluate_z (point, s, cs, n, b, d0, v0, d1, v1);
  }
  else if (fabs(n.y) >= fabs(n.z))
    embed_evaluate_y (point, s, cs, n, b, d0, v0, d1, v1);
  else
    embed_evaluate_z (point, s, cs, n, b, d0, v0, d1, v1);
#endif // dimension == 3
}

/*****************************************************************************
  embed_prolong_wall2 function
*****************************************************************************/

#define linear_lagrange(x,x1,a1,x2,a2) \
( (a1-a2)*x/(x1-x2) + (x1*a2-x2*a1)/(x1-x2))

#define quadratic_lagrange(x,x1,a1,x2,a2,x3,a3) \
  (a1*(x - x2)*(x - x3)/(x1 - x2)/(x1 - x3) +  \
   a2*(x - x1)*(x - x3)/(x2 - x1)/(x2 - x3) +  \
   a3*(x - x1)*(x - x2)/(x3 - x1)/(x3 - x2))


foreach_dimension()
static inline double embed_prolong_wall2_x (Point point, scalar s,
            coord n, coord p, double sb)
{
  coord c = {0., 0., 0.}; // Extrapolate to cell center  
  double d[2], v[2] = {nodata,nodata};
  embed_evaluate (point, s, cs, n, c, &d[0], &v[0], &d[1], &v[1]);
  
  /**
  If only the boundary condition is defined, this is a degenerate
  case, we use the boundary value as the cell-centered value. */

  if (v[0] == nodata)
    return sb;

  /**
  For non-degenerate cases, the cell-centered value is obtained using
  either first-order (linear) or second-order (quadratic)
  estimates. */

  if (v[1] != nodata){
        return quadratic_lagrange (0.,-p.x/n.x,sb,d[0],v[0],d[1],v[1]);
      }
  else
    return (linear_lagrange (0.,-p.x/n.x,sb,d[0],v[0])); 
}

double embed_prolong_wall2 (Point point, scalar s, coord n, coord p, double sb)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return embed_prolong_wall2_x (point, s, n, p, sb);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return embed_prolong_wall2_x (point, s, n, p, sb);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return embed_prolong_wall2_y (point, s, n, p, sb);
  return embed_prolong_wall2_z (point, s, n, p, sb);
#endif // dimension == 3
  return nodata;
}
