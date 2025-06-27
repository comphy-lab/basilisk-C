#define quadratic_lagrange(x,x1,a1,x2,a2,x3,a3) \
  (a1*(x - x2)*(x - x3)/(x1 - x2)/(x1 - x3) +  \
   a2*(x - x1)*(x - x3)/(x2 - x1)/(x2 - x3) +  \
   a3*(x - x1)*(x - x2)/(x3 - x1)/(x3 - x2))


foreach_dimension()
static inline double embed_extrapolate2_x (Point point, scalar s,
            coord n, coord p, double sb)
{
  // assert (cs[] > 0. && cs[] < 1.);
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

double embed_extrapolate2 (Point point, scalar s, coord n, coord p, double sb)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return embed_extrapolate2_x (point, s, n, p, sb);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return embed_extrapolate2_x (point, s, n, p, sb);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return embed_extrapolate2_y (point, s, n, p, sb);
  return embed_extrapolate2_z (point, s, n, p, sb);
#endif // dimension == 3
  return nodata;
}
