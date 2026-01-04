/**
# Distance from a point to an ellipse

See section 2.9 of
[DistancePointEllipseEllipsoid.pdf](https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf). */

static double RobustLength (double v0, double v1)
{
  return fabs(v0) > fabs(v1) ?
    fabs(v0)*sqrt(1. + sq(v1/v0)) :
    fabs(v1)*sqrt(1. + sq(v0/v1));
}

static double GetRoot (double r0, double z0, double z1, double g)
{
  double n0 = r0*z0 ;
  double s0 = z1 - 1., s1 = g < 0. ? 0. : RobustLength (n0, z1) - 1.;
  double s = 0.;
  for (int i = 0 ; i < 100; ++i) {
    s = (s0 + s1)/2.;
    if (s == s0 || s == s1 ) break;
    double ratio0 = n0/(s + r0), ratio1 = z1/(s + 1.) ;
    g = sq (ratio0) + sq (ratio1) - 1.;
    if (g > 0.) s0 = s;
    else  if (g < 0.) s1 = s;
    else break;
  }
  return s ;
}

double DistancePointEllipse (double e0, double e1, double y0, double y1,
			     double * x0, double * x1)
{
  bool sym0 = false, sym1 = false;
  if (y0 < 0.)
    y0 = - y0, sym0 = true;
  if (y1 < 0.)
    y1 = - y1, sym1 = true;
  
  double distance ;
  if (y1 > 0.) {
    if (y0 > 0.) {
      double z0 = y0/e0, z1 = y1/e1, g = sq(z0) + sq(z1) - 1.;
      if (g != 0.) {
	double r0 = sq(e0/e1), sbar = GetRoot (r0, z0, z1, g);
	*x0 = r0*y0/(sbar + r0);
	*x1 = y1/(sbar + 1.);
	distance = sqrt (sq(*x0 - y0) + sq(*x1 - y1));
      }
      else {
	*x0 = y0;
	*x1 = y1;
	distance = 0.;
      }
    }
    else {
      // y0 == 0
      *x0 = 0.;
      *x1 = e1;
      distance = fabs(y1 - e1);
    }
  }
  else {
    // y1 == 0
    double numer0 = e0*y0, denom0 = sq(e0) - sq(e1);
    if (numer0 < denom0) {
      double xde0 = numer0/denom0;
      *x0 = e0*xde0;
      *x1 = e1*sqrt(1. - xde0*xde0);
      distance = sqrt(sq(*x0 - y0) + sq(*x1));
    }
    else {
      *x0 = e0;
      *x1 = 0.;
      distance = fabs(y0 - e0);
    }
  }

  if (sym0) *x0 = - *x0;
  if (sym1) *x1 = - *x1;
  
  return sign(sq(y0/e0) + sq(y1/e1) - 1.)*distance;
}
