/**
# Computation of distance from Height Functions
*/

#include "heights.h"

double minimum (double a, double b, double c, double x)
{
  double f;
  int nmax = 10;
  do {
    f = (x + (2.*a*x + b)*(x*(a*x + b) + c))/
      (a*(6.*x*(a*x + b) + 2.*c) + b*b + 1.);
    x -= f;
  }
  while (nmax-- && fabs (f) > 1e-5);
  return x;
}

#define XMIN 1.5

foreach_dimension()
coord height_closest_y (Point point, vector h)
{
  coord p = { nodata, nodata };
  if (h.y[-1] != nodata && h.y[] != nodata && h.y[1] != nodata) {
    double a = (h.y[1] + h.y[-1] - 2.*h.y[])/2.;
    double b = (h.y[1] - h.y[-1])/2.;
    double c = height (h.y[]);
    double xm = nodata, ym = nodata, dm = HUGE;
    for (double x = -1; x <= 1.; x += 1.) {
      double xp = minimum (a, b, c, x), yp = xp*(a*xp + b) + c;
      double dp = sq(xp) + sq(yp);
      if (dp < dm)
	xm = xp, ym = yp, dm = dp;
    }
    if (fabs (xm) < XMIN)
      p.x = xm*Delta, p.y = ym*Delta;
  }
  return p;
}

scalar da[], db[];

coord height_closest (Point point, vector h, int * s)
{
  coord p;
  coord qx = height_closest_x (point, h);
  coord qy = height_closest_y (point, h);
  if (qx.x != nodata) {
    *s = (orientation(h.x[]) ? height(h.x[]) < 0. : height(h.x[]) > 0.) ? 1 : -1;
    if (qy.x != nodata) {
#if 1
      double a = fabs (h.x[0,1] - h.x[0,-1])/2.;
      double b = fabs (h.y[1] - h.y[-1])/2.;
      a /= sqrt (1. + a*a);
      b /= sqrt (1. + b*b);
      a = pow (1. - a, 4);
      b = pow (1. - b, 4);
      da[] = a, db[] = b;
#if 0
      foreach_dimension()
	p.x = (a*qx.x + b*qy.x)/(a + b);
#else
      foreach_dimension()
	p.x = (qx.x + qy.x)/2.;
#endif
#else
      double a = fabs(qx.x)/(XMIN*Delta), b = fabs(qy.x)/(XMIN*Delta);
      foreach_dimension()
	p.x = ((1. - a)*qx.x + (1. - b)*qy.x)/(2. - a - b);
	// p.x = a < b ? qx.x : qy.x;
#endif
    }
    else
      p = qx;
  }
  else {
    if (qy.x != nodata)
      *s = (orientation(h.y[]) ? height(h.y[]) < 0. : height(h.y[]) > 0.) ? 1 : -1;
    p = qy;
  }
  return p;
}

void height_distance (vector h, scalar d, double weight = 1.)
{
  foreach() {
    int s;
    coord p = height_closest (point, h, &s);
    if (p.x != nodata)
      d[] = (1. - weight)*d[] + weight*s*sqrt(sq(p.x) + sq(p.y));
  }
}

#undef XMIN
