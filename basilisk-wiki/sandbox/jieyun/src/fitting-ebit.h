/** 
# Different fitting functions used in the EBIT method
*/

/**
## General purpose geometric functions
 
 Test if two points x3, x4 are located at the same side of the
 straight line (x1, x2). For the following special cases, it will return false:
 1. x1 == x2
 2. (x1, x2, x3) or (x1, x2, x4) are collinear
*/
static bool is_same_side (coord x1, coord x2, coord x3, coord x4) {
  coord dx1, dx3, dx4;
  foreach_dimension() {
    dx1.x = x1.x - x2.x;
    dx3.x = x3.x - x2.x;
    dx4.x = x4.x - x2.x;
  }
  double ds = (dx3.x*dx1.y - dx3.y*dx1.x)*(dx4.x*dx1.y - dx4.y*dx1.x);
  return ds > 0.;
}

/** 
 Test if two line segments (x1, x2) and (x3, x4) intersect.
 The following configurations are considered as intersecting:
 1. x1 or x2 is exaclty located on (x3, x4), including the end points, and vice versa
 if x1 == x2 && x3 == x4, it will return true, fix this.
 if (x1, x2) and (x3, x4) are collinear, it always returns true, not correct.
*/
bool is_intersect (coord x1, coord x2, coord x3, coord x4) {
  if (is_same_side (x1, x2, x3, x4) || is_same_side (x3, x4, x1, x2))
    return false;
  else
    return true;
}


/** Function for calculating the intersection of two line segments
  seg. 1: (x0, ymin) - (x0, ymax); seg 2: (x1, y1) - (x2, y2)
  return -1 if there is no intersection. */
static double my_intersect (double x1, double y1, double x2, double y2, double x0, \
  double ymin = 0., double ymax = 1.) {
  double issame, y0;
  issame = (x2 - x0)*(x1 - x0);

  if (issame > 0.) {
    y0 = -HUGE;
  }
  else {
    double m = (y2 - y1)/(x2 - x1);
    y0 = m*(x0 - x1) + y1;
    if (y0 > ymax || y0 < ymin)
      y0 = -HUGE;
  }

  return y0;
}


/** Test if a point xp is inside the polygon formed by vertices xv. */
bool is_inside (coord *xv, int nv, coord xp, double epsy = 1.e-32) {
  double xmin = HUGE, ymin = HUGE, xmax = -HUGE, ymax = -HUGE;
  for (int iv = 0; iv < nv; iv++) {
    if (xv[iv].x < xmin)
      xmin = xv[iv].x;
    if (xv[iv].y < ymin)
      ymin = xv[iv].y;

    if (xv[iv].x > xmax)
      xmax = xv[iv].x;
    if (xv[iv].y > ymax)
      ymax = xv[iv].y;
  }

  if (xp.x > xmax || xp.x < xmin || xp.y > ymax || xp.y < ymin)
    return false;
  // chose a point ouside the polygon, test the corner cases later.
  double x2 = max(xp.x + 1., 2*xmax - xmin);
  // epsy: In order to avoid some corner cases, such as (xp, xo) pass through the vertex of polygon
  coord xo = {x2, xp.y + epsy};
  int nint = 0;
  for (int iv = 0; iv < nv; iv++) {
    int ivp = (iv + 1) % nv;
    if (is_intersect (xv[iv], xv[ivp], xp, xo))
      if (xv[iv].y < xo.y || xv[ivp].y < xo.y) // for corner case
        nint++;
  }

  if (nint % 2 == 0)
    return false;
  else
    return true;
}


/** 
## Parabolic fit */

#include "parabola.h"
int fit_parabolic (coord *xp, int np, double *r) {
  *r = -1.;

  coord m, fc = {0., 0., 0.};
  if (np == 4) {
    m.x = xp[2].y - xp[3].y;
    m.y = xp[3].x - xp[2].x;
  }
  else {
    m.x = xp[2].y - xp[0].y;
    m.y = xp[0].x - xp[2].x;
  }
  double ds = sqrt(sq(m.x) + sq(m.y)) + 1.e-32;

  m.x /= ds;
  m.y /= ds;

  for (int i = 0; i < np; i++) {
    foreach_dimension()
      fc.x += xp[i].x/np;
  }

  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  for (int i = 0; i < np; i++) {
    coord fcc = {xp[i].x, xp[i].y, xp[i].z};
    parabola_fit_add (&fit, fcc, 1.);
  }

  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL);

  *r = 1./kappa;

  return 0;
}


/**
## Functions for polynomial fit */

/** Evaluate the value of a no-order polynomial */
double evaluate_poly (double *a, int no, double x) {
  double y = 0.;
  for (int i = no - 1; i >= 0; i--)
    y = y*x + a[i];

  return y;
}


/** Roots findind subroutine for a cubic polynomial,
 a[i] is the coefficient of the x^i term, [xmin, xmax] is the
 range where we want to search the roots x_i. Return the the number of roots. */
int roots_cubic (double *ax, double xmin, double xmax, double *x) {
  int nroots = 0;
  double maxa = 1.e-12, epsa = 1.e-6;
  for (int i = 1; i < 4; i++)
    maxa = max(maxa, fabs(ax[i]));

  double a[4], xr[3];
  int nr = 0;
  for (int i = 0; i < 4; i++)
    a[i] = ax[i];

  if (fabs(a[3]/maxa) <= epsa) {
    if(fabs(a[2]/maxa) <= epsa) {
      // linear equation
      xr[0] = -a[0]/(a[1] + 1.e-32);
      nr = 1;
    }
    else {
      // quadratic polynomial
      double delta = sq(a[1]) - 4.*a[2]*a[0];
      if (delta >= 0.) {
        delta = sqrt(delta);
        xr[0] = (-a[1] + delta)/(2.*a[2]);
        xr[1] = (-a[1] - delta)/(2.*a[2]);
        nr = 2;
      }
      else
        return 0;
    }
  }
  else {
    double a3 = ax[3];
    for (int i = 0; i < 4; i++)
      a[i] = ax[i]/a3;

    double q = a[1]/3. - sq(a[2])/9.;
    double r = (a[1]*a[2] - 3.*a[0])/6. - cube(a[2])/27.;
    double delta = sq(r) + cube(q);

    if (delta <= 0.) {
      // three real roots
      double q3 = sqrt(cube(-q));
      double theta = (q == 0.) ? 0. : acos(r/q3);
      double p23 = 2.*pi/3.;

      theta /= 3.;
      double phi[3] = {theta, theta - p23, theta + p23};
      for (int i = 0; i < 3; i++) {
        xr[i] = 2.*sqrt(-q)*cos(phi[i]) - a[2]/3.;
        nr = 3;
      }
    }
    else {
      // only one real root
      double aa = pow(fabs(r) + sqrt(delta), 1./3.);
      double t1 = aa - q/aa;
      if (r < 0)
        t1 *= -1.;
      xr[0] = t1 - a[2]/3.;
      nr = 1;
    }
  }

  for (int i = 0; i < nr; i++) {
    double xx = xr[i];
    if (xx <= xmax && xx >= xmin) {
      x[nroots] = xx;
      nroots++;
    }
  }
  return nroots;
}


/** 
## Local spline fit proposed by [Dritschel, 1988](#dritchel1988). */

/** Local spline fitting method */
int fit_local_spline (coord x1, coord x2, double kp1, double kp2, double *ax, double *ay) {
  coord dt, dn;
  double ei = 0;

  foreach_dimension() {
    dt.x = x2.x - x1.x;
    ei += sq(dt.x);
  }
  dn.x = -dt.y;
  dn.y = dt.x;
  ei = sqrt(ei);

  double alphai = -ei*(kp1/3. + kp2/6.);
  double betai = 0.5*ei*kp1;
  double gammai = ei*(kp2 - kp1);

  ax[0] = x1.x;
  ax[1] = dt.x + alphai*dn.x;
  ax[2] = betai*dn.x;
  ax[3] = gammai*dn.x;

  ay[0] = x1.y;
  ay[1] = dt.y + alphai*dn.y;
  ay[2] = betai*dn.y;
  ay[3] = gammai*dn.y;

  return 1;
}


/** Return the curvature of a local spline fitting. */
inline double kappa_local_spline (double kp1, double kp2, double xi) {
  return kp1 + (kp2 - kp1)*xi;
}


/** Return the unit tangential of a local spline. */
coord tn_local_spline (coord x1, coord x2, double kp1, double kp2, double xi) {
  coord dt, dn;
  double ei = 0;

  foreach_dimension() {
    dt.x = x2.x - x1.x;
    ei += sq(dt.x);
  }
  dn.x = -dt.y;
  dn.y = dt.x;
  ei = sqrt(ei);

  double alphai = -ei*(kp1/3. + kp2/6.);
  double betai = 0.5*ei*kp1;
  double gammai = ei*(kp2 - kp1);
  double dedxi = alphai + 2.*betai*xi + 3.*gammai*sq(xi);
  coord tn = {dt.x + dedxi*dn.x, dt.y + dedxi*dn.y};

  // normalize the tangential
  ei *= sqrt(1. + sq(dedxi));
  tn.x /= ei;
  tn.y /= ei;
  return tn;
}


/** Return the intersection point between a local spline and a
 horizontal straight line, x = f(eta_0) = x0, return y = g(eta_0). */
double intersection_cubic (double *ax, double *ay, double x0) {
  double _ax[4];
  for (int i = 0; i < 4; i++)
    _ax[i] = ax[i];

  _ax[0] -= x0;

  int nr;
  double phi[3];
  nr = roots_cubic (_ax, 0., 1., phi);
  if (nr == 1)
    return evaluate_poly (ay, 4, phi[0]);
  else
    return -1.;
}


double get_para_cubic (double *ax, double x0) {
  double _ax[4];
  for (int i = 0; i < 4; i++)
    _ax[i] = ax[i];

  _ax[0] -= x0;

  int nr;
  double phi[3];
  nr = roots_cubic (_ax, 0., 1., phi);
  if (nr == 1)
    return phi[0];
  else
    return -1.;
}


/** 
## Circle fit */

/** Get the center and radius of circle defined by three points, 
  return negative radius for colinear points*/
void get_circle (double x1, double y1, double x2, double y2, double x3, double y3, \
  double *xc, double *yc, double *r) {
  double a1, a2, xxc, yyc;

  a1 = (y3 - y1)*(x3 - x2) - (y3 - y2)*(x3 - x1);
  a2 = (x3 - x1)*(x3 - x2) + (y3 - y1)*(y3 - y2);

  if (fabs(a1) < 1.e-12) {
    // colinear three points or any two points are too close
    *xc = 0.;
    *yc = 0.;
    *r = -1.;
  }
  else {
    a2 /= a1;
    xxc = 0.5*((x1 + x2) + a2*(y2 - y1));
    yyc = 0.5*((y1 + y2) + a2*(x1 - x2));
    *xc = xxc;
    *yc = yyc;
    *r = sqrt(sq(x3 - xxc) + sq(y3 - yyc));
  }
}

/** Calculate the curvature by fitting three points by a circle.
 Sign of the curvature is based on n = (-t_y, t_x), t = p2 - p1.
*/
double get_kappa_circle (coord p1, coord p2, coord p3) {
  double kappa, eps_ds = 1.e-12;
  coord t13 = {p3.x - p1.x, p3.y - p1.y};
  coord t23 = {p3.x - p2.x, p3.y - p2.y};
  coord t12 = {p2.x - p1.x, p2.y - p1.y};
  coord n12 = {p1.y - p2.y, p2.x - p1.x};

  double a1 = t13.y*t23.x - t23.y*t13.x;
  double a2 = t13.x*t23.x + t13.y*t23.y;

  double ds23 = sqrt(sq(t23.x) + sq(t23.y));
  double ds13 = sqrt(sq(t13.x) + sq(t13.y));
  double ds12 = sqrt(sq(t12.x) + sq(t12.y));
  double ds = ds23*ds13;

  if (ds12 < eps_ds || ds13 < eps_ds || ds23 < eps_ds) {
    // Any two points are too close
    return 0.;
  }
  else if (fabs(a1)/(ds + 1.e-12) < 1.e-12) {
    // Three points are colinear
    return 0.;
  }
  else {
    a2 /= a1;
    coord pc = {0.5*((p1.x + p2.x) + a2*t12.y), 0.5*((p1.y + p2.y) - a2*t12.x)};
    kappa = 1./sqrt(sq(p3.x - pc.x) + sq(p3.y - pc.y));
    double signk = (pc.x - 0.5*(p1.x + p2.x))*n12.x \
      + (pc.y - 0.5*(p1.y + p2.y))*n12.y;
    kappa *= sign(signk);
  }

  return kappa;
}


/** Fit more than two points by a circle, which passes through the first two points.
 Not a robust method.
*/
int fit_circle (coord *xp, int np, double *xc, double *yc, double *r) {
  double x1, x2, y1, y2, xm, ym, nx, ny, dd;

  *xc = 0.;
  *yc = 0.;
  *r = -1.;

  x1 = xp[0].x; x2 = xp[1].x;
  y1 = xp[0].y; y2 = xp[1].y;

  xm = 0.5*(x1 + x2);
  ym = 0.5*(y1 + y2);

  nx = y1 - y2;
  ny = x2 - x1;

  dd = sqrt(nx*nx + ny*ny);
  if (dd < 1.e-16)
    return -1;

  nx /= dd;
  ny /= dd;

  dd = 0.5*dd;
  double a1=0., a2=0., a3=0.;

  for (int ip = 0; ip < np; ip++) {
    double dx, dy, tmp1, tmp2;
    dx = xp[ip].x - xm;
    dy = xp[ip].y - ym;
    tmp1 = dx*dx + dy*dy - dd*dd;
    tmp2 = dx*nx + dy*ny;
    a1 += tmp1*tmp2;
    a2 += tmp2*tmp2;
    a3 += (tmp2*tmp2/(dx*dx + dy*dy + 1.e-32));
  }
  // all points are colinear
  if (a3 < 1.e-32)
    return -1;
  
  a3 = 0.5*a1/a2;
  *xc = xm + a3*nx;
  *yc = ym + a3*ny;
  *r = sqrt(a3*a3 + dd*dd);
  return 0;
}


double intersection_circle (double xc, double yc, double r, double x0) {
  double delta, y0, y1, y2;

  delta = sq(r) - sq(x0 - xc);

  if (delta < 0.) {
    return -1.;
  }
  else {
    delta = sqrt(delta);
    y1 = yc + delta;
    y2 = yc - delta;
    
    y0 = -1.;
    if (y1 > 0. && y1 < 1.) {
      y0 = y1;
    }
    else if(y2 > 0. && y2 < 1.) {
      y0 = y2;
    }
    return y0;
  }
}

int intersection_circle_two (double xc, double yc, double r, double x0, double *y) {
  double delta;

  y[0] = nodata;
  y[1] = nodata;
  delta = sq(r) - sq(x0 - xc);

  if (delta < 0.)
    return 0;
  else {
    delta = sqrt(delta);
    y[0] = yc + delta;
    y[1] = yc - delta;

    return delta == 0. ? 1 : 2;
  }
}


/**
## Circular interface initialization

*/

/** Initialize a circular interface using bisection method, we need the intersection points on cell edge. 
 Which function in vofi should I use to replace this part? */ 
#define circle_sign(x, y, xc, yc, ra) (sq(x - xc) + sq(y - yc) - sq(ra))

void init_circle (double xc, double yc, double ra, \
  scalar f, face vector fs, int itmax = 40) {
  coord dir = {0., 1.};
  foreach_face() {
    double x1, y1, x2, y2, x3, y3, eps, s1, s2;
    x1 = x - 0.5*Delta*dir.x;
    y1 = y - 0.5*Delta*dir.y;

    x2 = x + 0.5*Delta*dir.x;
    y2 = y + 0.5*Delta*dir.y;

    s1 = circle_sign(x1, y1, xc, yc, ra);
    s2 = circle_sign(x2, y2, xc, yc, ra);
    if (s1*s2 < 0.) {
      eps = 1.e-16*Delta;
      for (int it = 0; it <= itmax; it++) {
        x3 = 0.5*(x1 + x2);
        y3 = 0.5*(y1 + y2);
        if (x3 - x1 < eps && y3 - y1 < eps)
          break;
        s1 = circle_sign(x1, y1, xc, yc, ra);
        s2 = circle_sign(x3, y3, xc, yc, ra);
        if (s1*s2 < 0.) {
          x2 = x3; y2 = y3;
        }
        else {
          x1 = x3; y1 = y3;
        }
      }
      fs.x[] = ((x3 - x)*dir.x+ (y3 - y)*dir.y)/Delta + 0.5;
    }
  }

  boundary ((scalar *) {fs});
}

/** Initialize a circular interface using analytical formulae. */ 
#if dimension <= 2
void init_circle_analytical (double xc, double yc, double ra, \
  scalar f, face vector fs, double rb = -1.) {

  rb = (rb < 0.) ? ra : rb;
  coord cc = {xc, yc}, rc = {ra, rb};

  foreach_face() {
    coord xo = {x, y};
    double del = sq(rc.y) - sq(rc.y/rc.x)*sq(xo.x - cc.x);
    double y3 = xo.y - 0.5*Delta;
    if (del >= 0.) {
      double x1 = sqrt(fabs(del));
      double y1 = (cc.y + x1 - y3)/Delta;
      double y2 = (cc.y - x1 - y3)/Delta;
      fs.x[] = (y1 >= 0. && y1 <= 1.) ? y1 : ((y2 >= 0. && y2 <= 1.) ? y2 : 0.);
    }
  }
}
#endif

/**
## References

~~~bib
@article{dritchel1988,
  title={Contour surgery: a topological reconnection scheme 
    for extended integrations using contour dynamics},
  author={Dritschel, David G.},
  journal={Journal of Computational Physics},
  volume={77},
  pages={240-266},
  year={1988},
  publisher={Elsevier}
}
~~~
*/
