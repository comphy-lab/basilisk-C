/** Different fitting function for intersection calculation */

/** Function for calculating the intersection of two line segments
  seg. 1: (x0, 0) - (x0, 1); seg 2: (x1, y1) - (x2, y2)
  return -1 if there is no intersection. */
static double my_intersect(double x1, double y1, double x2, double y2, double x0){
  double issame, y0;
  issame = (x2 - x0) * (x1 - x0);

  if (issame > 0.){
    y0 = -1.;
  }
  else{
    double m = (y2 - y1) / (x2 - x1);
    y0 = m * (x0 - x1) + y1;
    if (y0 > 1. || y0 < 0.)
      y0 = -1.;
  }

  return y0;
}

/** Get the center and radius of circle defined by three points, 
  return negative radius for colinear points*/
void get_circle(double x1, double y1, double x2, double y2, double x3, double y3, \
  double *xc, double *yc, double *r){
  double alpha, a1, a2, xxc, yyc;

  a1 = (y3 - y1) * (x3 - x2) - (y3 - y2) * (x3 - x1);
  a2 = (x3 - x1) * (x3 - x2) + (y3 - y1) * (y3 - y2);

  if (fabs(a1) < 1.e-12){
    // colinear three points or any two points are too close
    *xc = 0.;
    *yc = 0.;
    *r = -1.;
  }
  else{
    a2 /= a1;
    xxc = 0.5 * ((x1 + x2) + a2 * (y2 - y1));
    yyc = 0.5 * ((y1 + y2) + a2 * (x1 - x2));
    *xc = xxc;
    *yc = yyc;
    *r = sqrt(sq(x3 - xxc) + sq(y3 - yyc));
  }
}


/** Fit more than two points by a circle, which passes through the first two points*/
int fit_circle(double *xp, double *yp, int np, double *xc, double *yc, double *r){
  double x1, x2, y1, y2, xm, ym, nx, ny, dd;

  *xc = 0.;
  *yc = 0.;
  *r = -1.;

  x1 = xp[0]; x2 = xp[1];
  y1 = yp[0]; y2 = yp[1];

  xm = 0.5 * (x1 + x2);
  ym = 0.5 * (y1 + y2);

  nx = y1 - y2;
  ny = x2 - x1;

  dd = sqrt(nx * nx + ny * ny);
  if (dd < 1.e-16)
    return -1;

  nx /= dd;
  ny /= dd;

  dd = 0.5 * dd;
  double a1=0., a2=0., a3=0.;

  for(int ip=0; ip<np; ip++){
    double dx, dy, tmp1, tmp2, tmp3;
    dx = xp[ip] - xm;
    dy = yp[ip] - ym;
    tmp1 = dx * dx + dy * dy - dd * dd;
    tmp2 = dx * nx + dy * ny;
    a1 += tmp1 * tmp2;
    a2 += tmp2 * tmp2;
    a3 += (tmp2 * tmp2 / (dx * dx + dy * dy + 1.e-32));
    // printf("xp:%g, yp:%g\n", xp[ip], yp[ip]);
    // printf("dx:%g, dy:%g\n", dx, dy);
    // printf("tmp1:%g, tmp2:%g\n", tmp1, tmp2);
  }
  // printf("a1:%g, a2:%g, a3:%g\n", a1, a2, a3);
  // all points are colinear
  if (a3 < 1.e-32)
    return -1;
  
  a3 = 0.5 * a1 / a2;
  *xc = xm + a3 * nx;
  *yc = ym + a3 * ny;
  *r = sqrt(a3 * a3 + dd * dd);
  return 0;
}


double intersection_circle(double xc, double yc, double r, double x0){
  double delta, y0, y1, y2;

  delta = sq(r) - sq(x0 - xc);

  if (delta < 0.){
    return -1.;
  }
  else{
    delta = sqrt(delta);
    y1 = yc + delta;
    y2 = yc - delta;
    
    y0 = -1.;
    if (y1 > 0. && y1 < 1.){
      y0 = y1;
    }
    else if(y2 > 0. && y2 < 1.){
      y0 = y2;
    }
    return y0;
  }
}

// Initialize a circular interface eaxctly, we need the intersection points on cell edge.
// Which function in vofi should I use to replace this part?
double circle_sign(double x, double y, double xc, double yc, double ra){
  return sq(x - xc) + sq(y - yc) - sq(ra);
}

void init_circle(double xc, double yc, double ra, scalar f, face vector fs){
  double dir[2] = {0., 1.};
  foreach_dimension(){
    foreach_face(x){
      double x1, y1, x2, y2, x3, y3, eps, s1, s2;
      x1 = x - 0.5 * Delta * dir[0];
      y1 = y - 0.5 * Delta * dir[1];

      x2 = x + 0.5 * Delta * dir[0];
      y2 = y + 0.5 * Delta * dir[1];

      s1 = circle_sign(x1, y1, xc, yc, ra);
      s2 = circle_sign(x2, y2, xc, yc, ra);
      if (s1 * s2 < 0.){
        eps = 1.e-16 * Delta;
        for (int it=0; it<=40; it++){
          x3 = 0.5 * (x1 + x2);
          y3 = 0.5 * (y1 + y2);
          if (x3 - x1 < eps && y3 - y1 < eps)
            break;
          s1 = circle_sign(x1, y1, xc, yc, ra);
          s2 = circle_sign(x3, y3, xc, yc, ra);
          if (s1 * s2 < 0.){
            x2 = x3; y2 = y3;
          }
          else{
            x1 = x3; y1 = y3;
          }
        }
        fs.x[] = ((x3 - x) * dir[0]+ (y3 - y) * dir[1]) / Delta + 0.5;
        // printf("x:%g, y:%g, dxp:%g, dyp:%g\n", x / Delta, y / Delta, (x2 - x1) / Delta, (y2 - y1) / Delta);
      }
    }
    dir[0] = 1. - dir[0];
    dir[1] = 1. - dir[1];
  }
}

#if dimension <= 2
void init_circle_analytical(double xc, double yc, double ra, scalar f, face vector fs){
  double dir[2] = {0., 1.};
  struct { double x, y; } cc = {xc, yc};
  foreach_dimension(){
    foreach_face(x){
      double x1, y1, x2, y2, x3, y3, s1, s2;
      struct { double x, y; } xo = {x, y};
      s1 = fs.x[];
      x1 = sq(ra) - sq(xo.x - cc.x);
      x2 = x1;
      y3 = xo.y - 0.5 * Delta;
      if (x1 >= 0.) {
        x1 = sqrt(fabs(x1));
        y1 = (cc.y + x1 - y3) / Delta;
        y2 = (cc.y - x1 - y3) / Delta;
        fs.x[] = (y1 >= 0. && y1 <= 1.) ? y1 : ((y2 >= 0. && y2 <= 1.) ? y2 : 0.);
      }
    }
  }
}
#endif