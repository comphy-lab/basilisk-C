/**
# Distance from a point to an ellipsoid

See section 3.10 of
[DistancePointEllipseEllipsoid.pdf](https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf). 

Here, the result is the signed distance.
*/



#include "../popinet/distance_point_ellipse.h"

static double RobustLength3D (double v0, double v1, double v2)
{
  double vk = max( max(fabs(v0), fabs(v1)) , fabs(v2) );
  return vk*sqrt(sq(v0/vk)+sq(v1/vk)+sq(v2/vk));
}

static double GetRoot3D(double r0, double r1, 
  double z0, double z1, double z2, double g){
  double n0 = r0*z0, n1 = r1*z1;
  double s0 = z2 - 1., s1 =  (g < 0 ? 0 : RobustLength3D(n0,n1,z2 -1));
  double s = 0.;
  for (int i = 0; i < 149; i++){
      s = (s0 + s1)/2.;
      if(s == s0 || s == s1) break;
      double ratio0 = n0 / (s + r0), ratio1 = n1/(s + r1), ratio2 = z2/(s +1.);
      g = sq(ratio0) + sq(ratio1) + sq(ratio2) -1.;
      if(g>0) s0 = s;
      else if (g <0) s1 = s;
      else break;
  }
  return s;
}

double DistancePointEllipsoid(double e0, double e1, double e2, double y0,
  double y1, double y2, double * x0, double * x1, double * x2){
/*
Calculations must be done in the first quadrant
*/
bool sym0 = false, sym1 = false, sym2 = false;
  if (y0 < 0.)
    y0 = - y0, sym0 = true;
  if (y1 < 0.)
    y1 = - y1, sym1 = true;
  if (y2 < 0.)
    y2 = - y2, sym2 = true;

  double distance;
  if (y2 > 0.) {
    if (y1 > 0.) {
      if (y0 > 0.) {
        double z0 = y0/e0, z1 = y1/e1, z2 = y2/e2;
        double g = sq(z0) + sq(z1) + sq(z2) - 1.;
        if (g != 0.) {
          double r0 = sq(e0/e2), r1 = sq(e1/e2);
          double sbar = GetRoot3D (r0, r1, z0, z1, z2, g);
          *x0 = r0*y0/(sbar + r0);
          *x1 = r1*y1/(sbar + r1);
          *x2 = y2/(sbar + 1.);
          distance = sqrt (sq(*x0 - y0) + sq(*x1 - y1)+sq(*x2 - y2));
        }
        else {
          *x0 = y0;
          *x1 = y1;
          *x2 = y2;
          distance = 0.;
        }
      }
      else {
        // y0 == 0
        *x0 = 0.;
        distance = DistancePointEllipse(e1, e2, y1, y2, x1, x2);
      }
    }
    else{ 
      // y1 ==0
      if(y0>0){
        *x1 = 0;
        distance = DistancePointEllipse(e0, e2, y0, y2, x0, x2);
      }
      else{
        // y0 == 0
        *x0 = 0, *x1 = 0, *x2 = e2;
        distance = fabs(y2 - e2);
      }
    }
  }
  else{
    // y2 ==0
    double denom0 = sq(e0) - sq(e2), denom1 = sq(e1) - sq(e2);
    double numer0 = e0*y0, numer1 = e1*y1;
    bool computed = false;
    if(numer0 < denom0 && numer1 < denom1){
      double xde0   = numer0/denom0 , xde1   = numer1/denom1 ;
      double xde0sq = xde0*xde0     , xde1sq = xde1*xde1  ;  
      double discr  = 1. - xde0sq -xde1sq;
      if( discr > 0){
        *x0 = e0*xde0;
        *x1 = e1*xde1;
        *x2 = e2*sqrt(discr);
        distance = sqrt(sq(*x0 - y0) + sq(*x1 - y1) + sq(*x2));
        computed = true;
      }
    }
    if(!computed){
      *x2 = 0.; 
      distance = DistancePointEllipse(e0, e1, y0, y1, x0, x1);
    }
  }

  if (sym0) *x0 = - *x0;
  if (sym1) *x1 = - *x1;
  if (sym2) *x2 = - *x2;
  return sign(sq(y0/e0) + sq(y1/e1) + sq(y2/e2) - 1.)*distance;
}