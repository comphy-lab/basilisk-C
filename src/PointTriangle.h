#define vecdot(a,b) ((a).x*(b).x + (a).y*(b).y + (a).z*(b).z)
#define vecdotproduct(a,b) ((coord){(a).y*(b).z - (a).z*(b).y,	\
                                    (a).z*(b).x - (a).x*(b).z,	\
	                            (a).x*(b).y - (a).y*(b).x})
#define vecdiff(a,b)  ((coord){(a).x - (b).x, (a).y - (b).y, (a).z - (b).z})
#define vecdist2(a,b) (sq((a).x - (b).x) + sq((a).y - (b).y) + sq((a).z - (b).z))

/**
 * Returns the squared distance between point P and triangle P0, P1,
 * P2. Squared distance is returned in d2.  s and t returns the
 * closest point in parametric form in terms of edges P0P1 and P0P2.
 */
double PointTriangleDistance (const coord * P,
			      const coord * P0,
			      const coord * P1,
			      const coord * P2,
			      double * s, double * t)
{
  // From http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
  double d2;
  coord diff = vecdiff (*P0, *P);
  coord edge0 = vecdiff(*P1, *P0);
  coord edge1 = vecdiff(*P2, *P0);
  double a00 = vecdot(edge0, edge0);
  double a01 = vecdot(edge0, edge1);
  double a11 = vecdot(edge1, edge1);
  double b0 = vecdot(diff, edge0);
  double b1 = vecdot(diff, edge1);
  double c = vecdot(diff, diff);
  double det = fabs(a00*a11 - a01*a01);
  *s = a01*b1 - a11*b0;
  *t = a01*b0 - a00*b1;
  
  if (*s + *t <= det)
  {
    if (*s < 0.0)
    {
      if (*t < 0.0)  // region 4
      {
        if (b0 < 0.0)
        {
          *t = 0.0;
          if (-b0 >= a00)
          {
            *s = 1.0;
            d2 = a00 + 2.0*b0 + c;
          }
          else
          {
            *s = -b0/a00;
            d2 = b0**s + c;
          }
        }
        else
        {
          *s = 0.0;
          if (b1 >= 0.0)
          {
            *t = 0.0;
            d2 = c;
          }
          else if (-b1 >= a11)
          {
            *t = 1.0;
            d2 = a11 + 2.0*b1 + c;
          }
          else
          {
            *t = -b1/a11;
            d2 = b1**t + c;
          }
        }
      }
      else  // region 3
      {
        *s = 0.0;
        if (b1 >= 0.0)
        {
          *t = 0.0;
          d2 = c;
        }
        else if (-b1 >= a11)
        {
          *t = 1.0;
          d2 = a11 + 2.0*b1 + c;
        }
        else
        {
          *t = -b1/a11;
          d2 = b1**t + c;
        }
      }
    }
    else if (*t < 0.0)  // region 5
    {
      *t = 0.0;
      if (b0 >= 0.0)
      {
        *s = 0.0;
        d2 = c;
      }
      else if (-b0 >= a00)
      {
        *s = 1.0;
        d2 = a00 + 2.0*b0 + c;
      }
      else
      {
        *s = -b0/a00;
        d2 = b0**s + c;
      }
    }
    else  // region 0
    {
      if (det == 0.) // degenerate triangle
	d2 = HUGE;
      else {
	double invDet = 1.0/det;
	*s *= invDet;
	*t *= invDet;
	d2 = *s*(a00**s + a01**t + 2.0*b0) + *t*(a01**s + a11**t + 2.0*b1) + c;
      }
    }
  }
  else
  {
    if (*s < 0.0)  // region 2
    {
      double tmp0 = a01 + b0;
      double tmp1 = a11 + b1;
      if (tmp1 > tmp0)
      {
        double numer = tmp1 - tmp0;
        double denom = a00 - 2.0*a01 + a11;
        if (numer >= denom)
        {
          *s = 1.0;
          *t = 0.0;
          d2 = a00 + 2.0*b0 + c;
        }
        else
        {
          *s = numer/denom;
          *t = 1.0 - *s;
          d2 = *s*(a00**s + a01**t + 2.0*b0) + *t*(a01**s + a11**t + 2.0*b1) + c;
        }
      }
      else
      {
        *s = 0.0;
        if (tmp1 <= 0.0)
        {
          *t = 1.0;
          d2 = a11 + 2.0*b1 + c;
        }
        else if (b1 >= 0.0)
        {
          *t = 0.0;
          d2 = c;
        }
        else
        {
          *t = -b1/a11;
          d2 = b1**t + c;
        }
      }
    }
    else if (*t < 0.0)  // region 6
    {
      double tmp0 = a01 + b1;
      double tmp1 = a00 + b0;
      if (tmp1 > tmp0)
      {
        double numer = tmp1 - tmp0;
        double denom = a00 - 2.0*a01 + a11;
        if (numer >= denom)
        {
          *t = 1.0;
          *s = 0.0;
          d2 = a11 + 2.0*b1 + c;
        }
        else
        {
          *t = numer/denom;
          *s = 1.0 - *t;
          d2 = *s*(a00**s + a01**t + 2.0*b0) + *t*(a01**s + a11**t + 2.0*b1) + c;
        }
      }
      else
      {
        *t = 0.0;
        if (tmp1 <= 0.0)
        {
          *s = 1.0;
          d2 = a00 + 2.0*b0 + c;
        }
        else if (b0 >= 0.0)
        {
          *s = 0.0;
          d2 = c;
        }
        else
        {
          *s = -b0/a00;
          d2 = b0**s + c;
        }
      }
    }
    else  // region 1
    {
      double numer = a11 + b1 - a01 - b0;
      if (numer <= 0.0)
      {
        *s = 0.0;
        *t = 1.0;
        d2 = a11 + 2.0*b1 + c;
      }
      else
      {
        double denom = a00 - 2.0*a01 + a11;
        if (numer >= denom)
        {
          *s = 1.0;
          *t = 0.0;
          d2 = a00 + 2.0*b0 + c;
        }
        else
        {
          *s = numer/denom;
          *t = 1.0 - *s;
          d2 = *s*(a00**s + a01**t + 2.0*b0) + *t*(a01**s + a11**t + 2.0*b1) + c;
        }
      }
    }
  }
  
  // Account for numerical round-off error
  if (d2 < 0.0)
  {
    d2 = 0.0;
  }

  return d2;
}

int PointTriangleOrientation (const coord * P,
			      const coord * P0,
			      const coord * P1,
			      const coord * P2)
{
  coord diff = vecdiff (*P0, *P);
  coord edge0 = vecdiff (*P1, *P0);
  coord edge1 = vecdiff (*P2, *P0);
  coord n = vecdotproduct (edge0, edge1);
  return sign (vecdot (diff, n));
}

/**
   Returns the squared distance between p and [p0:p1].
*/
double PointSegmentDistance (const coord * p, const coord * p0, const coord * p1,
			     coord * segmentClosest, double * segmentParameter)
{
  // The direction vector is not unit length.  The normalization is deferred
  // until it is needed.
  coord direction = vecdiff (*p1, *p0);
  coord diff = vecdiff (*p, *p1);
  double t = vecdot(direction, diff);
  if (t >= (double)0)
    {
      *segmentParameter = (double)1;
      *segmentClosest = *p1;
    }
  else
    {
      diff = vecdiff (*p, *p0);
      t = vecdot(direction, diff);
      if (t <= (double)0)
        {
	  *segmentParameter = (double)0;
	  *segmentClosest = *p0;
        }
      else
        {
	  double sqrLength = vecdot(direction, direction);
	  if (sqrLength > (double)0)
            {
	      t /= sqrLength;
	      *segmentParameter = t;
	      (*segmentClosest).z = 0.;
	      foreach_dimension()
		(*segmentClosest).x = (*p0).x + t * direction.x;
            }
	  else
            {
	      *segmentParameter = (double)0;
	      *segmentClosest = *p0;
            }
        }
    }

  diff = vecdiff (*p, *segmentClosest);
  return vecdot(diff, diff);
}

int PointSegmentOrientation (const coord * P,
			     const coord * P0,
			     const coord * P1)
{
  coord diff = vecdiff (*P0, *P);
  coord edge0 = vecdiff (*P0, *P1);
  coord n = vecdotproduct (diff, edge0);
  return sign(n.z);
}
