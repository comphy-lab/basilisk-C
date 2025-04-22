/**
# Okada fault model

This is an implementation of the formulae of [Okada, 1985](#okada1985). */

/* formulae (25)-(30) */
static void rectangular_source (const double U[3], double cosd, double sind,
				double mulambda, double d,
				double psi, double eta, double q,
				double u[3])
{
  double R = sqrt (psi*psi + eta*eta + q*q);
  double X = sqrt (psi*psi + q*q);
  double dtilde = eta*sind - q*cosd;
  double ytilde = eta*cosd + q*sind;
  double atanp = fabs (q) > 1e-6 ? atan (psi*eta/(q*R)) : 0.;

  mulambda = mulambda/(1. + mulambda);
  double logReta = R + eta > 1e-6 ? log (R + eta) : - log (R - eta);
  double Reta = fabs (R + eta) > 1e-6 ? R + eta : HUGE;
  double I1, I2, I3, I4, I5;
  if (fabs (cosd) > 1e-6) {
    /* formula (28) */
    I5 = fabs (psi) < 1e-6 ? 0. :
      mulambda*2./cosd*atan ((eta*(X + q*cosd) + 
			      X*(R + X)*sind)/(psi*(R + X)*cosd));
    I4 = mulambda/cosd*(log (R + dtilde) - sind*logReta);
    I3 = mulambda*(1./cosd*ytilde/(R + dtilde) - logReta) + sind/cosd*I4;
    I2 = mulambda*(- logReta) - I3;
    I1 = mulambda*(-1./cosd*psi/(R + dtilde)) - sind/cosd*I5;
  }
  else {
    /* formula (29) */
    double R1 = R + dtilde;
    I1 = - mulambda/2.*psi*q/(R1*R1);
    I3 = mulambda/2.*(eta/R1 + ytilde*q/(R1*R1) - logReta);
    I2 = mulambda*(- logReta) - I3;
    I4 = - mulambda*q/R1;
    I5 = - mulambda*psi*sind/R1;
  }
    
  /* strike-slip, formula (25) */  
  if (U[0] != 0.) {
    double U1pi = U[0]/(2.*M_PI);
    u[0] -= U1pi*(psi*q/(R*Reta) + atanp + I1*sind);
    u[1] -= U1pi*(ytilde*q/(R*Reta) + q*cosd/Reta + I2*sind);
    u[2] -= U1pi*(dtilde*q/(R*Reta) + q*sind/Reta + I4*sind);
  }

  /* dip-slip, formula (26) */  
  if (U[1] != 0.) {
    double U2pi = U[1]/(2.*M_PI);
    u[0] -= U2pi*(q/R - I3*sind*cosd);
    u[1] -= U2pi*(ytilde*q/(R*(R + psi)) + cosd*atanp - I1*sind*cosd);
    u[2] -= U2pi*(dtilde*q/(R*(R + psi)) + sind*atanp - I5*sind*cosd);
  }

  /* tensile, formula (27) */  
  if (U[2] != 0.) {
    double U3pi = U[2]/(2.*M_PI);
    u[0] += U3pi*(q*q/(R*Reta) - I3*sind*sind);
    u[1] += U3pi*(-dtilde*q/(R*(R + psi)) - 
		  sind*(psi*q/(R*Reta) - atanp) - I1*sind*sind);
    u[2] += U3pi*(ytilde*q/(R*(R + psi)) + 
		  cosd*(psi*q/(R*Reta) - atanp) - I5*sind*sind);
  }
}

/* formula (24) */
static void okada_rectangular_source (const double U[3], 
				      double L, double W, double d, 
				      double delta, double mulambda,
				      double x, double y,
				      double u[3])
{
  double cosd = cos (delta), sind = sin (delta);
  double p = y*cosd + d*sind;
  double q = y*sind - d*cosd;

  /**
  There seems to be a problem with dimensions here. `x`, `p`, and `q`
  should be the dimensionless coordinates, not the dimensional
  ones... See also [tsunami.c](/src/examples/tsunami.c). */
  
  u[0] = u[1] = u[2] = 0.;
  rectangular_source (U, cosd, sind, mulambda, d,
		      x, p, q,
		      u);
  rectangular_source (U, cosd, sind, mulambda, d,
		      x - L, p - W, q,
		      u);

  double u1[3] = {0., 0., 0.};
  rectangular_source (U, cosd, sind, mulambda, d,
		      x, p - W, q,
		      u1);
  rectangular_source (U, cosd, sind, mulambda, d,
		      x - L, p, q,
		      u1);
  u[0] -= u1[0];
  u[1] -= u1[1];
  u[2] -= u1[2];
}

static double dtheta (double theta1, double theta2)
{
  double d = theta1 - theta2;
  if (d > 180.) d -= 360.;
  if (d < -180.) d += 360.;
  return d;
}

typedef struct {
  double depth, x, y;
  double strike, dip, rake;
  double length, width, U;
  double vU[3];
} Fault;

void okada (scalar d,
	    double x = 0, double y = 0, double depth = 0,
	    double strike = 0, double dip = 0, double rake = 0,
	    double length = 0, double width = 0, double U = 0,
	    double mu = 1, double lambda = 1,
	    double R = 6371220., /* Earth radius (metres) */
	    bool flat = false, bool centroid = false,
	    Fault * faults = NULL)
{
  Fault lfaults[2] = {0};
  if (faults == NULL) {
    lfaults[0] = (Fault) {depth, x, y,
                          strike, dip, rake,
                          length, width, U};
    faults = lfaults;
  }
  foreach (cpu) {
    d[] = 0.;
    for (Fault * f = faults; f && f->depth > 0.; f++) {
      double depth = f->depth, dtr = pi/180.;
      if (centroid)
	depth -= f->width*fabs (sin (f->dip*dtr))/2.;
      if (f->rake != nodata) {
	f->vU[0] = f->U*cos (f->rake*dtr);
	f->vU[1] = f->U*sin (f->rake*dtr);
      }
      double sina = sin ((90. - f->strike)*dtr);
      double cosa = cos ((90. - f->strike)*dtr);
      double sind = sin (f->dip*dtr);
      /* depth of Okada origin */
      depth = sind > 0. ? depth + f->width*sind : depth;
      /* origin to the centroid */
      double x0 = f->length/2., y0 = f->width/2.*cos (f->dip*dtr);

      double xc, yc;
      if (flat) {
	xc = x - f->x;
	yc = y - f->y;
      }
      else {
	xc = R*cos(y*dtr)*dtheta(x, f->x)*dtr;
	yc = R*dtheta(y, f->y)*dtr;
      }
      double x1 =   cosa*xc + sina*yc;
      double y1 = - sina*xc + cosa*yc;
      double oka[3];
      okada_rectangular_source (f->vU, f->length, f->width, depth, 
				f->dip*dtr,
				mu/lambda,
				x0 + x1, y0 + y1,
				oka);
      d[] += oka[2];
    }
  }
}

/**
## User interface

Use function fault() to alter water depth *h* (where *h > dry*) according
to the fault parameters:

* *x, y*: coordinates of the fault centroid (see boolean *flat* for
   coordinate type).
* *depth*: depth of the top edge of the fault (see also *centroid*).
* *strike, dip, rake*: fault parameters in degrees.  (0 <= *strike* <
  360, -90 <= *dip* <= 90, -90 <= *rake* <= 90 where *rake* = 90 degs
  and *dip* > 0 is reverse faulting. NB: Okada defines normal faulting
  by *rake* = 90 deg and *dip* < 0 whereas the seismological
  convention now is generally 0 <= *dip* <= 90 and *rake* = -90 for
  normal faulting).
* *mu, lambda*: only the ratio is important and default is  *mu/lambda* = 1.
* *length, width, U*:  length and width of the fault plane and slip on
  the fault plane (generally in meters).
* *R*: is the radius of the earth (for when x, y are in longitude and
  latitude i.e. *flat* = false).
* *iterate*: is the function to use to iterate.
* *flat*: true assumes x, y cartesian, false assumes longitude and
  latitude (default).
* *centroid*: assumes that depth is measured to the centroid of the
  fault, not the top edge.
* *faults*: if non-NULL, this defines an array of several fault
  parameters which will be used instead of the parameters for a single
  fault above. This is useful to efficiently define a deformation
  composed of many Okada subfaults. Note that the array must be
  terminated by a "dummy fault" of depth smaller than or equal to zero. 
*/

void fault (double x = 0, double y = 0, double depth = 0,
	    double strike = 0, double dip = 0, double rake = 0,
	    double length = 0, double width = 0, double U = 0,
	    double mu = 1, double lambda = 1,
	    double R = 6371220., /* Earth radius (metres) */
	    int (* iterate) (void) = NULL,
	    bool flat = false, bool centroid = false,
	    Fault * faults = NULL)
{
  scalar hold[];
  // save the initial water depth
  scalar_clone (hold, h);
  foreach()
    hold[] = h[];

  int nitermax = 20;
  do {
    okada (h, x, y, depth, strike, dip, rake, length, width, U,
	   mu, lambda, R, flat, centroid, faults);
    // h[] now contains the Okada vertical displacement
    foreach() {
      // deformation is added to hold[] (water depth) only in wet areas
      h[] = hold[] > dry ? max (0., hold[] + h[]) : hold[];
      eta[] = zb[] + h[];
    }
  } while (iterate && (* iterate) () && nitermax--);
}

/**
## References

~~~bib
@article{okada1985,
  title={Surface deformation due to shear and tensile faults in a half-space},
  author={Okada, Yoshimitsu},
  journal={Bulletin of the seismological society of America},
  volume={75},
  number={4},
  pages={1135--1154},
  year={1985},
  publisher={The Seismological Society of America}
}
~~~
*/
