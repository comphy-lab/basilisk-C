/* Implementation of the formulae of Okada, 1985, "Surface deformation
   due to shear and tensile faults in a half-space", Bulletin of the
   Seismological Society of America, 75:4, 1135-1154, */

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
  double Reta = fabs (R + eta) > 1e-6 ? R + eta : 1e30;
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

struct Okada {
  scalar d;
  double x, y, depth;
  double strike, dip, rake;
  double mu, lambda;
  double length, width, vU[3], U;
  double R;
  int (* iterate) (void);
  bool flat;
};

void okada (struct Okada p)
{
  // default settings
  if (p.mu == 0.)     p.mu = 1.;
  if (p.lambda == 0.) p.lambda = 1.;
  if (p.R == 0.)      p.R = 6371220.; /* Earth radius (metres) */

  double dtr = pi/180.;
  if (p.rake != nodata) {
    p.vU[0] = p.U*cos (p.rake*dtr);
    p.vU[1] = p.U*sin (p.rake*dtr);
  }
  double sina = sin ((90. - p.strike)*dtr);
  double cosa = cos ((90. - p.strike)*dtr);
  double sind = sin (p.dip*dtr);
  /* depth of the bottom edge */
  double depth = sind > 0. ? p.depth + p.width*sind : p.depth;
  /* origin to the centroid */
  double x0 = p.length/2., y0 = p.width/2.*cos (p.dip*dtr);
  
  foreach() {
    if ( p.flat ) {
      x -= p.x;
      y -= p.y;
    }
    else {
      x = p.R*cos(y*dtr)*dtheta(x, p.x)*dtr;
      y = p.R*dtheta(y, p.y)*dtr;
    }
    double x1 =   cosa*x + sina*y;
    double y1 = - sina*x + cosa*y;
    double oka[3];
    okada_rectangular_source (p.vU, p.length, p.width, depth, 
			      p.dip*dtr,
			      p.mu/p.lambda,
			      x0 + x1, y0 + y1,
			      oka);
    val(p.d,0,0) = oka[2];
  }
}

void fault (struct Okada p)
{
  scalar hold[];
  // save the initial water depth
  scalar_clone (hold, h);
  foreach()
    hold[] = h[];
  boundary ({hold});

  p.d = h;
  do {
    okada (p);
    // h[] now contains the Okada vertical displacement
    foreach() {
      // deformation is added to hold[] (water depth) only in wet areas
      h[] = hold[] > dry ? max (0., hold[] + h[]) : hold[];
      eta[] = zb[] + h[];
    }
  } while (p.iterate && p.iterate());
}
