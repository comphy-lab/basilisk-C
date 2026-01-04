/* A collection of Riemann solvers for the Saint-Venant system 
 *
 * References:
 *    [1] Kurganov, A., & Levy, D. (2002). Central-upwind schemes for the
 *    Saint-Venant system. Mathematical Modelling and Numerical
 *    Analysis, 36(3), 397-425.
 */

#define SQRT3 1.73205080756888
#define epsilon 1e-30

void kinetic (double hm, double hp, double um, double up, double Delta,
	      double * fh, double * fq, double * dtmax)
{
  double ci = sqrt(G*hm/2.);
  double Mp = max(um + ci*SQRT3, 0.);
  double Mm = max(um - ci*SQRT3, 0.);
  double cig = ci/(6.*G*SQRT3);
  *fh = cig*3.*(Mp*Mp - Mm*Mm);
  *fq = cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
  if (Mp > 0.) {
    double dt = CFL*Delta/Mp;
    if (dt < *dtmax)
      *dtmax = dt;
  }

  ci = sqrt(G*hp/2.);
  Mp = min(up + ci*SQRT3, 0.);
  Mm = min(up - ci*SQRT3, 0.);
  cig = ci/(6.*G*SQRT3);
  *fh += cig*3.*(Mp*Mp - Mm*Mm);
  *fq += cig*2.*(Mp*Mp*Mp - Mm*Mm*Mm);
  if (Mm < - epsilon) {
    double dt = CFL*Delta/-Mm;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

void kurganov (double hm, double hp, double um, double up, double Delta,
	       double * fh, double * fq, double * dtmax)
{
  double cp = sqrt(G*hp), cm = sqrt(G*hm);
  double ap = max(up + cp, um + cm); ap = max(ap, 0.);
  double am = min(up - cp, um - cm); am = min(am, 0.);
  double qm = hm*um, qp = hp*up;
  double a = max(ap, -am);
  if (a > epsilon) {
    *fh = (ap*qm - am*qp + ap*am*(hp - hm))/(ap - am); // (4.5) of [1]
    *fq = (ap*(qm*um + G*sq(hm)/2.) - am*(qp*up + G*sq(hp)/2.) + 
	    ap*am*(qp - qm))/(ap - am);
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
  else
    *fh = 0., *fq = 0.;
}

void hllc (double hm, double hp, double um, double up, double Delta,
	   double * fh, double * fq, double * dtmax)
{
  double cm = sqrt (G*hm), cp = sqrt (G*hp);
  double ustar = (um + up)/2. + cm - cp;
  double cstar = (cm + cp)/2. + (um - up)/4.;
  double SL = hm == 0. ? up - 2.*cp : min (um - cm, ustar - cstar);
  double SR = hp == 0. ? um + 2.*cm : max (up + cp, ustar + cstar);

  if (0. <= SL) {
    *fh = um*hm;
    *fq = hm*(um*um + G*hm/2.);
  }
  else if (0. >= SR) {
    *fh = up*hp;
    *fq = hp*(up*up + G*hp/2.);
  }
  else {
    double fhm = um*hm;
    double fum = hm*(um*um + G*hm/2.);
    double fhp = up*hp;
    double fup = hp*(up*up + G*hp/2.);
    *fh = (SR*fhm - SL*fhp + SL*SR*(hp - hm))/(SR - SL);
    *fq = (SR*fum - SL*fup + SL*SR*(hp*up - hm*um))/(SR - SL);
  }

  double a = max(fabs(SL), fabs(SR));
  if (a > epsilon) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}
