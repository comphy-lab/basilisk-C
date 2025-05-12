/** 
# Legendre Polynomials of Fractional Degrees

## General Definitions and Complete Elliptic Integrals

Associated *Legendre* polynomials of the first and second kind of degree $\nu$, 
respectively $P_\nu$ and $Q_\nu$, can be written with the help of 
*hypergeometric* functions of order $\mu = 0$ 
[ref. [DMLF (1)](https://dlmf.nist.gov/14.3)]:

$$
P_\nu(\cos \theta)
=
{}_2F_1 \left(\nu + 1, - \nu, 1, \dfrac{1 - \cos \theta}{2} \right)
$$

It can be shown that hypergeometric functions and 
[***complete elliptic integrals***](https://en.wikipedia.org/wiki/Elliptic_integral) 
of the first kind $K(m)$ and of the second kind $E(m)$ are related. 
Those last special functions have been implemented in `C programming language` 
by [*Eric Dennison*, (1998)](https://tiggerntatie.github.io/emagnet/offaxis/iloopoffaxis.htm), 
and his open-source code will be extremely useful to compute efficiently 
fractional *Legendre* polynomials:
*/


// ELLIPTIC.C:  Functions for computing complete and incomplete elliptic
//              integrals of the first and second kind, and an example
//              of usage for computing magnetic fields due to a circular
//              current filament.

// DISCLAIMER

// Although every effort has been expended to ensure correctness, this 
// software is not guaranteed to be correct.

// You may do anything you want with this software, as long as you maintain
// the existing comments intact, including the names of the original
//  authors.  If you modify this software, please indicate that fact in 
// your comments.

// E. Dennison (Nov 17, 1998)


#define DBL_EPSILON 2.2204460492503131e-016
#define DBL_MAX 	1.7976931348623158e+308
#define DBL_MIN 	2.2250738585072014e-308

#define MAX(x,y) (x>y?x:y)
#define MAX3(x,y,z) MAX(MAX(x,y),z)
#define MIN(x,y) (x>y?y:x)
#define MIN3(x,y,z) MIN(MIN(x,y),z)


// NOTE: constants declared in float.h for the x86 CPU are:

// DBL_EPSILON	2.2204460492503131e-016,smallest such that 1.0+DBL_EPSILON!=1.0
// DBL_MAX 	1.7976931348623158e+308  maximum possible value for type double
// DBL_MIN 	2.2250738585072014e-308  minimum possible positive value for double

// If you are compiling this code for a non-Intel CPU, your values for these
// constants may be different.





/* F(k): complete elliptic integral of the first kind */

#define F(k,ierr) drf(0.0,1.0-pow(k,2),1.0,&ierr)


// drf.c:   Compute the complete or incomplete elliptic integral of the
//          first kind.

// Description:

// For x,y and z non-negative and at most one of them zero, drf(x,y,z)
// = integral from zero to infinity of 

 
//           -1/2     -1/2     -1/2
// (1/2)(t+x)    (t+y)    (t+z)    dt.

// If x, y or z is zero, the integral is complete.

// *piErr returns non-zero if any arguments are invalid.

// Credits:

// This function is adapted by E. Dennison from FORTRAN code written by:

//   Carlson, B.C.
//   Notis, E.M
//   Pexton, R.L.



double drf(double x, double y, double z, int* piErr)
{
int    iErr=0;
double mu,
       xn,yn,zn,
       xndev,yndev,zndev,
       xnroot,ynroot,znroot,
       lambda,
       epsilon,
       e2,e3,
       result,
       s;
	
    const double c1 = 1.0/24.0;
    const double c2 = 3.0/44.0;
    const double c3 = 1.0/14.0;
    const double errtol = pow(DBL_EPSILON*4.0,1.0/6.0);
    const double lolim = 5.0*DBL_MIN;
    const double hilim = DBL_MAX/5.0;

    if (piErr)
    {
        if (MIN3(x,y,z)<0.0)
        {
            iErr = 1;
        }
        else if (MIN3(x+y,x+z,y+z)<lolim)
        {
            iErr = 2;
        }
        else if (MAX3(x,y,z)>hilim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            *piErr = iErr;
        }
        result = 0.0;
    }

else
    {
        xn = x;
        yn = y;
        zn = z;

        while (1)
        {
            mu = (xn+yn+zn)/3.0;
            xndev = 2.0-(mu+xn)/mu;
            yndev = 2.0-(mu+yn)/mu;
            zndev = 2.0-(mu+zn)/mu;
            epsilon = MAX3(fabs(xndev),fabs(yndev),fabs(zndev));
            if (epsilon < errtol) break;
            xnroot = sqrt(xn);
            ynroot = sqrt(yn);
            znroot = sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
        }
        e2 = xndev*yndev - pow(zndev,2);
        e3 = xndev*yndev*zndev;
        s = 1.0 + (c1*e2-0.1-c2*e3)*e2 + c3*e3;

        if (piErr)
        {
            *piErr = 0;
        }
        result = s/sqrt(mu);
    }
    return result;
}
/* END drf() */


/* E(k): complete elliptic integral of the second kind */


#define E(k,ierr)   (drf(0.0,1.0-pow(k,2),1.0,&ierr)\
                    -(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))


// FastE(k): fast, complete elliptic integral of the second kind.
//           Use this macro if the complete elliptic integral of
//           the first kind was previously computed for the same
//           value of k.


#define FastE(F,k,ierr)	((F)-(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))


// drd.c:   Compute the complete or incomplete elliptic integral of the
//          second kind.

// Description:

// For x and y non-negative, x+y and z positive, drf(x,y,z) = integral 
// from zero to infinity of 

//            -1/2     -1/2     -3/2
//  (3/2)(t+x)    (t+y)    (t+z)    dt.

// If x or y is zero, the integral is complete.

// *piErr returns non-zero if any arguments are invalid.

// Credits:

// This function is adapted by E. Dennison from FORTRAN code written by:

//   Carlson, B.C.
//   Notis, E.M
//   Pexton, R.L.



double drd(double x, double y, double z, int* piErr)
{
    int     iErr=0;
    double  mu,
            xn,yn,zn,
            xndev,yndev,zndev,
            xnroot,ynroot,znroot,
            lambda,
            epsilon,
            ea,eb,ec,ed,ef,
            sigma,
            power4,
            result,
            s1,s2;
	
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/6.0;
    const double c3 = 9.0/22.0;
    const double c4 = 3.0/26.0;
    const double errtol = pow(DBL_EPSILON/3.0,1.0/6.0);
    double uplim;
    const double lolim = 2.0/pow(DBL_MAX,2.0/3.0);
    double tuplim = pow(DBL_MIN,1.0/3.0);
    tuplim = pow(0.1*errtol,1.0/3.0)/tuplim;
    uplim = pow(tuplim,2.0);

    if (piErr)
    {
        if (MIN(x,y)<0.0)
        {
            iErr = 1;
        }
        else if (MAX3(x,y,z)>uplim)
        {
            iErr = 2;
        }
        else if (MIN(x+y,z)<lolim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            *piErr = iErr;
        }
        result = 0.0;
    }

else
    {
        xn = x;
        yn = y;
        zn = z;
        sigma = 0.0;
        power4 = 1.0;
        while (1)
        {
            mu = (xn+yn+3.0*zn)*0.2;
            xndev = (mu-xn)/mu;
            yndev = (mu-yn)/mu;
            zndev = (mu-zn)/mu;
            epsilon = MAX3(fabs(xndev),fabs(yndev),fabs(zndev));
            if (epsilon < errtol) break;
            xnroot = sqrt(xn);
            ynroot = sqrt(yn);
            znroot = sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            sigma = sigma+power4/(znroot*(zn+lambda));
            power4 = power4*0.25;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
        }
        ea = xndev*yndev;
        eb = zndev*zndev;
        ec = ea-eb;
        ed = ea-6.0*eb;
        ef = ed+ec+ec;
        s1 = ed*(-c1+0.25*c3*ed-1.5*c4*zndev*ef);
        s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea));
        if (piErr)
        {
            *piErr = 0;
        }
        result = 3.0*sigma+power4*(1.0+s1+s2)/(mu*sqrt(mu));
    }
    return result;
}
/* END drd() */


/**
Following the prescriptions of [Abramowitz \& Stegun, (1964)](#abramowitz1964) 
[formulae $8.13.8/10/11/12$], [DMLF (2)](https://dlmf.nist.gov/14.5#v), 
[WolframAlpha](https://www.wolframalpha.com/input?i=Hypergeometric2F1\%5B3\%2F2\%2C-1\%2F2\%2C1\%2C\%281-Cos\%5Bx\%5D\%29\%2F2\%5D) 
(with $(1-x)/2$ as an argument, **not** the convention $\sqrt{(1-x)/2}$) and 
[`maths.stackexchange.com`](https://math.stackexchange.com/questions/102255/associated-legendre-polynomials-of-fractional-order), 
we now have for $\nu = 1/2$:

$$
\begin{array}{lcl}
  P_{1/2}(\cos \theta) 
  &=& 
  \dfrac{2}{\pi} 
  \left(
    2 E \left(\sqrt{\dfrac{1 - \cos \theta}{2}} \right)
    - K \left(\sqrt{\dfrac{1 - \cos \theta}{2}} \right)
  \right) \\
  Q_{1/2}(\cos \theta) 
  &=&
  K \left(\sqrt{\dfrac{1 + \cos \theta}{2}} \right)
  - 2 E \left(\sqrt{\dfrac{1 + \cos \theta}{2}} \right)  
\end{array}
$$

and:

$$
P_{-1/2}(\cos \theta) 
= 
\dfrac{2}{\pi} 
K \left(\sqrt{\dfrac{1 - \cos \theta}{2}} \right)
\quad ; \quad 
Q_{-1/2}(\cos \theta) 
= 
K \left(\sqrt{\dfrac{1 + \cos \theta}{2}} \right)
$$

Using the complete elliptic integrals code from *Eric Dennison*, (1998), 
these expressions are then readily implemented:
*/

/* Computation of different Legendre polynomials of fractional degrees */

// Degree 1/2
// ----------

// P_{1/2}(\cos θ)
double legendreP_half (double theta)
{
  double k = sqrt((1. - cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double ek = FastE(fk,k,ierr);
  double P_half = (2./pi)*(2.*ek - fk);
  
  return P_half;
}

// Q_{1/2}(\cos θ)
double legendreQ_half (double theta)
{
  double k = sqrt((1. + cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double ek = FastE(fk,k,ierr);
  double Q_half = (fk - 2.*ek);
  
  return Q_half;
}

// Degree -1/2
// -----------

// P_{-1/2}(\cos θ)
double legendreP_minus_half (double theta)
{
  double k = sqrt((1. - cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double P_minus_half = (2./pi)*fk;
  
  return P_minus_half;
}

// Q_{-1/2}(\cos θ)
double legendreQ_minus_half (double theta)
{
  double k = sqrt((1. + cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double Q_minus_half = fk;
  
  return Q_minus_half;
}


/**
## Interesting Recurrence Relations

According to the formula $8.5.3$ of [Abramowitz \& Stegun, (1964)](#abramowitz1964), 
also the up-to-date [NIST version](https://univ.jeanpaulcalvi.com/Posters/ConfAuchWeb/abramovitz2.pdf)
$P_\nu$ and $Q_\nu$ polynomials verify the same recurrence relations. 
In particular, for $P_\nu$:

$$
P_{\nu+1}(x) = \dfrac{2\nu + 1}{\nu +1} x \, P_\nu(x) 
  - \dfrac{\nu}{\nu+1} P_{\nu-1}(x) 
\quad\quad (*)
$$

This last expression is of great interest, as it offers a simple mean to compute 
higher degrees, and because derivatives $\partial_\theta P_\nu(\cos \theta)$ are 
also computed thanks to these higher degrees $\nu + 1$. Indeed, the formula 
$8.5.4$ of [Abramowitz \& Stegun, (1964)](#abramowitz1964) gives:

$$
(x^2-1)\dfrac{\mathrm{d} P_{\nu}(x)}{\mathrm{d}x} 
= \nu \left[x \, P_\nu(x) - P_{\nu-1}(x)\right]
\quad\quad (**)
$$


## Some Examples

From $(*)$ and the general definitions given above of $P_{1/2}$, $Q_{1/2}$, 
$P_{-1/2}$, $Q_{-1/2}$, we can readily compute $P_{3/2}$ and $Q_{3/2}$:

$$
P_{3/2}(\cos \theta)
=
- \dfrac{2}{3 \pi}\left( 1 + 4 \cos\theta \right)
K\left(\sqrt{\dfrac{1-\cos\theta}{2}}\right)
+\dfrac{16}{3 \pi}\cos(\theta) E\left(\sqrt{\dfrac{1-\cos \theta}{2}}\right)
$$

and:

$$
Q_{3/2}(\cos \theta) = \dfrac{4}{3} \cos \theta \, 
\left(
K \left( \sqrt{\dfrac{1+ \cos \theta}{2}} \right)
- 2 E \left( \sqrt{\dfrac{1+\cos \theta}{2}} \right)
\right)
- \dfrac{1}{3}K \left( \sqrt{\dfrac{1+\cos \theta}{2}} \right)
$$
*/

// Degree 3/2
// ----------

// P_{3/2}(\cos θ)
double legendreP_three_half (double theta)
{
  double k = sqrt((1. - cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double ek = FastE(fk,k,ierr);
  double A = -(2./(3*pi))*(1. + 4.*cos(theta)) ;
  double B = (16./(3*pi))*cos(theta);
  double P_three_half = A*fk + B*ek;
  
  return P_three_half;
}

// Q_{3/2}(\cos θ)
double legendreQ_three_half (double theta)
{
  double k = sqrt((1. + cos(theta))/2.);
  int ierr;
  double fk = F(k,ierr);
  double ek = FastE(fk,k,ierr);
  double A = (4./3.)*cos(theta) ;
  double Q_three_half = A*(fk - 2.*ek) - (1./3.)*fk;
  
  return Q_three_half;
}


/** 
It is also possible to combine $(*)$ and $(**)$ to express certain derivatives 
only with *Legendre* polynomials of *positive* degree:

$$
\dfrac{\mathrm{d}P_{1/2}(\cos \theta)}{\mathrm{d}\theta}
=
\dfrac{3}{2}\left(
  -\cot \theta \, P_{1/2}(\cos \theta)
  + \csc \theta \, P_{3/2}(\cos \theta)
\right)
$$

Hence:
*/

// ∂P_{1/2}(\cos θ)/∂θ
double legendreP_prime_half (double theta)
{
  double A = (-3./(tan(theta))) * legendreP_half(theta);
  double B = (3./(sin(theta))) * legendreP_three_half(theta);

  return (0.5 * (A + B));
}


/**
## Other useful functions implemented
*/

// Degree 5/2 (with recursion formula)
// ----------

// P_{5/2}(\cos θ)
double legendreP_five_half (double theta)
{
  double P_half = legendreP_half(theta) ;
  double P_three_half = legendreP_three_half(theta);
  double P_five_half = (8./5.)*cos(theta)*P_three_half - (3./5.)*P_half;
  
  return P_five_half;
}

// Q_{5/2}(\cos θ)
double legendreQ_five_half (double theta)
{
  double Q_half = legendreQ_half(theta) ;
  double Q_three_half = legendreQ_three_half(theta);
  double Q_five_half = (8./5.)*cos(theta)*Q_three_half - (3./5.)*Q_half;
  
  return Q_five_half;
}

/* ========================================================================== */

/* Computation of some Legendre Derivatives */

// ∂P_{1/2}(\cos π-θ)/∂θ
double legendreP_prime_half_minus_pi (double theta)
{
  double A = -(3./2.)/(sin(theta));
  double B = cos(theta)*legendreP_half(pi - theta);
  double C = legendreP_three_half(pi - theta);

  return (C * (A + B));
}

// ∂Q_{1/2}(\cos θ)/∂θ
double legendreQ_prime_half (double theta)
{
  double A = (-3./(tan(theta))) * legendreQ_half(theta);
  double B = (3./(sin(theta))) * legendreQ_three_half(theta);

  return (0.5 * (A + B));
}

// ∂²P_{1/2}(\cos θ)/∂θ²
double legendreP_2prime_half (double theta)
{
  double A = (3./8.) * sq(1/sin(theta));
  double B = (7. + 3.*cos(2.*theta)) * legendreP_half(theta) ;
  double C = 10.*(-2.*cos(theta)*legendreP_three_half(theta) 
                  + legendreP_five_half(theta));

  return (A * (B + C));
}

// ∂²Q_{1/2}(\cos θ)/∂θ²
double legendreQ_2prime_half (double theta)
{
  double A = (3./8.) * sq(1/sin(theta));
  double B = (7. + 3.*cos(2.*theta)) * legendreQ_half(theta) ;
  double C = 10.*(-2.*cos(theta)*legendreQ_three_half(theta) 
                  + legendreQ_five_half(theta));

  return (A * (B + C));
}



/** 
## APPLICATION: *Sierou \& Lister*, (2004) -- Appendix A

*Legendre* polynomials of fractional degrees are of the utmost importance 
for surface tension driven recoiling cones of liquid in the presence of a 
*far-field dipolar flow*; see [Sierou \& Lister, (2004)](#sierou2004), Appendix 
A, for an explanation of the notations used in the code hereafter.
*/



/* Computation of Zero-Order Coefficients for the velocity potentials */
// Cf. 'Appendix A' from Sierou & Lister (2004)

// For domains A & O (if working with a cone recoiling towards NEGATIVE `z`)
// -------------------------------------------------------------------------

double phi_0_A (double th_alph, double mu_0)
{
  double num = cos(th_alph)*legendreQ_half(th_alph) - legendreQ_three_half(th_alph);
  double denomA = legendreP_three_half(th_alph)*legendreQ_half(th_alph);
  double denomB = legendreP_half(th_alph)*legendreQ_three_half(th_alph);

  return ( (num / (denomA - denomB)) * mu_0 );
}


double phi_0_OQ (double th_alph, double mu_0)
{
  double num = cos(th_alph)*legendreP_half(th_alph) - legendreP_three_half(th_alph);
  double denomA = legendreP_three_half(th_alph)*legendreQ_half(th_alph);
  double denomB = legendreP_half(th_alph)*legendreQ_three_half(th_alph);

  return ( (num / (denomA - denomB)) * mu_0 );
}


/* Computation of First-Order Coeffs for the shape and velocity potentials */

double theta_1 (double th_alph, double phi_0)
{
  return (- phi_0 * legendreP_prime_half(th_alph)) ;
}


double phi_1_A (double th_alph, double phi_0){
  // post-pinchoff formulation (therefore, the minus sign)
  return (
    - (
      1./tan(th_alph) 
    + sq(phi_0)/(8.*sq(sin(th_alph)))*(
      (4.*cos(2*th_alph) + 5)*sq(legendreP_half(th_alph))
      - 18.*cos(th_alph)*legendreP_half(th_alph)*legendreP_three_half(th_alph)
      + 9. * sq(legendreP_three_half(th_alph))
      )
    )
  );
}


// For domains B & O (if working with a cone recoiling towards POSITIVE `z`)
// -------------------------------------------------------------------------

double phi_0_B (double th_bet, double mu_0)
{
  double num = legendreP_prime_half(th_bet);
  double denomA = legendreP_half(pi-th_bet)*legendreP_prime_half(th_bet);
  double denomB = legendreP_half(th_bet)*legendreP_prime_half_minus_pi(th_bet);

  return ( (num / (denomA - denomB)) * mu_0 );
}


double phi_0_OP (double th_bet, double mu_0)
{
  double num = legendreP_prime_half_minus_pi(th_bet);
  double denomA = legendreP_half(pi-th_bet)*legendreP_prime_half(th_bet);
  double denomB = legendreP_half(th_bet)*legendreP_prime_half_minus_pi(th_bet);

  return ( (num / (denomA - denomB)) * mu_0 );
}



/* ========================================================================== */
/* Computation of O(1/R) coeff. for the self-similar gradient of pressure */

// ∂p/∂R
double fpR (double R, double th, double th_alph, double phi_0){
  double first_term = phi_1_A (th_alph, phi_0);
  double prefactor = sq(phi_0)/8.;
  double second_term = sq(legendreP_half(th)) 
    + 4.*sq(legendreP_prime_half(th));
  return (
    -(1/sq(R))*(first_term - prefactor*second_term)
  );
}

// ∂p/∂θ
double fpth (double R, double th, double th_alph, double phi_0){
  double first_term = 2 * legendreP_half(th) *legendreP_prime_half(th);
  double prefactor = -sq(phi_0)/(8. * sq(R));
  double second_term = 8.*legendreP_prime_half(th)*legendreP_2prime_half(th);
  return (
    prefactor*(first_term + second_term)
  );
}

// Projection onto e_z:
double fpz (double R, double th, double th_alph, double phi_0){
  double fR = fpR (R, th, th_alph, phi_0);
  double fth = fpth (R, th, th_alph, phi_0);
  return(
    fR*cos(th) - fth*sin(th)
  );
}

// Projection onto e_r:
double fpr (double R, double th, double th_alph, double phi_0){
  double fR = fpR (R, th, th_alph, phi_0);
  double fth = fpth (R, th, th_alph, phi_0);
  return(
    fR*sin(th) + fth*cos(th)
  );
}

/* ========================================================================== */
/* Computation of the self-similar velocity */

// ∂φ/∂R 
// (interior)
double uR_0_A (double R, double th, double phi_0){
  return (
    (phi_0/2.)*legendreP_half(th)/sqrt(R)
  );
}

double uR_1_A (double R, double th, double th_alph, double phi_0){
  return (
      - phi_1_A (th_alph, phi_0)/sq(R)
  );
}

double uR_0_B (double R, double th, double phi_0){
  return (
    (phi_0/2.)*legendreP_half(pi - th)/sqrt(R)
  );
}


// (exterior)
double uR_0_OQ (double R, double th, double phi_0){
  return (
    (phi_0/2.)*legendreQ_half(th)/sqrt(R)
  );
}

double uR_0_OP (double R, double th, double phi_0){
  return (
    (phi_0/2.)*legendreP_half(th)/sqrt(R)
  );
}

double uR_1_OP (double R, double th, double th_alph, double phi_0){
  return (
      - uR_1_A (R, th, th_alph, phi_0)
  );
}



// ∂φ/∂θ 
//(interior)
double uth_0_A (double R, double th, double phi_0){
  return (
    phi_0*legendreP_prime_half(th)/sqrt(R)
  );
}

double uth_0_B (double R, double th, double phi_0){
  return (
    phi_0*legendreP_prime_half_minus_pi(th)/sqrt(R)
  );
}

//(exterior)
double uth_0_OQ (double R, double th, double phi_0){
  return (
    phi_0*legendreQ_prime_half(th)/sqrt(R)
  );
}

double uth_0_OP (double R, double th, double phi_0){
  return (
    phi_0*legendreP_prime_half(th)/sqrt(R)
  );
}

/** Finally, we need to project fields defined in spherical-polar 
coordinates onto the cylindrical coordinate system: */

// Projection onto e_z:
double proj_ez (double fR, double fth, double th){
  return (
    fR*cos(th) - fth*sin(th)
  );
}

// Projection onto e_r:
double proj_er (double fR, double fth, double th){
  return (
    fR*sin(th) + fth*cos(th)
  );
}




/**
## References

~~~bib
@book{abramowitz1964,
  address = {New York},
  author = {Abramowitz, Milton and Stegun, Irene A.},
  edition = {ninth Dover printing, tenth GPO printing},
  interhash = {d4914a420f489f7c5129ed01ec3cf80c},
  intrahash = {23ec744709b3a776a1af0a3fd65cd09f},
  publisher = {Dover},
  title = {Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables},
  URL ={https://personal.math.ubc.ca/~cbm/aands/abramowitz_and_stegun.pdf},
  year = {1964}
}

@article{sierou2004,
  author = {Sierou, A. and Lister, J. R.},
  title = {Self-similar recoil of inviscid drops},
  journal = {Physics of Fluids},
  volume = {16},
  number = {5},
  pages = {1379-1394},
  year = {2004},
  month = {05},
  issn = {1070-6631},
  doi = {10.1063/1.1689031},
  url = {https://doi.org/10.1063/1.1689031},
  eprint = {https://pubs.aip.org/aip/pof/article-pdf/16/5/1379/19153172/1379\_1\_online.pdf},
}
~~~
*/







