/**
# TO DO

* The gradient use for high-order reconstruction are not defined in a well-balanced way.
* This artery module is only valid for 1D grid.
* This artery module has not been tested with multigrid.
* Implement a generic pressure ($a^m$ - $a^n$)
* Implement a viscous term $f_r$
*/

/**
# A solver for the 1D blood flow equations

One-dimensional (1D) blood flow models have been successfully used to
describe the flow of blood in large arteries of the systemic
circulation. They are valuable and efficient tools to capture pulse
wave propagation in large arterial networks and obtain satisfactory
evaluations of average quantities such as the flow rate $q$, the
cross-sectional area $a$ or the pressure $p$.

The 1D blood equations are:
$$
\left\{
\begin{aligned}
& \partial_t a + \partial_x q = 0 \\
& \partial_t q +\partial_x (q^2/a) = - a/\rho\partial_x p - f_r,
\end{aligned}
\right.
$$ 
where $a$ is the cross-sectional area of the artery, $q$ is the
flow rate, $p$ is the pressure, $\rho$ is the blood density and
$f_r$ is the viscous friction term, defined as:
$$
f_r = \phi \nu q/a,
$$
with $\nu$ the dynamic viscosity of blood. We typically set
$\phi = 8\pi$.

To close the 1D system of equations, we introduce the following
pressure law, describing the elastic behavior of the arterial wall:
$$
p = p_0 + k [\sqrt{a} - \sqrt{a_0}],
$$
where $a_0$ is the cross-sectional area at rest of the artery, $k$ is
the arterial wall rigidity and $p_0$ is the pressure applied on the
exterior of the artery.

## Variables and parameters

The primary fields are the arterial cross-sectional area $a$ and the
flow rate $q$. The secondary fields are the vessel properties and
correspond to the cross-sectional area at rest $a_0$ and the wall
rigidity $k$. The order of the declaration of these fields is
important as it conditions the order in which the variables will be
refined. */

scalar k[], a0[], a[], q[];

/**
The following parameters are used to define the zero and to
characterize algorithm convergence. Their values should not be
modified. */

double dry = 1e-30;
double epsilon = 1.e-14;
int maxit = 100;

/**
## Flux

As the flux can be different for arteries and veins depending on the
chosen pressure-area relationship, the user must define a function
which, given the left and right states of each field, provides the
solution to the homogeneous Riemann problem (without viscous or
viscoelastic terms). In case of varying geometrical and mechanical
properties, this flux function must involve a well-balancing
reconstruction procedure and a well-balanced approximation of the
topography source term. */

void (* riemann) (double km, double kp, double a0m, double a0p,
		  double am, double ap, double qm, double qp,
		  double Delta, double * fal, double * far,
		  double * fql, double * fqr, double * dtmax) = NULL;

/**
## Time-integration

The time-integration is performed using the generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
This generic time-integration scheme needs to know which fields are
updated. This list will be updated in the *defaults* event below. */

scalar * evolving;

/**
It also requires a function to compute the derivatives of the
conserved variables. Note that contrary to the
[saint-venant](saint-venant.h) solver, we do not need to overload the
*advance* function. */

double update_artery (scalar * evolving, scalar * updates, double dtmax);

/**
#### Initialisation and cleanup

We use the defaults event defined in
[predictor-corrector](predictor-corrector.h) to setup the initial
conditions. */

event defaults (i = 0) {

  evolving = list_concat ({a}, {q});

  foreach()
    for (scalar s in evolving)
      s[] = 0.;
  
  update = update_artery;
  
#if TREE
  for (scalar s in {k, a0, a, q}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
    }
  
  /**
  We also define special restriction and prolongation for *eta*. */
  
#endif
}

/**
At the end of the simulation we need to free any memory allocated in
the *defaults* event. */

event cleanup (i = end, last) {
  free (evolving);
}

/**
The event below will be executed after all the other initial events to
take into account user-defined field initialisations. */

event init (i = 0);

/**
#### Update function used in [predictor-corrector](predictor-corrector.h)
*/

double update_artery (scalar * evolving, scalar * updates, double dtmax) {
  
  /**
  We first recover the currently evolving cross-sectional area and
  flow rate. */

  scalar a = evolving[0], q = evolving[1];

  /**
  The variables *Fa* and *Fq* will contain the fluxes for $a$ and $q$
  respectively. The index $l$ and $r$ respectively correspond to the
  left and right of each mesh face. */

  face vector Fal[], Far[], Fql[], Fqr[];

  /**
  The gradients, used for high order spatial reconstruction of the
  conserved fields, are stored in locally-allocated
  fields. First-order reconstruction is used for the gradient
  fields. */

  vector gk[], ga0[], ga[], gq[];
  for (scalar s in {gk, ga0, ga, gq}) {
#if TREE
    s.prolongation = refine_linear;
#endif
  }
  gradients ({k, a0, a, q}, {gk, ga0, ga, gq});

  /**
  We are ready to compute the fluxes through each face of the
  domain. */

  foreach_face (reduction (min:dtmax)) {
    
    /**
    We use the central values of each scalar/vector quantity and the
    gradients to compute the left and right states of each field at
    each mesh face. The index $m$ and $p$ correspond respectively to
    the left and right of each mesh face. */

    double dx = Delta/2.;
    double kp = k[] - dx*gk.x[];
    double a0p = a0[] - dx*ga0.x[];
    double ap = a[] - dx*ga.x[];
    double qp = q[] - dx*gq.x[];
    double km = k[-1] + dx*gk.x[-1];
    double a0m = a0[-1] + dx*ga0.x[-1];
    double am = a[-1] + dx*ga.x[-1];
    double qm = q[-1] + dx*gq.x[-1];
    
    /**
    We then call the generic Riemann solver and store the computed
    fluxes in temporary fields. */
    
    double fal, far, fql, fqr;  
    riemann (km, kp, a0m, a0p, am, ap, qm, qp, Delta*cm[]/fm.x[],
	     &fal, &far, &fql, &fqr, &dtmax);
    foreach_dimension() {
      Fal.x[] = fm.x[]*fal;
      Far.x[] = fm.x[]*far;
      Fql.x[] = fm.x[]*fql;
      Fqr.x[] = fm.x[]*fqr;
    }
  }

  /**
  The update for each scalar quantity is the divergence of the
  fluxes. */

  scalar da = updates[0], dq = updates[1];
  foreach() {
    da[] = 0.;
    dq[] = 0.;
    foreach_dimension() {
      da[] += (Far.x[] - Fal.x[1])/(cm[]*Delta);
      dq[] += (Fqr.x[] - Fql.x[1])/(cm[]*Delta);
    }
  }
  
  return dtmax;
}

/**
## HLL fluxes

#### Help functions for the HLL GLU flux
*/

double get_asr(double alpha, double asl)
{
  return (1. - (1. - alpha)*asl)/alpha;
}

double phi1(double kl, double kr, double alpha, double asl)
{
  double asr = get_asr (alpha, asl);  
  return kr*sqrt(asr) - kl*sqrt(asl);
}

double dphi1(double kl, double kr, double alpha, double asl)
{
  double asr = get_asr (alpha, asl); 
  return 0.5*(kr/sqrt(asr)*(1. - 1./alpha) - kl/sqrt(asl));
}

double phi2(double kl, double kr, double alpha, double beta, double asl)
{
  double asr = get_asr (alpha, asl);
  if (fabs(beta) < dry ) return phi1 (kl, kr, alpha, asl);
  else return 0.5*beta*(1./(asr*asr) - 1./(asl*asl)) + phi1 (kl, kr, alpha, asl);
}

double dphi2(double kl, double kr, double alpha, double beta, double asl)
{
  double asr = get_asr (alpha, asl);  
  if (fabs(beta) < dry ) return dphi1 (kl, kr, alpha, asl);
  else return - beta*((1. - 1./alpha)*pow(asr, -3.) - pow(asl, -3.))
	 + dphi1 (kl, kr, alpha, asl);
}

/**
Newton-Raphson algorithm coupled to a bisection algorithm
*/

double newtonraphsonbisection(double (*phi) (double, double, double,
					     double, double),
			      double (*dphi) (double, double, double,
					      double, double),
			      double x1, double x2, double kl, double kr,
			      double alpha, double beta, double scm)
{
  double x = 0.5*(x1 + x2);
  double fx = phi(kl,kr,alpha,beta,x) - scm;
  double f1 = phi(kl,kr,alpha,beta,x1) - scm;
  double f2 = phi(kl,kr,alpha,beta,x2) - scm;

  if (fabs(f1) <= dry) return x1;
  else if (fabs(f2) <= dry) return x2;
  else if (fabs(fx) <= dry) return x;
  
  double dx = x2 - x1, dxo = dx;
  double dfx = dphi(kl,kr,alpha,beta,x);

  int it = 0;
  while (fabs(dx) > epsilon && it < maxit) {
    it ++;
    if ( f1 > 0. || f2 < 0. ) printf("Error: f1>0 or f2<0");    
    dfx = dphi(kl,kr,alpha,beta,x);
    dxo = dx;
    if (((x - x2)*dfx - fx)*((x - x1)*dfx -fx) > 0.
	|| fabs(2.*fx) >fabs(dxo*dfx)) {
      dx = 0.5*(x2 - x1);
      x = x1 + dx;
    }
    else {
      dx = fx / dfx;
      x = x - dx;
    }
    fx = phi(kl,kr,alpha,beta,x) - scm;
    if (fx < 0.)  {
      x1 = x;
      f1 = fx;
    }
    else if (fx > 0.) {
      x2 = x;
      f2 = fx;
    }
  }
  return x;
}

double solve_reduced_bernoulli(double kl, double kr,
			       double alpha, double dp0sa)
{  
  double x2 = 0., x1 = 1./(1.-alpha); 
  if (-dp0sa > kl/sqrt(1.-alpha)) return x1;
  else if (-dp0sa < -kr/sqrt(alpha)) return x2 ;
  else {
    if (kr - kl < dp0sa) x1 = 1.;
    else if (kr - kl > dp0sa) x2 = 1.;
    else return 1.;
    return newtonraphsonbisection (phi2, dphi2, x1, x2, kl, kr,
				   alpha, 0., dp0sa);
  }
}

double solve_full_bernoulli(double kl, double kr,
			    double alpha, double beta, double dp0sa)
{
  double x1 = 0., x2 = 0.;
  double acl = pow(2.*beta/kl, 0.4);
  double acr = pow(2.*beta/kr, 0.4);

  if (1. < alpha*acr + (1.-alpha)*acl) {
    x1 = max((1. - alpha*acr)/(1.-alpha), dry);
    x2 = min(1. / (1.-alpha) - dry, acl);
  }
  else {
    x2 = acl;
    x1 = (1. - alpha*acr)/(1.-alpha);
  }

  if (phi2(kl,kr,alpha,beta,x2) < dp0sa) return x2;
  else if (phi2(kl,kr,alpha,beta,x1) > dp0sa) return x1;
  else {
    if (kr - kl < dp0sa) x1 = 1.;
    else if (kr - kl > dp0sa) x2 = 1.;
    else return 1.;
    return newtonraphsonbisection (phi2, dphi2, x1, x2, kl, kr,
				   alpha, beta, dp0sa);
  }
}

/**
#### High-order variable reconstruction
*/

double order1 (double s0, double s1, double s2)
{  
  return generic_limiter(0., 0.);
}

/**
#### HLL fluxes for the conservative system of equations

The user can choose among three pre-defined well-balanced Riemann
solvers. These are based on the HLL flux for the homogeneous
conservative blood flow equations. */

void hll(double k, double al, double ar, double ql, double qr,
	 double Delta, double * fal, double * far,
	 double * fql, double * fqr, double * dtmax)
{
  double ul = ql/al;
  double ur = qr/ar;
  double cl = sqrt(0.5*k*sqrt(al));
  double cr = sqrt(0.5*k*sqrt(ar));
  double SL = min(ul - cl, ur - cr);
  double SR = max(ul + cl, ur + cr);

  double fa = 0., fq = 0.;
  if (0. <= SL) {
    fa = al*ul;
    fq = al*(ul*ul + k*sqrt(al)/3.);
  }
  else if (0. >= SR) {
    fa = ar*ur;
    fq = ar*(ur*ur + k*sqrt(ar)/3.);
  }
  else {
    double fal = al*ul;
    double fql = al*(ul*ul + k*sqrt(al)/3.);
    double far = ar*ur;
    double fqr = ar*(ur*ur + k*sqrt(ar)/3.);
    fa = (SR*fal - SL*far + SL*SR*(ar - al))/(SR - SL);
    fq = (SR*fql - SL*fqr + SL*SR*(qr - ql))/(SR - SL);
  }

  *fal = fa;
  *far = fa;
  *fql = fq;
  *fqr = fq;
  
  double a = max(fabs(SL), fabs(SR));
  if (a > dry) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

/**
Next we define the hll flux combined with the hydrostatic
reconstruction. This scheme preserves the following equilibrium
states:
$$
\left\{
\begin{aligned}
& q/a = \mathrm{cst} \\
& p = \mathrm{cst}.
\end{aligned}
\right.
$$ */

void hll_hr(double km, double kp, double a0m, double a0p,
	    double am, double ap, double qm, double qp,
	    double Delta, double * fal, double * far,
	    double * fql, double * fqr, double * dtmax)
{
  double zl = km*sqrt(a0m), zr = kp*sqrt(a0p), zlr = min(zl, zr);
  double klr = max(km, kp);
  double al = pow(max(0., zlr + km*sqrt(am) - zl)/klr, 2.);
  double ql = qm/am*al;
  double ar = pow(max(0., zlr + kp*sqrt(ap) - zr)/klr, 2.);
  double qr = qp/ap*ar;
  hll(klr, al, ar, ql, qr, Delta, fal, far, fql, fqr, dtmax);
  *fql += - (klr*al*sqrt(al)/3. - km*am*sqrt(am)/3.);
  *fqr += - (klr*ar*sqrt(ar)/3. - kp*ap*sqrt(ap)/3.);
}

/**
We then introduce the hll flux combined with the low-Shapiro
hydrostatic reconstruction. This scheme preserves the following
equilibrium states:
$$
\left\{
\begin{aligned}
& q = \mathrm{cst} \\
& p = \mathrm{cst}.
\end{aligned}
\right.
$$ */

void hll_hrls(double km, double kp, double a0m, double a0p,
	      double am, double ap, double qm, double qp,
	      double Delta, double * fal, double * far,
	      double * fql, double * fqr, double * dtmax)
{
  double zl = km*sqrt(a0m), zr = kp*sqrt(a0p), zlr = min(zl, zr);
  double klr = max(km, kp);
  double al = pow(max(0., zlr + km*sqrt(am) - zl)/klr, 2.);
  double ql = qm;
  double ar = pow(max(0., zlr + kp*sqrt(ap) - zr)/klr, 2.);
  double qr = qp;
  hll(klr, al, ar, ql, qr, Delta, fal, far, fql, fqr, dtmax);
  *fql += - (klr*al*sqrt(al)/3. - km*am*sqrt(am)/3.);
  *fqr += - (klr*ar*sqrt(ar)/3. - kp*ap*sqrt(ap)/3.);
}

/**
Finally, we introduce the GLU hll flux. This scheme preserves exactly
the following equilibrium states:
$$
\left\{
\begin{aligned}
& q = \mathrm{cst} \\
& p = \mathrm{cst},
\end{aligned}
\right.
$$
and provides a very accurate estimation of the following equilibrium
states:
$$
\left\{
\begin{aligned}
& q = \mathrm{cst} \\
& \frac{1}{2}\rho u^2 + p = \mathrm{cst}.
\end{aligned}
\right.
$$
*/

void hll_glu(double kl, double kr, double a0l, double a0r,
	     double al, double ar, double ql, double qr,
	     double Delta, double * fal, double * far,
	     double * fql, double * fqr, double * dtmax)
{
  double ul = ql/al, ur = qr/ar;
  double cl = sqrt(0.5*kl*sqrt(al));
  double cr = sqrt(0.5*kr*sqrt(ar));
  double SL = min(0., min(ul - cl, ur - cr));
  double SR = max(0., max(ul + cl, ur + cr));
  
  double ahll = ((SR*ar - SL*al) - (qr - ql))/(SR-SL);
  double qhll = ((SR*qr - SL*ql) - (ar*(ur*ur + kr*sqrt(ar)/3.)
				    - al*(ul*ul + kl*sqrt(al)/3.)))/(SR-SL);
  double dp = kr*sqrt(ar) - kl*sqrt(al);
  double dp0 = kr*sqrt(a0r) - kl*sqrt(a0l);
         
  //double as = 2.*al*ar/(al + ar);
  double as = (al + ar + sqrt(al*ar))/3.;
  double dxs = (kr*ar*sqrt(ar) - kl*al*sqrt(al))/3. - as*(dp - dp0);
  double qs = qhll + dxs/(SR-SL);
  double asl = 0., asr = 0.;
  double SLasl = 0., SRasr = 0.;
 
  if (fabs(ahll) <= dry) {
    qs = 0.;
    SLasl = 0.;
    SRasr = 0.;
  }
  else {
    if (fabs(SL) <= dry) {
      SLasl = 0.;
      SRasr = (SR-SL)*ahll;
    }
    else if (fabs(SR) <= dry) {
      SLasl = - (SR - SL)*ahll;
      SRasr = 0.;
    }
    else {
      double alpha = SR/(SR - SL);
      double beta = qs*qs/pow(ahll, 5./2.);
      double dp0sa = dp0/sqrt(ahll);
      if (fabs(qs) <= dry) {
	asl = solve_reduced_bernoulli (kl, kr, alpha, dp0sa);
      }
      else {  
        asl = solve_full_bernoulli (kl, kr, alpha, beta, dp0sa);
      }
      asr = get_asr (alpha, asl);
      SLasl = SL*asl*ahll;
      SRasr = SR*asr*ahll;
    }

    *fal = ql + (SLasl - SL*al);
    *fql = al*(ul*ul + kl*sqrt(al)/3.) + SL*(qs - ql);
    *far = qr + (SRasr - SR*ar);
    *fqr = ar*(ur*ur + kr*sqrt(ar)/3.) + SR*(qs - qr);

    double a = max(fabs(SL), fabs(SR));
    if (a > dry) {
      double dt = CFL*Delta/a;
      if (dt < *dtmax) *dtmax = dt;
    }
  }
}