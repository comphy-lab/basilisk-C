/**
# A hydrostatic reconstruction solver for the 1D blood flow equations

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
rigidity $k$. We also introduce the terciary fields $zb=k\sqrt{a_0}$,
$h=k\sqrt{a}$ and $eta=h-zb$, which are equilibrium related
quantities. The order of the declaration of these fields is important
as it conditions the order in which the variables will be refined (*k*
before *zb*, *zb* before *a* and *a* before *eta*). */

scalar k[], zb[], a[], eta[];
scalar q[];

/**
The following parameters are used to define the zero and to
characterize algorithm convergence. Their values should not be
modified. */

double dry = 1.e-30;

/**
The *SEPS* constant is used to avoid division by zero. */

#undef SEPS
#define SEPS 1.e-30

/**
## Time-integration

The time-integration is performed using the generic second-order
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
#### Setup

This generic time-integration scheme needs to know which fields are
updated. This list will be updated in the *defaults* event below. */

scalar * evolving = NULL;

/**
We need to overload the default advance function of the
predictor-corrector scheme. */

trace
static void advance_bloodflow (scalar * output, scalar * input, 
			       scalar * updates, double dt)
{
  // recover scalar and vector fields from lists
  scalar ai = input[0], ao = output[0], da = updates[0];
  scalar qi = input[1], qo = output[1], dq = updates[1];

  foreach() {
    ao[] = ai[] + dt*da[];
    qo[] = qi[] + dt*dq[];

    eta[] = k[]*sqrt (ao[]) - zb[];
  }
    
  // fixme: on trees eta is defined as eta = h - zb and not ho -
  // zb in the refine_eta() and restriction_eta() functions below
  boundary ((scalar *) {ao, qo, eta});
}

/**
When using an adaptive discretisation (i.e. a tree)., we need to make
sure that *eta* is maintained as *h-zb* whenever cells are refined or
restricted. */

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = k[]*sqrt (a[]) - zb[];
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = k[]*sqrt (a[]) - zb[];
}
#endif

/**
#### HLL flux
*/

void hll(double k, double al, double ar, double ql, double qr,
	 double Delta, double * fa, double * fq, double * dtmax)
{
  double ul = ql/al;
  double ur = qr/ar;
  double cl = sqrt (0.5*k*sqrt (al));
  double cr = sqrt (0.5*k*sqrt (ar));
  double SL = min (ul - cl, ur - cr);
  double SR = max (ul + cl, ur + cr);

  if (0. <= SL) {
    *fa = al*ul;
    *fq = al*(ul*ul + k*sqrt (al)/3.);
  }
  else if (0. >= SR) {
    *fa = ar*ur;
    *fq = ar*(ur*ur + k*sqrt (ar)/3.);
  }
  else {
    double fal = al*ul;
    double fql = al*(ul*ul + k*sqrt (al)/3.);
    double far = ar*ur;
    double fqr = ar*(ur*ur + k*sqrt (ar)/3.);
    *fa = (SR*fal - SL*far + SL*SR*(ar - al))/(SR - SL);
    *fq = (SR*fql - SL*fqr + SL*SR*(qr - ql))/(SR - SL);
  }

  double a = max (fabs (SL), fabs (SR));
  if (a > dry) {
    double dt = CFL*Delta/a;
    if (dt < *dtmax)
      *dtmax = dt;
  }
}

/**
#### Update function used in [predictor-corrector](predictor-corrector.h)
*/

trace
double update_bloodflow (scalar * evolving, scalar * updates, double dtmax) {
  
  /**
  We first recover the currently evolving cross-sectional area and
  flow rate. */

  scalar q = evolving[1];

  /**
  The variables *Fa* and *Fq* will contain the fluxes for $a$ and $q$
  respectively. The field *S* is necessary to store the topography
  source term. */

  face vector Fa[], Fq[], S[];

  /**
  The gradients, used for high-order spatial reconstruction of the
  conserved fields, are stored in locally-allocated
  fields. First-order reconstruction is used for the gradient
  fields. Note that to preserve the equilibrium states, the spatial
  reconstruction accounts for the equilibrium states the method is
  designed to preserve.  */

vector gk[], gzb[], geta[], gq[];
  for (scalar s in {gk, gzb, geta, gq}) {
    s.gradient = zero;
#if TREE
    s.prolongation = refine_linear;
#endif
  }
  gradients ({k, zb, eta, q}, {gk, gzb, geta, gq});

  /**
  ###### Left/right state reconstruction
  */

  /**
  We are ready to compute the fluxes through each face of the
  domain. */

  foreach_face (reduction (min:dtmax)) {
    
    /**
    We use the central values of each scalar/vector quantity and the
    gradients to compute the left and right states of each field at
    each cell face. The index $l$ and $r$ correspond respectively to
    the left and right of each cell face. */

    double dx = Delta/2.;

    double kr = k[] - dx*gk.x[];
    double kl = k[-1] + dx*gk.x[-1];

    double zbr = zb[] - dx*gzb.x[];
    double zbl = zb[-1] + dx*gzb.x[-1];

    double etar = eta[] - dx*geta.x[];
    double etal = eta[-1] + dx*geta.x[-1];

    double qr = q[] - dx*gq.x[];
    double ql = q[-1] + dx*gq.x[-1];

    /**
    We then compute the left and right states of *a*. */

    double ar = sq ((zbr + etar)/kr);
    double al = sq ((zbl + etal)/kl);
    
    /**
    As part of the hydrostatic reconstruction procedure, we define
    left and right states (m,p) at each cell face of the conserved
    variables *a* and *q*. */
    
    double klr = max (kl, kr);
    double zblr = min (zbl, zbr);
    
    double ap = sq (max (0., zblr + etar)/klr);
    double am = sq (max (0., zblr + etal)/klr);

    double qp = qr; // qr/ar*ap
    double qm = ql; // ql/al*am
   
    /**
    We then call the generic HLL Riemann solver and store the computed
    fluxes in temporary fields. */
    
    double fa, fq;  
    hll (klr, am, ap, qm, qp, Delta*cm[]/fm.x[], &fa, &fq, &dtmax);

    /**
    ###### Topography source term
    */
    
    double sr = klr*ap*sqrt (ap)/3. - kr*ar*sqrt (ar)/3.;
    double sl = klr*am*sqrt (am)/3. - kl*al*sqrt (al)/3.;
    
    foreach_dimension() {
      Fa.x[] = fm.x[]*fa;
      Fq.x[] = fm.x[]*(fq - sr); // Used to update the face 0 of the cell
      S.x[]  = fm.x[]*(fq - sl); // Used to update the face 1 of the cell
    }
  }

  /**
  ###### Updates for evolving quantities

  The update for each scalar quantity is the divergence of the
  fluxes. */

  scalar da = updates[0], dq = updates[1];
  foreach() {
    da[] = 0.;
    dq[] = 0.;
    foreach_dimension() {
      da[] += (Fa.x[] - Fa.x[1])/(cm[]*Delta);
      dq[] += (Fq.x[] - S.x[1,0])/(cm[]*Delta);
    }

    /**
    An extra cell centered term must be added to the momentum update
    to preserve consistency when second-order reconstruction is
    used. This term is 0 for first-order reconstruction. In the
    following, to match the indices above, the indices (l,r) now refer
    respectively to the left of the right face and to the right of the
    left face of each cell. */

    double dx = Delta/2.;

    double kr = k[] - dx*gk.x[];
    double kl = k[] + dx*gk.x[];

    double zbr = zb[] - dx*gzb.x[];
    double zbl = zb[] + dx*gzb.x[];

    double etar = eta[] - dx*geta.x[];
    double etal = eta[] + dx*geta.x[];

    double ar = sq ((zbr + etar)/kr);
    double al = sq ((zbl + etal)/kl);
    
    double kc = sqrt (kl*kr);
    double zbc = 0.5*(zbl + zbr);
    double acr = sq ((zbc + etar)/(kc + SEPS));
    double acl = sq ((zbc + etal)/(kc + SEPS));

    dq[] -= (
    	     (kc*acl*sqrt(acl)/3. - kl*al*sqrt (al)/3.) -
    	     (kc*acr*sqrt(acr)/3. - kr*ar*sqrt (ar)/3.)
    	     )/(cm[]*Delta);
  }
  return dtmax;
}

/**
#### Initialisation and cleanup

We use the defaults event defined in
[predictor-corrector](predictor-corrector.h) to setup the initial
conditions. */

event defaults (i = 0)
{
  evolving = list_concat ({a}, {q});

  foreach()
    for (scalar s in evolving)
      s[] = 0.;

  advance = advance_bloodflow;
  update = update_bloodflow;
  
#if TREE
  for (scalar s in {k, zb, a, eta, q}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
    }
  
  /**
  We also define special reconstructions for *eta*. */

  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
#endif
}

/**
The event below will be executed after all the other initial events to
take into account user-defined field initialisations. */

event init (i = 0)
{
  foreach()
    eta[] = k[]*sqrt (a[]) - zb[];
}

/**
At the end of the simulation we need to free any memory allocated in
the *defaults* event. */

event cleanup (i = end, last)
{
  free (evolving);
}

/**
#### Well-balanced reconstruction of $a$ */

#include "elevation.h"
