/**
# Vertically-staggered non-hydrostatic solver

Staggered non-hydrostatic solver described in section 3.6.1 of [Popinet
(2020)](/Bibliography#popinet2020), updated to work with the current multilayer
framework ([hydro.h](hydro.h), [implicit.h](implicit.h)). It should be a drop-in replacement for
[nh.h](nh.h): just use

~~~c
#include "layered/hydro.h"
#include "layered/nh-staggered.h"  // instead of layered/nh.h
~~~

As in [nh.h](nh.h), the vertical velocity *w* is advected as a tracer (which
conserves the layer-integrated vertical momentum $hw$ through advection and
remapping), and the non-hydrostatic pressure *phi* is staggered vertically with
respect to *w*. The main differences with the (Keller box, collocated) [nh.h](nh.h)
is the linear system for $\phi_k$ in each column is tridiagonal (can be solved
with Thomas algorithm) rather than Hessenberg.

In contrast with the original 2019 version, the semi-implicit free-surface
coupling of [implicit.h](implicit.h) is retained: the coupled system for $\phi^{n+1}$ and
$\eta^{n+1}$ is solved with the multigrid solver, as in [nh.h](nh.h). Wave breaking
can be parameterised with the *breaking* parameter (off by default).

/!\ With few layer, it is better to use the keller scheme; this version is
intended to be used for many layers.

*/

#define NH 1
#include "implicit.h"

scalar w, phi;
mgstats mgp;
double breaking = HUGE;

event defaults (i = 0)
{
  hydrostatic = false;
  mgp.nrelax = 4;
  assert (nl > 0);
  w = new scalar[nl];
  phi = new scalar[nl];
  reset ({w, phi}, 0.);
  if (!linearised)
    tracers = list_append (tracers, w);
}

/**
## Viscous term

Vertical diffusion of the vertical velocity *w* (no-flux at the
surface, zero at the bottom). */

event viscous_term (i++)
{
  if (nu > 0.)
    foreach()
      vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.);
}

#define PHIT 6. // extrapolation coefficient of the top Dirichlet BC

/**
This applies the vertical pressure gradient correction to the
vertical velocity
$$
w^{n+1}_l = w^\star_l - \frac{\Delta t}{h_l}\,  [\phi]_l
$$
with the top Dirichlet condition $\phi_{surface} = 0$ imposed through
a quadratic extrapolation (PHIT) of $\phi$ below the surface. */


static void correct_w (double dt)
{
  foreach()
    foreach_layer() {
      double phit;
      if (point.l == nl - 1) {
	double phib = point.l > 0 ? phi[0,0,-1] : phi[];
	phit = (phib - PHIT*phi[])/3.; // top BC (Dirichlet) Eq. (23)
      }
      else
	phit = phi[0,0,1];
      if (h[] > dry)
	w[] -= dt*(phit - phi[])/h[];
  }
}

/**
## Face coefficients

The coefficients of the staggered horizontal pressure gradient at a
face.  `i` is the face offset (0: left face, 1: right face
of the cell). */

// todo? store coeffs1 for each mgsolve for efficiency?

static void coeffs1 (Point point, int i,
		     double * a, double * b0, double * b1, double * c)
{
  double zr = zb[i], zl = zb[i-1];
  foreach_layer() {
    int l = point.l;
    b0[l] = - (zr - zl - 2.*h[i-1])/(2.*Delta);
    b1[l] = - (zr - zl + 2.*h[i])/(2.*Delta);
    if (l == nl - 1) {
      a[l] = (zr + h[i] - zl - h[i-1])/3./(2.*Delta);
      b0[l] -= PHIT*a[l];
      b1[l] -= PHIT*a[l];
      c[l] = 0.;
    }
    else {
      a[l] = 0.;
      c[l] = (zr + h[i] - zl - h[i-1])/(2.*Delta);
    }
    zr += h[i], zl += h[i-1];
  }
}

/**
The staggered face "acceleration" $a_x = -[\nabla(h\phi) -
\phi\nabla z]/h$ for each layer, at face *i*. */


static void face_gradient (Point point, scalar phi, int i, double * ax)
{
  double ac[nl], bc0[nl], bc1[nl], cc[nl];
  coeffs1 (point, i, ac, bc0, bc1, cc);
  foreach_layer() {
    int l = point.l;
    double phibL = l > 0      ? phi[i-1,0,-1] : phi[i-1];
    double phibR = l > 0      ? phi[i,0,-1]   : phi[i];
    double phitL = l < nl - 1 ? phi[i-1,0,1]  : phi[i-1];
    double phitR = l < nl - 1 ? phi[i,0,1]    : phi[i];
    double hl = h[i-1], hr = h[i];
    ax[l] = (hl + hr > dry) ?
      - fm.x[i]*(ac[l]*(phibL + phibR) +
		 bc0[l]*phi[i-1] + bc1[l]*phi[i] +
		 cc[l]*(phitL + phitR))/((hl + hr)/2.) : 0.;
  }
}

/**
This macro computes the height-weighted face value
$pg = -h_f a_x$ for each layer at face *i* (the quantity added to the
face acceleration *ha* and fluxes *hu*), and executes *code* for each
layer. It plays the role of the `hpg` macro of [hydro.h](hydro.h) for the
staggered discretisation. */


#define shpg(pg, phi, i, code) do {					\
  double _ax[nl];							\
  face_gradient (point, phi, i, _ax);					\
  double pg = 0.;							\
  foreach_layer() {                                                     \
    pg = - hf.x[i]*_ax[point.l];                                        \
    code;								\
  }									\
} while (0)


/**
## Tridiagonal system

For the staggered scheme, the linear system for the column of
$\phi_l$ is tridiagonal. The row for layer $l$ is
$$
a_l \phi_{l-1} + b_l \phi_l + c_l \phi_{l+1} = d_l
$$

with the vertical $[\phi]_l$ terms, the horizontal terms and, the semi-implicit
free-surface coupling $+\theta\nabla\cdot(h_f\nabla\eta)$. The top Dirichlet
condition is included in the coefficients of the last row.

The horizontal terms do not include the metric factors $cm/fm$ (valid for
Cartesian metrics). */


static void matrix (Point point, scalar phi, scalar eta, scalar rhs,
		    double * a, double * b, double * c, double * d)
{
  double al[nl], bl0[nl], bl1[nl], cl[nl];
  double ar[nl], br0[nl], br1[nl], cr[nl];
  coeffs1 (point, 0, al, bl0, bl1, cl); // left face
  coeffs1 (point, 1, ar, br0, br1, cr); // right face
  foreach_layer () {

    int l = point.l;
    d[l] = rhs[];
    c[l] = 1.0;
    if (l == 0)
      a[l] = 0.; // bottom BC (Neumann)
    else
      a[l] = h[]/h[0,0,-1];
    b[l] = - a[l] - c[l];
    if (l == nl - 1) {
      // top BC (Dirichlet, phi_surface = 0)
      if (l == 0)
	b[l] -= (PHIT - 1.)*c[l]/3.;
      else {
	b[l] -= PHIT*c[l]/3.;
	a[l] += c[l]/3.;
      }
    }

    // horizontal non-hydrostatic coupling
    a[l] -= theta_H*h[]*(ar[l] - al[l])/Delta;
    b[l] -= theta_H*h[]*(br0[l] - bl1[l])/Delta;
    c[l] -= theta_H*h[]*(cr[l] - cl[l])/Delta;
    double phib1  = l > 0      ? phi[1,0,-1]  : phi[1];
    double phib_1 = l > 0      ? phi[-1,0,-1] : phi[-1];
    double phit1  = l < nl - 1 ? phi[1,0,1]   : phi[1];
    double phit_1 = l < nl - 1 ? phi[-1,0,1]  : phi[-1];
    d[l] += theta_H*h[]*
      ((ar[l]*phib1 + br1[l]*phi[1] + cr[l]*phit1) -
       (al[l]*phib_1 + bl0[l]*phi[-1] + cl[l]*phit_1))/Delta;

/**
## Semi-implicit coupling coefficient:

The implicit barotropic gradient enters the stored (theta-weighted) flux
twice: once through the theta-splitting of the gravity term
(a_baro(eta^{n+1}) = theta_H*g*grad(eta^{n+1})), and once through the
theta-weighting of the flux itself (hu = theta_H*(hf*uf + dt*ha), where
the "- dt*ha" of the implicit rebuild is cancelled by hydro.h's
"hu += dt*ha"). The eta-dependent part of the corrected flux is thus

    theta_H^2*dt*hf*a_baro(eta^{n+1}),

and the phi-row must anticipate exactly this divergence for the solved
(phi,eta) to leave the corrected state divergence-free.
*/
    /* d[l] += theta_H*(1-theta_H)*h[]*(hf.x[1]*a_baro (eta, 1) - */
    /*     		 hf.x[]*a_baro (eta, 0))/(Delta*cm[]); */
    d[l] += sq(theta_H)*h[]*(hf.x[1]*a_baro (eta, 1) -
        		 hf.x[]*a_baro (eta, 0))/(Delta*cm[]);
  }
}

/**
## Relaxation operator

The tridiagonal system of each column is solved exactly with the
Thomas algorithm, and the value of $\eta$ is then updated exactly
given the new $\phi$ (as in [nh.h](nh.h)). */

trace
static void relax_nh (scalar * phil, scalar * rhsl, int lev, void * data)
{
  scalar phi = phil[0], rhs = rhsl[0];
  scalar eta = phil[1], rhs_eta = rhsl[1];
  face vector alpha = *((vector *)data);
#if GAUSS_SEIDEL || _GPU
  for (int parity = 0; parity < 2; parity++)
    foreach_level_or_leaf (lev)
      if (level == 0 || ((point.i + parity) % 2) != (point.j % 2))
#else
  foreach_level_or_leaf (lev)
#endif
  {
    double ta[nl], tb[nl], tc[nl], td[nl];
    matrix (point, phi, eta, rhs, ta, tb, tc, td);
    for (int l = 1; l < nl; l++) {
      tb[l] -= ta[l]*tc[l-1]/tb[l-1];
      td[l] -= ta[l]*td[l-1]/tb[l-1];
    }
    ta[nl-1] = td[nl-1]/tb[nl-1];
    for (int l = nl - 2; l >= 0; l--)
      ta[l] = (td[l] - tc[l]*ta[l+1])/tb[l];
    foreach_layer () {
      phi[] = ta[point.l];
    }

    /**
    The value of $\eta$ is updated using the same discretisation of
    the pressure gradient as the residual and the correction below. */

    double n = 0.;
    foreach_dimension() {
      shpg (pg, phi, 0, n -= pg);
      shpg (pg, phi, 1, n += pg);
    }
    n *= theta_H*sq(dt);
    double d = rigid ? 0. : - cm[]*Delta;
    n -= cm[]*Delta*rhs_eta[];
    eta[] = 0.;
    foreach_dimension() {
      n += alpha.x[0]*a_baro (eta, 0) - alpha.x[1]*a_baro (eta, 1);
      diagonalize (eta)
	d -= alpha.x[0]*a_baro (eta, 0) - alpha.x[1]*a_baro (eta, 1);
    }
    eta[] = n/d;
  }
}

/**
## Residual computation

The residual for $\phi$ is computed from the same tridiagonal system
as the relaxation ("naive" discretisation, only 1st order on
trees). The residual for $\eta$ uses the face field $g$ combining the
staggered non-hydrostatic pressure gradient and the implicit
barotropic term, as in [nh.h](nh.h). */

trace
static double residual_nh (scalar * phil, scalar * rhsl,
			   scalar * resl, void * data)
{
  scalar phi = phil[0], rhs = rhsl[0], res = resl[0];
  scalar eta = phil[1], rhs_eta = rhsl[1], res_eta = resl[1];
  face vector alpha = *((vector *)data);
  double maxres = 0.;
  face vector g = new face vector[nl];
  foreach_face() {
    double pgh = theta_H*a_baro (eta, 0);
    shpg (pg, phi, 0,
	  g.x[] = - 2.*(pg + hf.x[]*pgh));
  }
  foreach (reduction(max:maxres)) {
    double a[nl], b[nl], c[nl], d[nl];
    matrix (point, phi, eta, rhs, a, b, c, d);
    foreach_layer() {
      int l = point.l;
      res[] = d[l] - b[l]*phi[];
      if (l > 0)
	res[] -= a[l]*phi[0,0,-1];
      if (l < nl - 1)
	res[] -= c[l]*phi[0,0,1];
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
    }
    res_eta[] = rhs_eta[] - (rigid ? 0. : eta[]);
    foreach_layer()
      foreach_dimension()
	res_eta[] += theta_H*sq(dt)/2.*(g.x[1] - g.x[])/(Delta*cm[]);
  }
  delete ((scalar *){g});
  return maxres;
}

/**
## Coupled solution

The rhs is the divergence form of the predicted fluxes and
vertical momentum, scaled by $1/\Delta t$
*/

event pressure (i++)
{
  scalar rhs = new scalar[nl];
  double h1 = 0., v1 = 0.;
  foreach (reduction(+:h1) reduction(+:v1)) {
    double zl = zb[-1], zr = zb[1];
    foreach_layer() {
      double ufl = hf.x[]  > dry ? hu.x[]/hf.x[]    : 0.;
      double ufr = hf.x[1] > dry ? hu.x[1]/hf.x[1]  : 0.;
      rhs[] = h[]*w[] - h[]*(ufl + ufr)*(zr + h[1] - zl - h[-1])/(4.*Delta);
      if (point.l > 0) {
	double ul = hf.x[0,0,-1] > dry ? hu.x[0,0,-1]/hf.x[0,0,-1] : 0.;
	double ur = hf.x[1,0,-1] > dry ? hu.x[1,0,-1]/hf.x[1,0,-1] : 0.;
	rhs[] -= h[]*w[0,0,-1];
	rhs[] += h[]*(ul + ur)*(zr - zl)/(4.*Delta);
      }
      foreach_dimension()
	rhs[] += h[]*(hu.x[1] - hu.x[])/Delta;
      rhs[] /= dt;
      h1 += dv()*h[];
      v1 += dv();
      zl += h[-1], zr += h[1];
    }
  }

  /**
  The fields used by the relaxation function above need to be
  restricted to all levels (the other fields `cm`, `fm`,
  `alpha_eta` were already restricted in
  [implicit.h](implicit.h)). */

  restriction ({zb, h, hf});

  /**
  We then call the multigrid solver to get both the non-hydrostatic
  pressure $\phi$ and the free-surface elevation $\eta$. */


  scalar res;
  if (res_eta.i >= 0)
    res = new scalar[nl];
  mgp = mg_solve ({phi,eta_r}, {rhs,rhs_eta}, residual_nh, relax_nh, &alpha_eta,
		  res = res_eta.i >= 0 ? (scalar *){res,res_eta} : NULL,
		  nrelax = 4, minlevel = 1,
		  tolerance = TOLERANCE*sq(h1/(dt*v1)));
  delete ({rhs});
  if (res_eta.i >= 0)
    delete ({res});


  /**
  The non-hydrostatic pressure gradient is added to the face-weighted
  acceleration *ha* and to the face fluxes *hu* (which store
  $\theta_H(hu)^{n+1}$, hence the factor $\theta_H\Delta t$). */

  face vector su[];
  foreach_face() {
    su.x[] = 0.;
    shpg (pg, phi, 0,
	  ha.x[] += pg;
	  su.x[] -= pg;
	  hu.x[] += theta_H*dt*pg );
  }


  /**
  The vertical pressure gradient is applied to the vertical velocity
  *w*, which is then limited by the breaking parameter in the same
  way as in [nh.h](nh.h). */

  correct_w (dt);
  foreach() {
    double wmax = HUGE;
    if (breaking < HUGE) {
      wmax = 0.;
      foreach_layer()
	wmax += h[];
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }

    foreach_layer() {
      if (h[] > dry && fabs (w[]) > wmax)
	w[] = (w[] > 0. ? 1. : -1.)*wmax;
    }

    /**
    The rhs for $\eta^{n+1}$ is updated for the second pass of the
    [semi-implicit free-surface solver](implicit.h). */

    foreach_dimension()
      rhs_eta[] += theta_H*sq(dt)*(su.x[1] - su.x[])/(Delta*cm[]);
  }
}

/**
## Cleanup */

event cleanup (i = end, last) {
  delete ({w, phi});
}
