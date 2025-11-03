/**
# Non-hydrostatic extension of the multilayer solver

This adds the non-hydrostatic terms of the [vertically-Lagrangian
multilayer solver for free-surface flows](hydro.h) described in
[Popinet, 2020](/Bibliography#popinet2020). The corresponding system
of equations is
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta) 
  {\color{blue} - \mathbf{{\nabla}} (h \phi)_k + \left[ \phi 
      \mathbf{{\nabla}} z \right]_k},\\
  {\color{blue} \partial_t (hw)_k + \mathbf{{\nabla}} \cdot 
  \left( hw \mathbf{u} \right)_k} & {\color{blue} = - [\phi]_k,}\\
  {\color{blue} \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k + 
  \left[ w - \mathbf{u} \cdot \mathbf{{\nabla}} z \right]_k} & 
  {\color{blue} = 0},
\end{aligned}
$$
where the terms in blue are non-hydrostatic.

The additional $w_k$ and $\phi_k$ fields are defined. The convergence
statistics of the multigrid solver are stored in *mgp*.

Wave breaking is parameterised usng the *breaking* parameter, which is
turned off by default (see Section 3.6.4 in [Popinet,
2020](/Bibliography#popinet2020)). 

Note that this version differs in many ways from that presented in
[Popinet, 2020](/Bibliography#popinet2020). The most important
differences are the [time-implicit integration](implicit.h) of the
barotropic free-surface evolution (explicit in the previous version)
and the exact projection of the baroclinic velocities (approximate in
the previous version). This results in the solution of a coupled
system for the non-hydrostatic pressure $\phi^{n+1}$ and the
free-surface elevation $\eta^{n+1}$. 

See also the [general introduction](README). */

#define NH 1
#include "implicit.h"

scalar w, phi;
mgstats mgp;
double breaking = HUGE;

/**
## Setup

The $w_k$ and $\phi_k$ scalar fields are allocated and the $w_k$ are
added to the list of advected tracers. */

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

Vertical diffusion is added to the vertical component of velocity
$w$. */

event viscous_term (i++)
{
  if (nu > 0.)
    foreach()
      vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.);
}

/**
## Assembly of the Hessenberg matrix

For the Keller box scheme, the linear system of equations verified by
the non-hydrostatic pressure $\phi$ is expressed as an [Hessenberg
matrix](https://en.wikipedia.org/wiki/Hessenberg_matrix) for each column.

The Hessenberg matrix $\mathbf{H}$ for a column at a particular *point* is
stored in a one-dimensional array with `nl*nl` elements. It encodes
the coefficients of the left-hand-side of the Poisson equation as
$$
\begin{aligned}
  (\mathbf{H}\mathbf{\phi} - \mathbf{d})_l & =
  - \text{rhs}_l +
  h_l \nabla\cdot g_l^{n + \theta} +\\
  & 4 (\phi_{l + 1 / 2} - \phi_{l - 1 / 2}) + 8 h_l \sum^{l - 1}_{k = 0}
  (- 1)^{l + k}  \frac{\phi_{k + 1 / 2} - \phi_{k - 1 / 2}}{h_k}\\
  g_l^{n + \theta} & = \nabla (h^{\star}_k \phi^{n + \theta}_{k - 1 / 2} + 
  h^{\star}_k \phi^{n + \theta}_{k + 1 / 2}) 
  - 2 [\phi \nabla \hat{z}]^{n + \theta}_k
  + 2 \theta gh^{\star}_k \nabla \eta^{n + 1}
\end{aligned}
$$
where $\mathbf{\phi}$ is the vector of $\phi_l$ for this column and
$\mathbf{d}$ is a vector dependent only on the values of $\phi$ in the
neighboring columns. Note that in contrast with [Popinet,
2020](/Bibliography#popinet2020), all the (metric) terms are
retained. */

static void box_matrix (Point point, scalar phi, scalar rhs,
			face vector hf, scalar eta,
			double H[nl*nl], double d[nl])
{
  coord dz, dzp;
  foreach_dimension()
    dz.x = zb[] - zb[-1], dzp.x = zb[1] - zb[];
  foreach_layer()
    foreach_dimension()
      dz.x += h[] - h[-1], dzp.x += h[1] - h[];
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = h[0,0,m]/(sq(Delta)*cm[]);
    d[l] = rhs[0,0,m];
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    foreach_dimension() {
      double s = Delta*slope_limited((dz.x - h[0,0,m] + h[-1,0,m])/Delta);
      double sp = Delta*slope_limited((dzp.x - h[1,0,m] + h[0,0,m])/Delta);
      d[l] -= a*(gmetric(0)*(h[-1,0,m] - s)*phi[-1,0,m] +
		 gmetric(1)*(h[1,0,m] + sp)*phi[1,0,m] +
		 2.*theta_H*Delta*(hf.x[0,0,m]*a_baro (eta, 0) -
				   hf.x[1,0,m]*a_baro (eta, 1)));
      H[l*nl + l] -= a*(gmetric(0)*(h[0,0,m] + s) +
			gmetric(1)*(h[0,0,m] - sp));
    }
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      foreach_dimension() {
        double s = Delta*slope_limited(dz.x/Delta);
        double sp = Delta*slope_limited(dzp.x/Delta);
	d[l] -= a*(gmetric(0)*(h[-1,0,m] + s)*phi[-1,0,m+1] +
		   gmetric(1)*(h[1,0,m] - sp)*phi[1,0,m+1]);
	H[l*(nl + 1) - 1] -= a*(gmetric(0)*(h[0,0,m] - s) +
				gmetric(1)*(h[0,0,m] + sp));
      }
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = h[0,0,nl-1-k];
      if (hk > dry) {
	H[l*nl + k] -= 8.*s*h[0,0,m]/hk;
	H[l*nl + k - 1] += 8.*s*h[0,0,m]/hk;
      }
    }
    foreach_dimension()
      dz.x -= h[0,0,m] - h[-1,0,m], dzp.x -= h[1,0,m] - h[0,0,m];
  }
}

/**
## Relaxation operator */

#include "hessenberg.h"

face vector hf;

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

    /**
    The updated values of $\phi$ in a column are obtained as
    $$
    \mathbf{\phi} = \mathbf{H}^{-1}\mathbf{b}
    $$
    were $\mathbf{H}$ and $\mathbf{b}$ are the Hessenberg matrix and
    vector constructed by the function above. */

    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, hf, eta, H, b);
    solve_hessenberg (H, b);
    int l = nl - 1;
    foreach_layer()
      phi[] = b[l--];


    /**
    The value of $\eta$ also needs to be updated since it is solved
    implicitly and depends on $\phi$. */
    
    double n = 0.;
    foreach_dimension() {
      hpg (pg, phi, 0,
	   n -= pg);
      hpg (pg, phi, 1,
	   n += pg);
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
## Residual computation */

trace
static double residual_nh (scalar * phil, scalar * rhsl,
			   scalar * resl, void * data)
{
  scalar phi = phil[0], rhs = rhsl[0], res = resl[0];
  scalar eta = phil[1], rhs_eta = rhsl[1], res_eta = resl[1];
  double maxres = 0.;

  face vector g = new face vector[nl];
  foreach_face() {
    double pgh = theta_H*a_baro (eta, 0);
    hpg (pg, phi, 0,
	 g.x[] = - 2.*(pg + hf.x[]*pgh));
  }

  foreach (reduction(max:maxres)) {

    /**
    The residual for $\phi$ is computed as
    $$
    \begin{aligned}
    \text{res}_l = & \text{rhs}_l -
    h_l \nabla\cdot g_l^{n + \theta} -
    4 (\phi_{l + 1 / 2} - \phi_{l - 1 / 2}) - 8 h_l \sum^{l - 1}_{k = 0}
    (- 1)^{l + k}  \frac{\phi_{k + 1 / 2} - \phi_{k - 1 / 2}}{h_k}
    \end{aligned}
    $$
    */

    coord dz;
    foreach_dimension()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer() {
      res[] = rhs[] + 4.*phi[];      
      foreach_dimension() {
	res[] -= h[]*(g.x[1] - g.x[])/(Delta*cm[]);
	res[] += h[]*(g.x[] + g.x[1])/(hf.x[] + hf.x[1] + dry)*
	  slope_limited((dz.x + hf.x[1] - hf.x[])/(Delta*cm[]));
	if (point.l > 0)
	  res[] -= h[]*(g.x[0,0,-1] + g.x[1,0,-1])/
	    (hf.x[0,0,-1] + hf.x[1,0,-1] + dry)*slope_limited(dz.x/(Delta*cm[]));
      }
      if (point.l < nl - 1)
        res[] -= 4.*phi[0,0,1];
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
	double hk = h[0,0,k];
	if (hk > dry)
	  res[] += 8.*s*(phi[0,0,k] - phi[0,0,k+1])*h[]/hk;
      }
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      foreach_dimension()
	dz.x += hf.x[1] - hf.x[];
    }

    /**
    The residual for $\eta$ is computed as
    $$
    \text{res}_\eta = \text{rhs}_\eta - \beta\eta + 
    \theta_H \Delta t^2 \sum_l \nabla \cdot (- \nabla_z\phi_l + 
    \theta_H h_l g \nabla \eta)
    $$
    where $\beta = 1$ for a free surface and $\beta = 0$ for a rigid lid. */

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

The coupled system for $\phi_k^{n+1}$ and $\eta^{n+1}$ is solved using
the multigrid solver. */

event pressure (i++)
{

  /**
  The r.h.s. is computed as
  $$
  \frac{2 h_k}{\theta \Delta t}  \left( 2 w^n_k + 
  \nabla \cdot (hu)_k^{\star} - [u \cdot \nabla \hat{z}]^\star_k + 
  4 \sum^{k - 1}_{l = 0} (- 1)^{k + l} w^n_l \right)
  $$
  Note that the discrete approximation below must verify [Galilean
  invariance](https://en.wikipedia.org/wiki/Galilean_invariance) i.e.
  $$
  \begin{aligned}
  \nabla \cdot (h(u + U_0))_k^{\star} - [(u + U_0)\cdot 
  \nabla \hat{z}]^\star_k & = \nabla \cdot (hu)_k^{\star} - [u \cdot 
  \nabla \hat{z}]^\star_k
  \end{aligned}
  $$
  with $U_0$ an arbitrary constant vector. This can be simplified as
  $$
  \nabla \cdot (h_k^{\star}U_0) - [U_0\cdot \nabla \hat{z}]^\star_k = 0
  $$
  or further
  $$
  \nabla h_k^\star = [\nabla \hat{z}]^\star_k
  $$
  which is obviously verified (analytically) since by definition
  $$
  h_k^\star = [\hat{z}]^\star_k
  $$
  Note that it is not so obvious that this is verified numerically, as
  this depends on the choices made for several approximations. In
  particular, in the expressions below, Galilean invariance implies
  the relations

~~~c
  hu.x[] == U0*hf.x[];
  hu.x[1] == U0*hf.x[1];  
~~~

  which depend on the [detail of the calculation](hydro.h#face_fields)
  of `hu`. Note also that the slope limiter will break Galilean
  invariance. */

  scalar rhs = new scalar[nl];
  double h1 = 0., v1 = 0.;
  foreach (reduction(+:h1) reduction(+:v1)) {
    coord dz;
    foreach_dimension()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer() {
      rhs[] = 2.*w[];
      foreach_dimension()
        rhs[] += (hu.x[1] - hu.x[])/(Delta*cm[]) -
	  u.x[]*slope_limited((dz.x + hf.x[1] - hf.x[])/(Delta*cm[]));
      if (point.l > 0)
	foreach_dimension()
	  rhs[] += u.x[0,0,-1]*slope_limited(dz.x/(Delta*cm[]));
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
	rhs[] += 4.*s*w[0,0,k];
      rhs[] *= 2.*h[]/(theta_H*dt);
      foreach_dimension()
	dz.x += hf.x[1] - hf.x[];
      h1 += dv()*h[];
      v1 += dv();
    }
  }

  /**
  We then call the multigrid solver, using the relaxation and residual
  functions defined above, to get both the non-hydrostatic pressure
  $\phi$ and free-surface elevation $\eta$. */

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
  acceleration and to the face fluxes. */

  face vector su[];
  foreach_face() {
    su.x[] = 0.;
    hpg (pg, phi, 0,
	 ha.x[] += pg;
	 su.x[] -= pg;
	 hu.x[] += theta_H*dt*pg );
  }

  /**
  The maximum allowed vertical velocity is computed as
  $$
  w_\text{max} = b \sqrt{g | H |_{\infty}}
  $$
  with $b$ the breaking parameter.

  The vertical pressure gradient is added to the vertical velocity as 
  $$
  w^{n + 1}_l = w^{\star}_l - \Delta t \frac{[\phi]_l}{h^{n+1}_l}
  $$
  and the vertical velocity is limited by $w_\text{max}$ as 
  $$
  w^{n + 1}_l \leftarrow \text{sign} (w^{n + 1}_l) 
  \min \left( | w^{n + 1}_l |, w_\text{max} \right)
  $$
  */

  foreach() {
    double wmax = HUGE;
    if (breaking < HUGE) {
      wmax = 0.;
      foreach_layer()
	wmax += h[];
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
    foreach_layer()
      if (h[] > dry) {
	if (point.l == nl - 1)
	  w[] += dt*phi[]/h[];
	else
	  w[] -= dt*(phi[0,0,1] - phi[])/h[];
	if (fabs(w[]) > wmax)
	  w[] = (w[] > 0. ? 1. : -1.)*wmax;
      }

    /**
    The r.h.s. for $\eta^{n+1}$ is updated. It will be used in a
    second pass (neglecting the non-hydrostatic terms) in the
    [semi-implicit free-surface solver](implicit.h). This should be a
    small correction which is only necessary to limit the accumulation
    of divergence noise for long integration times. */

    foreach_dimension()
      rhs_eta[] += theta_H*sq(dt)*(su.x[1] - su.x[])/(Delta*cm[]);
  }
}

/**
## Cleanup

The *w* and *phi* fields are freed. */
      
event cleanup (i = end, last) {
  delete ({w, phi});
}
