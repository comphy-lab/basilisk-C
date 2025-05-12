/**
# Non-hydrostatic extension of the multilayer solver

This adds the non-hydrostatic terms of the [vertically-Lagrangian
multilayer solver for free-surface flows](hydro.h) described in
[Popinet, 2019](/Bibliography#popinet2019). The corresponding system
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

The additional $w_k$ and $\phi_k$ fields are stored in two additional
lists: *wl* and *phil*. The convergence statistics of the multigrid
solver are stored in *mgp*. 

Wave breaking is parameterised usng the *breaking* parameter, which is
turned off by default (see Section 3.6.4 in [Popinet,
2019](/Bibliography#popinet2019)). */

#define NH 1
#include "poisson.h"

scalar * wl = NULL, * phil = NULL;
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
  
  assert (wl == NULL && phil == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar w = new scalar;
    scalar phi = new scalar;
    wl = list_append (wl, w);
    phil = list_append (phil, phi);
  }
  reset (wl, 0.);
  reset (phil, 0.);

  int l = 0;
  for (scalar w in wl) {
    tracers[l] = list_append (tracers[l], w);
    l++;
  }
}

/**
## Viscous term

Vertical diffusion is added to the vertical component of velocity
$w$. */

event viscous_term (i++)
{
  if (nu > 0.) {
    // fixme: ugly hack
    scalar lb = lambda_b, d = dut, u = u_b;
    lambda_b = dut = u_b = zeroc;
    foreach()
      vertical_viscosity (point, hl, (scalar *) wl, dt);
    boundary (wl);
    lambda_b = lb; dut = d; u_b = u;
  }
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
  h_l \partial_x \partial_x (h_l \phi_{l - 1 / 2}) + 
  h_l \partial_x \partial_x (h_l \phi_{l + 1 / 2}) +\\
  & 4 (\phi_{l + 1 / 2} - \phi_{l - 1 / 2}) + 8 h_l \sum^{l - 1}_{k = 0}
  (- 1)^{l + k}  \frac{\phi_{k + 1 / 2} - \phi_{k - 1 / 2}}{h_k}
\end{aligned}
$$
where $\mathbf{\phi}$ is the vector of $\phi_l$ for this column and
$\mathbf{d}$ is a vector dependent only on the values of $\phi$ in the
neighboring columns. */

static void box_matrix (Point point, scalar * phil, scalar * rhsl,
			double * H, double * d)
{
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    scalar h = hl[m], rhs = rhsl[m], phi = phil[m];
    d[l] = rhs[];
    foreach_dimension()
      d[l] -= h[]*(h[-1]*phi[-1] + h[1]*phi[1])/sq(Delta);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    H[l*nl + l] = - 2.*dimension*sq(h[])/sq(Delta) - 4.;
    if (l > 0) {
      scalar phip = phil[m+1];
      foreach_dimension()
	d[l] -= h[]*(h[-1]*phip[-1] + h[1]*phip[1])/sq(Delta);
      H[l*(nl + 1) - 1] = - 2.*dimension*sq(h[])/sq(Delta) + 4.;
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      scalar hk = hl[nl-1-k];
      if (hk[] > dry) {
	H[l*nl + k] -= 8.*s*h[]/hk[];
	H[l*nl + k - 1] += 8.*s*h[]/hk[];
      }
    }
  }
}

/**
## Relaxation operator

The updated values of $\phi$ in a column are obtained as
$$
\mathbf{\phi} = \mathbf{H}^{-1}\mathbf{b}
$$
were $\mathbf{H}$ and $\mathbf{b}$ are the Hessenberg matrix and
vector constructed by the function above. */

#include "hessenberg.h"

trace
static void relax_nh (scalar * phil, scalar * rhsl, int lev, void * data)
{
  foreach_level_or_leaf (lev) {
    double H[nl*nl], b[nl];
    box_matrix (point, phil, rhsl, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
    for (scalar phi in phil)
      phi[] = b[l--];
  }
}

/**
## Residual computation

The residual is computed as
$$
\begin{aligned}
  \text{res}_l = & \text{rhs}_l -
  h_l \partial_x \partial_x (h_l \phi_{l - 1 / 2}) -
  h_l \partial_x \partial_x (h_l \phi_{l + 1 / 2}) -\\
  & 4 (\phi_{l + 1 / 2} - \phi_{l - 1 / 2}) - 8 h_l \sum^{l - 1}_{k = 0}
  (- 1)^{l + k}  \frac{\phi_{k + 1 / 2} - \phi_{k - 1 / 2}}{h_k}
\end{aligned}
$$
*/

trace
static double residual_nh (scalar * phil, scalar * rhsl,
			   scalar * resl, void * data)
{
  double maxres = 0.;
#if 0 // TREE
  /* conservative coarse/fine discretisation (2nd order) */
  assert (false);
  face vector g[];
  foreach_face()
    g.x[] = face_gradient_x (a, 0);
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    scalar phi, rhs, res, h;
    int l = 0;
    for (phi,rhs,res,h in phil,rhsl,resl,hl) {
      res[] = rhs[] + 4.*phi[];
      foreach_dimension()
	res[] -= h[]*(h[1]*phi[1] - 2.*h[]*phi[] + h[-1]*phi[-1])/sq(Delta);
      if (l < nl - 1) {
        scalar phip = phil[l+1];
        res[] -= 4.*phip[];
	foreach_dimension()
	  res[] -= h[]*(h[1]*phip[1] - 2.*h[]*phip[] + h[-1]*phip[-1])/sq(Delta);
      }
      for (int k = l - 1, s = -1; k >= 0; k--, s = -s) {
	scalar hk = hl[k];
	if (hk[] > dry) {
	  scalar phik = phil[k], phikp = phil[k+1];
	  res[] += 8.*s*(phik[] - phikp[])*h[]/hk[];
	}
      }
#endif // !TREE    
#if EMBED
      if (p->embed_flux) {
	double c, e = p->embed_flux (point, a, alpha, &c);
	res[] += c - e*a[];
      }
#endif // EMBED    
      if (fabs (res[]) > maxres)
	maxres = fabs (res[]);
      l++;
    }
  }
  boundary (resl);
  return maxres;
}

/**
## Pressure solution

We first define a slope-limiting function and then compute the
right-hand-side of the Poisson equation as
$$
\text{rhs}_l = \frac{2 h_l}{\Delta t}  \left( \partial_x (hu)^{\star}_l 
  - [u^{\star} \partial_x z]_l +
    2 w^{\star}_l + 4 \sum^{l - 1}_{k = 0} (- 1)^{l + k} w^{\star}_k \right)
$$
*/

#define slope_limited(dz) (fabs(dz) < 1. ? (dz) : ((dz) > 0. ? 1. : -1.))
 
event pressure (i++)
{
  double h1 = 0., v1 = 0.;
  scalar * rhsl = list_clone (phil);
  foreach (reduction(+:h1) reduction(+:v1)) {
    scalar rhs, h, w;
    vector uf;
    int l = 0;
    coord dz;
    foreach_dimension()
      dz.x = zb[1] - zb[-1];
    for (rhs,h,w,uf in rhsl,hl,wl,ufl) {
      rhs[] = 2.*h[]*w[];
      foreach_dimension()
	rhs[] -= h[]*(uf.x[] + uf.x[1])/2.*
	slope_limited((dz.x + h[1] - h[-1])/(2.*Delta));
      if (l > 0) {
	vector ufm = ufl[l-1];
	foreach_dimension()
	  rhs[] += h[]*(ufm.x[] + ufm.x[1])/2.*slope_limited(dz.x/(2.*Delta));
      }
      foreach_dimension()
	rhs[] += h[]*((h[] + h[1])*uf.x[1] - (h[] + h[-1])*uf.x[])/(2.*Delta);
      for (int k = l - 1, s = -1; k >= 0; k--, s = -s) {
	scalar uk = wl[k];
	rhs[] += 4.*s*h[]*uk[];
      }
      rhs[] *= 2./dt;
      foreach_dimension()
	dz.x += h[1] - h[-1];
      h1 += dv()*h[];
      v1 += dv();
      l++;
    }
  }

  /**
  We then call the multigrid solver, using the relaxation and residual
  functions defined above, to get the non-hydrostatic pressure
  $\phi$. */
  
  restriction (hl);
  mgp = mg_solve (phil, rhsl, residual_nh, relax_nh, NULL,
		  nrelax = 4, res = NULL, minlevel = 1,
		  tolerance = TOLERANCE*sq(h1/(dt*v1)));
  delete (rhsl), free (rhsl);

  /**
  The non-hydrostatic pressure gradient is added to the acceleration
  and to the face velocities as 
  $$
  \begin{aligned}
  \alpha_{i + 1 / 2, l} \leftarrow & 2 \frac{\partial_x
  (h \phi)_{i + 1 / 2, l} - [\phi \partial_x z]_{i + 1 / 2,
  l}}{h_{i + 1, l}^{n + 1} + h_{i, l}^{n + 1}}\\
  u_{i + 1 / 2, l} \leftarrow & u_{i + 1 / 2, l} - \alpha_{i + 1 / 2, l}\\
  a_{i + 1 / 2, l} \leftarrow & a_{i + 1 / 2, l} - \alpha_{i + 1 / 2, l}
  \Delta t
  \end{aligned}
  $$
  */
  
  foreach_face() {
    double H = 0., Hm = 0.;
    for (scalar h in hl)
      H += h[], Hm += h[-1];
    if ((H > dry && Hm > dry) ||
	(H > dry && eta[] >= zb[-1]) ||
	(Hm > dry && eta[-1] >= zb[])) {
      scalar phi, h;
      vector uf, a;
      int l = 0;
      double dz = zb[] - zb[-1];
      for (phi,h,uf,a in phil,hl,ufl,al) {
	if (h[] + h[-1] > dry) {
	  double ax;
#if 0 // metric terms (do not seem to work yet)
	  if (l == nl - 1)
	    ax = sq(fm.x[])*(h[]*phi[] - h[-1]*phi[-1] +
			     (phi[] + phi[-1])*dz)
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  else {
	    scalar phip = phil[l+1];
	    ax = sq(fm.x[])*(h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1])
			     - ((phip[] + phip[-1])*(dz + h[] - h[-1]) -
				(phi[] + phi[-1])*dz))
	      /((cm[] + cm[-1])*Delta/2.*(h[] + h[-1]));
	  }
#else
	  if (l == nl - 1)
	    ax = fm.x[]*((h[]*phi[] - h[-1]*phi[-1])/Delta +
			 (phi[] + phi[-1])*slope_limited(dz/Delta))
	      /(h[] + h[-1]);
	  else {
	    scalar phip = phil[l+1];
	    ax = fm.x[]*
	      ((h[]*(phi[] + phip[]) - h[-1]*(phi[-1] + phip[-1]))/Delta
	       - ((phip[] + phip[-1])*slope_limited((dz + h[] - h[-1])/Delta) -
		  (phi[] + phi[-1])*slope_limited(dz/Delta)))
	      /(h[] + h[-1]);
	  }
#endif
	  uf.x[] -= dt*ax;
	  a.x[] -= ax;
	}
	l++, dz += h[] - h[-1];
      }
    }
  }
  boundary ((scalar *)ufl);
  boundary ((scalar *)al);

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
    double wmax = 0.;
    for (scalar h in hl)
      wmax += h[];
    wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    scalar phi, w, h;
    int l = 0;
    for (phi,w,h in phil,wl,hl) {
      if (h[] > dry) {
	if (l == nl - 1)
	  w[] += dt*phi[]/h[];
	else {
	  scalar phip = phil[l+1];
	  w[] -= dt*(phip[] - phi[])/h[];
	}
	if (fabs(w[]) > wmax)
	  w[] = (w[] > 0. ? 1. : -1.)*wmax;
      }
      l++;
    }
  }
  boundary (wl);
}

/**
## Cleanup

The *wl* and *phil* fields and lists are freed. */
      
event cleanup (i = end, last) {
  delete (wl), free (wl), wl = NULL;
  delete (phil), free (phil), phil = NULL;
}