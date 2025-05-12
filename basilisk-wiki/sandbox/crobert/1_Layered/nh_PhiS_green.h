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
2020](/Bibliography#popinet2020)). */

#define NH 1
#define blue 1
#define green1 0
#define green 0
#define greenGSL 0
#define redImproved 1
#include "poisson.h"

#if greenGSL
  #include <stdio.h>
  #include <gsl/gsl_linalg.h>
  #pragma autolink -lgsl -lgslcblas
#endif

scalar w, phi;
mgstats mgp;
double breaking = HUGE;
scalar phi_s[]; //Added for phi_s

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
  if (nu > 0.) {
    foreach()
      vertical_diffusion (point, (struct Vertical_Diffusion){h, w, dt, nu});
    boundary ({w});
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

static void box_matrix (Point point, scalar phi, scalar rhs,
			double * H, double * d)
{
  coord dz, ddz; //Added for blue and green terms
  foreach_dimension() {
    dz.x = (zb[1] - zb[-1])/2.;
    ddz.x = zb[1] - 2.*zb[0] + zb[-1];
  }
  foreach_layer()
    foreach_dimension() {
      dz.x += hf.x[1] - hf.x[];
      ddz.x += h[1] - 2.*h[0] + h[-1];
    }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    coord hx, hxx;
    foreach_dimension() {
      hx.x = hf.x[1,0,m] - hf.x[0,0,m];
      hxx.x = (h[1,0,m] - 2.*h[0,0,m] + h[-1,0,m]);
    }
    d[l] = rhs[0,0,m];
    foreach_dimension()
      d[l] -= h[0,0,m]*(h[-1,0,m]*phi[-1,0,m] + h[1,0,m]*phi[1,0,m])/sq(Delta)
        + blue*2.*h[0,0,m]*(phi[1,0,m] - phi[-1,0,m])/2.*(dz.x - hx.x)/sq(Delta)
        + green1*2.*h[0,0,m]*(phi[1,0,m] - phi[-1,0,m])/2.*(dz.x - hx.x)/sq(Delta);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    H[l*nl + l] = - 2.*dimension*sq(h[0,0,m])/sq(Delta) - 4.;
    foreach_dimension () 
      H[l*nl + l] += blue*2.*h[0,0,m]*(ddz.x - hxx.x)/sq(Delta);
    if (l < nl-1)
      foreach_dimension () {
        H[l*nl + l] -= green*(dz.x - hx.x)*(dz.x - 2*hx.x)/sq(Delta);
        H[l*(nl + 1) + 1] += green*(dz.x - hx.x)*(dz.x - 2*hx.x)/sq(Delta);
      }
    if (l > 0) {
      foreach_dimension()  
	d[l] -= h[0,0,m]*(h[-1,0,m]*phi[-1,0,m+1] + h[1,0,m]*phi[1,0,m+1])/sq(Delta)
          - blue*2*h[0,0,m]*(phi[1,0,m+1] - phi[-1,0,m+1])/2.*dz.x/sq(Delta)
          - green1*2.*h[0,0,m]*(phi[1,0,m+1] - phi[-1,0,m+1])/2.*dz.x/sq(Delta);
      H[l*(nl + 1) - 1] = - 2.*dimension*sq(h[0,0,m])/sq(Delta) + 4.;
        foreach_dimension () 
          H[l*(nl + 1) - 1] += - blue*2.*h[0,0,m]*ddz.x/sq(Delta);
      if (l > 1)
        foreach_dimension () {                     
          H[l*(nl + 1) - 1] -= green*dz.x*(dz.x + hx.x)/sq(Delta);  
          H[l*(nl + 1) - 2] += green*dz.x*(dz.x + hx.x)/sq(Delta);
        }
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = h[0,0,nl-1-k];
      if (hk > dry) {
	H[l*nl + k] -= 8.*s*h[0,0,m]/hk;
	H[l*nl + k - 1] += 8.*s*h[0,0,m]/hk;
      }
    }
    foreach_dimension() {
      dz.x -= hx.x;
      ddz.x -= hxx.x;
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
  scalar phi = phil[0], rhs = rhsl[0];
  foreach_level_or_leaf (lev) {
    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, H, b);
#if greenGSL 
  gsl_matrix_view m
    = gsl_matrix_view_array (H, nl, nl);
  gsl_vector_view d
    = gsl_vector_view_array (b, nl);
  gsl_vector *x = gsl_vector_alloc (nl);
  int s;
  gsl_permutation * p = gsl_permutation_alloc (nl);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_solve (&m.matrix, p, &d.vector, x);
  for(int l = 0; l < nl; l++)
    b[l] = gsl_vector_get(x,l);
  gsl_permutation_free (p);
  gsl_vector_free (x);
#else
    solve_hessenberg (H, b, nl);
#endif
    int l = nl - 1;
    foreach_layer()
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
  scalar phi = phil[0], rhs = rhsl[0], res = resl[0];
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
    coord dz, ddz; //Added for blue and green terms
    foreach_dimension() {
      dz.x = (zb[1] - zb[-1])/2.;
      ddz.x = zb[1] - 2*zb[] + zb[-1];
    }
    foreach_layer() {
      coord hx;
      foreach_dimension()
        hx.x = hf.x[1] - hf.x[]; //Added for blue and green terms     
      res[] = rhs[] + 4.*phi[];
      foreach_dimension()
	res[] -= h[]*(h[1]*phi[1] - 2.*h[]*phi[] + h[-1]*phi[-1])/sq(Delta) 
    + blue*2.*h[]*((phi[1] - phi[-1])/2.*dz.x + phi[]*ddz.x)/sq(Delta)
    + green1*2.*h[]*(phi[1] - phi[-1])/2.*dz.x/sq(Delta);
      if (point.l > 0)
        foreach_dimension()
          res[] -= green*(phi[0,0,-1] - phi[])*(dz.x)*(dz.x - hx.x)/sq(Delta);
      foreach_dimension() {
      dz.x += hx.x;
      ddz.x += h[1] - 2*h[] + h[-1];
      }
      if (point.l < nl - 1) {
        res[] -= 4.*phi[0,0,1];
	foreach_dimension()
	  res[] -= h[]*(h[1]*phi[1,0,1] - 2.*h[]*phi[0,0,1] + h[-1]*phi[-1,0,1])/sq(Delta)
      - blue*2.*h[]*((phi[1,0,1] - phi[-1,0,1])/2.*dz.x + phi[0,0,1]*ddz.x)/sq(Delta)
      - green1*2.*h[]*(phi[1,0,1] - phi[-1,0,1])/2.*dz.x/sq(Delta);
        if (point.l < nl - 2)
          foreach_dimension()
            res[] -= green*(phi[0,0,2] - phi[0,0,1])*dz.x*(dz.x + hx.x)/sq(Delta);
      }
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
	double hk = h[0,0,k];
	if (hk > dry)
	  res[] += 8.*s*(phi[0,0,k] - phi[0,0,k+1])*h[]/hk;
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
    }
  }
  boundary ({res});
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
The default maximum slope is set to 30 degrees. */


double max_slope = 0.577350269189626; // tan(30.*pi/180.)
#define slope_limited(dz) (fabs(dz) < max_slope ? (dz) :		\
			   ((dz) > 0. ? max_slope : - max_slope))
 
event pressure (i++)
{
  double h1 = 0., v1 = 0.;
  scalar rhs = new scalar[nl];
  foreach (reduction(+:h1) reduction(+:v1)) {
    coord dz, ddz;
    foreach_dimension() {
      dz.x = zb[1] - zb[-1];
      ddz.x = zb[1] - 2*zb[] + zb[-1]; //added for blue term
    }
    foreach_layer() {
      rhs[] = 2.*w[];
#if redImproved //Average over the horyzontal velocity
      foreach_dimension()
        rhs[] -= ((point.l == nl-1 ? u.x[] : u.x[0,0,1])/2.*
                  slope_limited((dz.x + h[1] - h[-1])/(2.*Delta)))
                  + u.x[]/2.*slope_limited((h[1] - h[-1])/(2.*Delta));
      if (point.l > 0)
        foreach_dimension()
          rhs[] += u.x[0,0,-1]/2.*slope_limited(dz.x/Delta);
#else
      foreach_dimension()
	rhs[] -= u.x[]*slope_limited((dz.x + h[1] - h[-1])/(2.*Delta));
      if (point.l > 0)
	foreach_dimension()
	  rhs[] += u.x[0,0,-1]*slope_limited(dz.x/(2.*Delta));
#endif
      foreach_dimension()
	rhs[] += ((h[] + h[1])*fm.x[1]*
		  (hf.x[1] > dry ? hu.x[1]/hf.x[1] : 0.) -
		  (h[] + h[-1])*fm.x[]*
		  (hf.x[] > dry ? hu.x[]/hf.x[] : 0.))/(2.*Delta);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
	rhs[] += 4.*s*w[0,0,k];
      rhs[] *= 2.*h[]/dt;
      foreach_dimension() {
        dz.x += h[1] - h[-1];
        ddz.x += h[1] - 2*h[] + h[-1]; //Added for blue term
      }
      h1 += dv()*h[];
      v1 += dv();
    }
    rhs[0,0,nl-1] -= 4.*phi_s[]; //Added for phi_s
    foreach_dimension() //Added for phi_s
      rhs[0,0,nl-1] -= h[0,0,nl-1]*(h[1,0,nl-1]*phi_s[1] - 2.*h[0,0,nl-1]*phi_s[] + h[-1,0,nl-1]*phi_s[-1])/sq(Delta) //Added for phi_s
        - blue*2.*h[]*((phi_s[1] - phi_s[-1])/2.*slope_limited(dz.x/(2*Delta)) + blue*phi_s[]*ddz.x/Delta)/Delta                                
        - green1*2.*h[]*(phi_s[1] - phi_s[-1])/2.*slope_limited(dz.x/(2*Delta))/Delta;
  }

  /**
  We then call the multigrid solver, using the relaxation and residual
  functions defined above, to get the non-hydrostatic pressure
  $\phi$. */
  
  restriction ({h});
  mgp = mg_solve ({phi}, {rhs}, residual_nh, relax_nh, NULL,
		  nrelax = 4, res = NULL, minlevel = 1,
		  tolerance = TOLERANCE*sq(h1/(dt*v1)));
  delete ({rhs});

  /**
  The non-hydrostatic pressure gradient is added to the face-weighted
  acceleration and to the face fluxes as
  $$
  \begin{aligned}
  \alpha_{i + 1 / 2, l} \leftarrow & \partial_x
  (h \phi)_{i + 1 / 2, l} - [\phi \partial_x z]_{i + 1 / 2, l}\\
  hu_{i + 1 / 2, l} \leftarrow & hu_{i + 1 / 2, l} - \alpha_{i + 1 / 2, l}\\
  ha_{i + 1 / 2, l} \leftarrow & ha_{i + 1 / 2, l} - \alpha_{i + 1 / 2, l}
  \Delta t
  \end{aligned}
  $$
  */
  
  foreach_face() {
    double dz = zb[] - zb[-1];
    foreach_layer() {
      if (hf.x[] > dry) {
	double ax;
	// fixme: metric terms are missing
#if 0
	if (point.l == nl - 1)
	  ax = fm.x[]*((h[]*phi[] - h[-1]*phi[-1])/Delta +
		       (phi[] + phi[-1])*slope_limited(dz/Delta))
	    /(h[] + h[-1]);
	else
	  ax = fm.x[]*
	    ((h[]*(phi[] + phi[0,0,1]) -
	      h[-1]*(phi[-1] + phi[-1,0,1]))/Delta
	     - ((phi[0,0,1] + phi[-1,0,1])*
		slope_limited((dz + h[] - h[-1])/Delta) -
		(phi[] + phi[-1])*slope_limited(dz/Delta)))
	    /(h[] + h[-1]);
	hu.x[] -= dt*hf.x[]*ax;
	a.x[] -= hf.x[]*ax;
#else
	if (point.l == nl - 1)
//	  ax = fm.x[]*((h[]*phi[] - h[-1]*phi[-1])/Delta +
//		       (phi[] + phi[-1])*slope_limited(dz/Delta))/2.; Modified for phi_s
	  ax = fm.x[]*((h[]*(phi[] + phi_s[]) -
			h[-1]*(phi[-1] + phi_s[-1]))/Delta
		       - ((phi_s[0] + phi_s[-1])*
			  slope_limited((dz + h[] - h[-1])/Delta) -
			  (phi[] + phi[-1])*slope_limited(dz/Delta)))/2.;                  
	else
	  ax = fm.x[]*((h[]*(phi[] + phi[0,0,1]) -
			h[-1]*(phi[-1] + phi[-1,0,1]))/Delta
		       - ((phi[0,0,1] + phi[-1,0,1])*
			  slope_limited((dz + h[] - h[-1])/Delta) -
			  (phi[] + phi[-1])*slope_limited(dz/Delta)))/2.;
	hu.x[] -= dt*ax;
	ha.x[] -= ax;
#endif
      }
      dz += h[] - h[-1];
    }
  }
  boundary ((scalar *){hu});
  boundary_flux ({ha});

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
    foreach_layer()
      wmax += h[];
    wmax = wmax > 0. ? breaking*sqrt((G == 0 ? 1. : G)*wmax) : 0.;
    foreach_layer() {
      if (h[] > dry) {
	if (point.l == nl - 1)
	  w[] -= dt*(phi_s[] - phi[])/h[]; //Modified for phi_s
	else
	  w[] -= dt*(phi[0,0,1] - phi[])/h[];
	if (fabs(w[]) > wmax)
	  w[] = (w[] > 0. ? 1. : -1.)*wmax;
      }
    }
  }
  boundary ({w});
}

/**
## Cleanup

The *w* and *phi* fields are freed. */
      
event cleanup (i = end, last) {
  delete ({w, phi});
}
