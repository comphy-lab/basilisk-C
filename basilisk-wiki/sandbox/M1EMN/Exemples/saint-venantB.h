/**
# A solver for the Saint-Venant equations

This is the original file, only one change is done (multilayer-multilayerB)
The
[Saint-Venant equations](http://en.wikipedia.org/wiki/Shallow_water_equations)
can be written in integral form as the hyperbolic system of
conservation laws
$$
  \partial_t \int_{\Omega} \mathbf{q} d \Omega =
  \int_{\partial \Omega} \mathbf{f} (
  \mathbf{q}) \cdot \mathbf{n}d \partial
  \Omega - \int_{\Omega} hg \nabla z_b
$$
where $\Omega$ is a given subset of space, $\partial \Omega$ its boundary and
$\mathbf{n}$ the unit normal vector on this boundary. For
conservation of mass and momentum in the shallow-water context, $\Omega$ is a
subset of bidimensional space and $\mathbf{q}$ and
$\mathbf{f}$ are written
$$
  \mathbf{q} = \left(\begin{array}{c}
    h\\
    hu_x\\
    hu_y
  \end{array}\right), 
  \;\;\;\;\;\;
  \mathbf{f} (\mathbf{q}) = \left(\begin{array}{cc}
    hu_x & hu_y\\
    hu_x^2 + \frac{1}{2} gh^2 & hu_xu_y\\
    hu_xu_y & hu_y^2 + \frac{1}{2} gh^2
  \end{array}\right)
$$
where $\mathbf{u}$ is the velocity vector, $h$ the water depth and
$z_b$ the height of the topography. See also [Popinet, 
2011](/src/references.bib#popinet2011) for a more detailed
introduction.

## User variables and parameters

The primary fields are the water depth $h$, the bathymetry $z_b$ and
the flow speed $\mathbf{u}$. $\eta$ is the water level i.e. $z_b +
h$. Note that the order of the declarations is important as $z_b$
needs to be refined before $h$ and $h$ before $\eta$. */

scalar zb[], h[], eta[];
vector u[];

/**
The only physical parameter is the acceleration of gravity *G*. Cells
are considered "dry" when the water depth is less than the *dry*
parameter (this should not require tweaking). */

double G = 1.;
double dry = 1e-10;

/**
By default there is only a single layer i.e. this is the classical
Saint-Venant system above. This can be changed by setting *nl* to a
different value. The extension of the Saint-Venant system to multiple
layers is implemented in [multilayer.h](). */

int nl = 1;
/** 
--------------#include "multilayerB.h"
--------------------------------------------*/

 /**
# Multilayer Saint-Venant system with mass exchanges

*This is a modification of the original file* 


The [Saint-Venant system](saint-venant.h) is extended to multiple
layers following [Audusse et al, 2011](references.bib#audusse2011) as
$$
\partial_th + \partial_x\sum_{l=0}^{nl-1}h_lu_l = 0
$$
with
$$
h_l = \mathrm{layer}_lh
$$
with $\mathrm{layer}_l$ the relative thickness of the layers satisfying
$$
\mathrm{layer}_l >= 0,\;\sum_{l=0}^{nl - 1}\mathrm{layer}_l = 1.
$$
The momentum equation in each layer is thus
$$
\partial_t(h\mathbf{u}_l) + \nabla\cdot\left(h\mathbf{u}_l\otimes\mathbf{u}_l + 
\frac{gh^2}{2}\mathbf{I}\right) = 
- gh\nabla z_b + \frac{1}{\mathrm{layer}_l}\left[\mathbf{u}_{l+1/2}G_{l+1/2} - 
\mathbf{u}_{l-1/2}G_{l-1/2}
+ \nu\left(\frac{u_{l+1} - u_l}{h_{l+1/2}} - 
\frac{u_{l} - u_{l-1}}{h_{l-1/2}}\right)\right]
$$
where $G_{l+1/2}$ is the relative vertical transport velocity between
layers and the second term corresponds to viscous friction between
layers.

These last two terms are the only difference with the [one layer
system](saint-venant.h). 

The horizontal velocity in each layer is stored in *ul* and the
vertical velocity between layers in *wl*. */

vector * ul = NULL;
scalar * wl = NULL;
double * layer;

/**
## Viscous friction between layers

Boundary conditions on the top and bottom layers need to be added to close the
system for the viscous stresses. We chose to impose a Neumann condition on the
top boundary i.e.
$$
\partial_z u |_t = \dot{u}_t
$$
and a Navier slip condition on the bottom i.e.
$$
u|_b = u_b + \lambda_b \partial_z u|_b
$$
By default the viscosity is zero and we impose free-slip on the top
boundary and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
$\lambda_b = 0$, $u_b = 0$. */


/** $\nu B$ *is the yeld stress* */
double B ;
double nu ;

(const) scalar lambda_b = zeroc, dut = zeroc, u_b = zeroc;

/**
For stability, we discretise the viscous friction term implicitly as
$$
\frac{(hu_l)_{n + 1} - (hu_l)_{\star}}{\Delta t} =
\frac{\nu}{\mathrm{layer}_l}  \left( \frac{u_{l + 1} - u_l}{h_{l + 1 / 2}} -
\frac{u_l - u_{l - 1}}{h_{l - 1 / 2}} \right)_{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Mu}_{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal\tmrsub{m}atrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

double nu_eq(double,double);

void vertical_viscosity (Point point, double h, vector * ul, double dt)
{
  if (nu == 0.)
    return;
  
  double nueq[nl];
  double a[nl], b[nl], c[nl], rhs[nl];

  foreach_dimension() {

    /**
    The *rhs* of the tridiagonal system is $h_lu_l = h\mathrm{layer}_lu_l$. */
      
    int l = 0;
    for (vector u in ul)
      rhs[l] = h*layer[l]*u.x[], l++;
/**
  *definition of viscosity (it is a function of shear: Bingham flow,
  and the viscosity depends on pressure in the case of $\mu(I)$:
 
 for example, in the Bingham case, you define the function
   $$\tau = \nu (\frac{\partial u}{\partial z} + B)$$
   where $\nu$ is the viscosity and $\nu B$ is the yeld stress, we define an equivalent viscosity:
$$\nu_{eq}=\nu (1 + \frac{B}{|\frac{\partial u}{\partial z}|})$$
Note the regularization introduced to avoid division by zero*
 


*/
 double  press = G*h;//eta[]-zb[];
 for (l = 1; l < nl; l++) {
     vector um = ul[l-1] ;
     vector up = ul[l] ;
     press -= h*layer[l];
     double shear = (h>0? (up.x[]-um.x[])/(h*layer[l]) : HUGE);
     nueq[l] =  nu_eq(shear,press);
    }
    nueq[0] = nueq[1];
      

    /**
    The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
    $$
    a_{l > 0} = - \left( \frac{\nu \Delta t}{h_{l - 1 / 2}} \right)_{n + 1}
    $$
    $$
    c_{l < \mathrm{nl} - 1} = - \left( \frac{\nu \Delta t}{h_{l + 1 / 2}}
    \right)_{n + 1}
    $$
    $$
    b_{0 < l < \mathrm{nl} - 1} = \mathrm{layer}_l h_{n + 1} - a_l - c_l
    $$
    */
 
    for (l = 1; l < nl - 1; l++) {
      a[l] = - (nueq[l]+nueq[l-1])*dt/(h*(layer[l-1] + layer[l]));
      c[l] = - (nueq[l+1]+nueq[l])*dt/(h*(layer[l] + layer[l+1]));
      b[l] = layer[l]*h - a[l] - c[l];
    }
    
    /**
    For the top layer the boundary conditions give the (ghost)
    boundary value
    $$
    u_{\mathrm{nl}} = u_{\mathrm{nl} - 1} + \dot{u}_t h_{\mathrm{nl} - 1},
    $$
    which gives the diagonal coefficient and right-hand-side
    $$
    b_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1} h_{n + 1}
    - a_{\mathrm{nl} - 1}
    $$
    $$
    \mathrm{rhs}_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1}  
    (hu_{\mathrm{nl} - 1})_{\star} + \nu \Delta t \dot{u}_t
    $$
    */

    a[nl-1] = - (nueq[l-1]+nueq[l-2])*dt/(h*(layer[nl-2] + layer[nl-1]));
    b[nl-1] = layer[nl-1]*h - a[nl-1];
    rhs[nl-1] += nueq[l-1]*dt*dut[];

    /**
    For the bottom layer, the boundary conditions give the (ghost)
    boundary value $u_{- 1}$
    $$
    u_{- 1} = \frac{2 h_0}{2 \lambda_b + h_0} u_b + \frac{2 \lambda_b - h_0}{2
    \lambda_b + h_0} u_0,
    $$
    which gives the diagonal coefficient and right-hand-side
    $$
    b_0 = \mathrm{layer}_0 h_{n + 1} - c_0 + 
    \frac{2 \nu \Delta t}{2 \lambda_b + h_0}
    $$
    $$
    \mathrm{rhs}_0 = \mathrm{layer}_0  (hu_0)_{\star} + \frac{2 \nu \Delta t}{2
    \lambda_b + h_0} u_b
    $$
    */

    c[0] = - 2.*dt*nueq[0]/(h*(layer[0] + layer[1]));
    b[0] = layer[0]*h - c[0] + 2.*nueq[0]*dt/(2.*lambda_b[] + h*layer[0]);
    rhs[0] += 2.*nueq[0]*dt/(2.*lambda_b[] + h*layer[0])*u_b[];
    
    /**
    We can now solve the tridiagonal system using the [Thomas
    algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (l = 1; l < nl; l++) {
      b[l] -= a[l]*c[l-1]/b[l-1];
      rhs[l] -= a[l]*rhs[l-1]/b[l-1];
    }
    vector u = ul[nl-1];
    u.x[] = a[nl-1] = rhs[nl-1]/b[nl-1];
    for (l = nl - 2; l >= 0; l--) {
      u = ul[l];
      u.x[] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
    }
  }
}

/**
## Fluxes between layers

The relative vertical velocity between layers $l$ and $l+1$ is defined
as (eq. (2.22) of [Audusse et al, 2011](references.bib#audusse2011))
$$
G_{l+1/2} = \sum_{j=0}^{l}(\mathrm{div}_j + \mathrm{layer}_j\mathrm{dh})
$$
with
$$
\mathrm{div}_l = \nabla\cdot(h_l\mathbf{u}_l)
$$
$$
\mathrm{dh} = - \sum_{l=0}^{nl-1} \mathrm{div}_l
$$
*/

void vertical_fluxes (vector * evolving, vector * updates,
		      scalar * divl, scalar dh)
{
  foreach() {
    double Gi = 0., sumjl = 0.;
    for (int l = 0; l < nl - 1; l++) {
      scalar div = divl[l];
      Gi += div[] + layer[l]*dh[];
      sumjl += layer[l];
      scalar w = div;
      w[] = dh[]*sumjl - Gi;
      foreach_dimension() {

	/**
	To compute the vertical advection term, we need an estimate of
	the velocity at $l+1/2$. This is obtained using simple
	upwinding according to the sign of the interface velocity
	$\mathrm{Gi} = G_{l+1/2}$ and the values of the velocity in
	the $l$ and $l+1$ layers. Note that the inequality of
	upwinding is consistent with equs. (5.110) of [Audusse et al,
	2011](references.bib#audusse2011) and (77) of [Audusse et al,
	2011b](references.bib#audusse2011b) but not with eq. (2.23) of
	[Audusse et al, 2011](references.bib#audusse2011). */

	scalar ub = evolving[l].x, ut = evolving[l + 1].x;
	double ui = Gi < 0. ? ub[] : ut[];
	
	/**
	The flux at $l+1/2$ is then added to the updates of the bottom
	layer and substracted from the updates of the top layer. */
	
	double flux = Gi*ui;
	scalar du_b = updates[l].x, du_t = updates[l + 1].x;
	du_b[] += flux/layer[l];
	du_t[] -= flux/layer[l + 1];

	/**
	To compute the vertical velocity we use the definition of the
	mass flux term (eq. 2.13 of [Audusse et
	al, 2011](references.bib#audusse2011)):
	$$
	\mathrm{w}(\mathbf{x},z_{l+1/2}) = 
          \partial_t z_{l+1/2} - G_{l+1/2} + \mathbf{u}_{l+1/2}
	  \cdot \nabla z_{l+1/2}
	$$
	We can write the vertical position of the interface as:
	$$
	z_{l+1/2} = z_{b} + \sum_{j=0}^{l} h_{j}
	$$
        so that the vertical velocity is:
	$$
	\mathrm{w}(\mathbf{x},z_{l+1/2}) = 
          \mathrm{dh}\sum_{j=0}^{l}\mathrm{layer}_{j} - G_{l+1/2} + 
          \mathbf{u}_{l+1/2} \cdot \left[\nabla z_{b} + \nabla h 
				   \sum_{j=0}^{l}\mathrm{layer}_{j}\right]
	$$
	*/
	
	w[] += ui*((zb[1] - zb[-1]) + (h[1] - h[-1])*sumjl)/(2.*Delta);
      }
    }
  }
}

/** 
-------end of ---#include "multilayerB.h"
--------------------------------------------*/

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in *predictor-corrector.h* needs to
know which fields are updated. The list will be constructed in the
*defaults* event below. */

scalar * evolving = NULL;

/**
We need to overload the default *advance* function of the
predictor-corrector scheme, because the evolving variables ($h$ and
$\mathbf{u}$) are not the conserved variables $h$ and
$h\mathbf{u}$. */

trace
static void advance_saint_venant (scalar * output, scalar * input, 
				  scalar * updates, double dt)
{
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  vector * uol = (vector *) &output[1];

  // new fields in ho[], uo[]
  foreach() {
    double hold = hi[];
    ho[] = hold + dt*dh[];
    eta[] = zb[] + ho[];
    if (ho[] > dry) {
      for (int l = 0; l < nl; l++) {
        vector uo = vector(output[1 + dimension*l]);
      	vector ui = vector(input[1 + dimension*l]),
	  dhu = vector(updates[1 + dimension*l]);
	foreach_dimension()
	  uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
      }

      /**
      In the case of [multiple
      layers](multilayer.h#viscous-friction-between-layers) we add the
      viscous friction between layers. */
    
      if (nl > 1)
	vertical_viscosity (point, ho[], uol, dt);
    }
    else // dry
      for (int l = 0; l < nl; l++) {
        vector uo = vector(output[1 + dimension*l]);
        foreach_dimension()
	  uo.x[] = 0.;
      }
  }
    
  // fixme: on trees eta is defined as eta = zb + h and not zb +
  // ho in the refine_eta() and restriction_eta() functions below
  scalar * list = list_concat ({ho, eta}, (scalar *) uol);
  boundary (list);
  free (list);
}

/**
When using an adaptive discretisation (i.e. a tree)., we need
to make sure that $\eta$ is maintained as $z_b + h$ whenever cells are
refined or restrictioned. */

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h[];
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[] + h[];
}
#endif

/**
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "riemann.h"

trace
double update_saint_venant (scalar * evolving, scalar * updates, double dtmax)
{

  /**
  We first recover the currently evolving height and velocity (as set
  by the predictor-corrector scheme). */

  scalar h = evolving[0], dh = updates[0];
  vector u = vector(evolving[1]);
  
  /**
  *Fh* and *Fq* will contain the fluxes for $h$ and $h\mathbf{u}$
  respectively and *S* is necessary to store the asymmetric topographic
  source term. */

  face vector Fh[], S[];
  tensor Fq[];

  /**
  The gradients are stored in locally-allocated fields. First-order
  reconstruction is used for the gradient fields. */

  vector gh[], geta[];
  tensor gu[];
  for (scalar s in {gh, geta, gu}) {
    s.gradient = zero;
    #if TREE
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h, eta, u}, {gh, geta, gu});

  /**
  We go through each layer. */
  
  for (int l = 0; l < nl; l++) {

    /**
    We recover the velocity field for the current layer and compute
    its gradient (for the first layer the gradient has already been
    computed above). */
    
    vector u = vector (evolving[1 + dimension*l]);
    if (l > 0)
      gradients ((scalar *) {u}, (vector *) {gu});
    
    /**
    The faces which are "wet" on at least one side are traversed. */

    foreach_face (reduction (min:dtmax)) {
      double hi = h[], hn = h[-1];
      if (hi > dry || hn > dry) {

	/**
	#### Left/right state reconstruction
      
	The gradients computed above are used to reconstruct the left
	and right states of the primary fields $h$, $\mathbf{u}$,
	$z_b$. The "interface" topography $z_{lr}$ is reconstructed
	using the hydrostatic reconstruction of [Audusse et al,
	2004](/src/references.bib#audusse2004) */
      
	double dx = Delta/2.;
	double zi = eta[] - hi;
	double zl = zi - dx*(geta.x[] - gh.x[]);
	double zn = eta[-1] - hn;
	double zr = zn + dx*(geta.x[-1] - gh.x[-1]);
	double zlr = max(zl, zr);
	
	double hl = hi - dx*gh.x[];
	double up = u.x[] - dx*gu.x.x[];
	double hp = max(0., hl + zl - zlr);
	
	double hr = hn + dx*gh.x[-1];
	double um = u.x[-1] + dx*gu.x.x[-1];
	double hm = max(0., hr + zr - zlr);
	
	/**
	#### Riemann solver
	
	We can now call one of the approximate Riemann solvers to get
	the fluxes. */
	
	double fh, fu, fv;
	kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? u.y[-1] + dx*gu.y.x[-1] : u.y[] - dx*gu.y.x[])*fh;
	
	/**
	#### Topographic source term
      
	In the case of adaptive refinement, care must be taken to ensure
	well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */
	
        #if TREE
	if (is_prolongation(cell)) {
	  hi = coarse(h);
	  zi = coarse(zb);
	}
	if (is_prolongation(neighbor(-1))) {
	  hn = coarse(h,-1);
	  zn = coarse(zb,-1);
	}
        #endif
	
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	
	/**
	#### Flux update */
      
	Fh.x[]   = fm.x[]*fh;
	Fq.x.x[] = fm.x[]*(fu - sl);
	S.x[]    = fm.x[]*(fu - sr);
	Fq.y.x[] = fm.x[]*fv;
      }
      else // dry
	Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
    }

    boundary_flux ({Fh, S, Fq});

    /**
    #### Updates for evolving quantities
  
    We store the divergence of the fluxes in the update fields. Note that
    these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */
  
    vector dhu = vector(updates[1 + dimension*l]);
    foreach() {
      double dhl =
	layer[l]*(Fh.x[1,0] - Fh.x[] + Fh.y[0,1] - Fh.y[])/(cm[]*Delta);
      dh[] = - dhl + (l > 0 ? dh[] : 0.);
      foreach_dimension()
	dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);

      /**
      For [multiple layers](multilayer.h#fluxes-between-layers) we
      need to store the divergence in each layer. */
      
      if (l < nl - 1) {
	scalar div = wl[l];
	div[] = dhl;
      }

      /**
      We also need to add the metric terms. They can be written (see
      eq. (8) of [Popinet, 2011](references.bib#popinet2011)) 
      $$
      S_g = h \left(\begin{array}{c}
      0\\
      \frac{g}{2} h \partial_{\lambda} m_{\theta} + f_G u_y\\
      \frac{g}{2} h \partial_{\theta} m_{\lambda} - f_G u_x
      \end{array}\right)
      $$
      with
      $$
      f_G = u_y \partial_{\lambda} m_{\theta} - u_x \partial_{\theta} m_{\lambda}
      $$
      */

      double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
      double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
      double fG = u.y[]*dmdl - u.x[]*dmdt;
      dhu.x[] += h[]*(G*h[]/2.*dmdl + fG*u.y[]);
      dhu.y[] += h[]*(G*h[]/2.*dmdt - fG*u.x[]);
    }
  }

  /**
  For [multiple layers](multilayer.h#fluxes-between-layers) we need to
  add fluxes between layers. */

  if (nl > 1)
    vertical_fluxes ((vector *) &evolving[1], (vector *) &updates[1], wl, dh);
    
  return dtmax;
}

/**
## Initialisation and cleanup

We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0)
{
  assert (ul == NULL && wl == NULL);
  assert (nl > 0);
  ul = vectors_append (ul, u);
  for (int l = 1; l < nl; l++) {
    scalar w = new scalar;
    vector u = new vector;
    ul = vectors_append (ul, u);
    wl = list_append (wl, w);
  }

  evolving = list_concat ({h}, (scalar *) ul);
  foreach()
    for (scalar s in evolving)
      s[] = 0.;
  boundary (evolving);
  
  /**
  By default, all the layers have the same relative thickness. */

  layer = malloc (nl*sizeof(double));
  for (int l = 0; l < nl; l++)
    layer[l] = 1./nl;
  
  /**
  We overload the default 'advance' and 'update' functions of the
  predictor-corrector scheme and setup the prolongation and restriction
  methods on trees. */

  advance = advance_saint_venant;
  update = update_saint_venant;
#if TREE
  for (scalar s in {h,zb,u,eta}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
  }
  eta.refine  = refine_eta;
  eta.restriction = restriction_eta;
#endif
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  foreach()
    eta[] = zb[] + h[];
  boundary (all);
}

/**
At the end of the simulation, we free the memory allocated in the
*defaults* event. */

event cleanup (i = end, last) {
  free (evolving);
  free (layer);
  free (ul), ul = NULL;
  free (wl), wl = NULL;
}

#include "elevation.h"
