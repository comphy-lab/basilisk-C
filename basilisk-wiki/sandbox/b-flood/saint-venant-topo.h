/**
# A solver for the Saint-Venant equations (on steep slopes)

This solver is a variant (for steep slopes) of the [standard
Saint-Venant solver](/src/saint-venant.h).

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

#if !LAYERS
int nl = 1;
#endif

#include "multilayer.h"

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
refined or restricted. */

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

  vector gh[], geta[], gz[];
  tensor gu[];
  for (scalar s in {gh, geta, gu}) { // fixme: and gz
    s.gradient = zero;
    #if TREE
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h, eta, u, zb}, {gh, geta, gu, gz});

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
    Water slope reset */
    
    foreach()
      foreach_dimension()
        if (sq(geta.x[]) > sq(gz.x[]) &&
	    ((eta[] - h[] - Delta/2.*(geta.x[] - gh.x[]) > eta[-1] ||
	      eta[] - h[] + Delta/2.*(geta.x[] - gh.x[]) > eta[1])))	  
	  geta.x[] = gh.x[] + gz.x[];
    
    /**
    The faces which are "wet" on at least one side are traversed. */

    foreach_face (reduction (min:dtmax)) {
      double hi = h[], hn = h[-1];
      if (hi > dry || hn > dry) {

	/**
	#### Left/right state reconstruction on abrupt topography
	  
	We reconstruct the left and right state using the CN reconstruction 
	cite Buttinger. */
	  
	double dx = Delta/2.;	  
	double zi = zb[];
	double zn = zb[-1];	  
	double hr = hi - dx*gh.x[];  
	double hl = hn + dx*gh.x[-1];	  
	double etar = eta[] - dx*geta.x[];
	double etal = eta[-1] + dx*geta.x[-1];
	
	// We define the topography term at the interface p/m
	double zr = zi - dx*(gz.x[]);
	double zl = zn + dx*(gz.x[-1]);
	
	// Firstly, we define the Audusse terms
	double zA = max(zr,zl);
    	
	// Now the CN terms
	double zCN = min( zA, min(etal, etar));
	double hCNr = max(0.,min(etar - zCN, hr));
	double hCNl = max(0.,min(etal - zCN, hl ));
	
	/***
	#### Velocity reconstruction
	
	To avoid high velocities near dry cells, we reconstruct
	velocities according to Bouchut. */
	
	double ul, ur, vl, vr;
	if (hi > dry ){
	  ur = u.x[] - (1. + dx*gh.x[]/hi )*dx*gu.x.x[];
	  vr = u.y[] - (1. + dx*gh.x[]/hi )*dx*gu.y.x[];
	}
	else{
	  ur = u.x[] - dx*gu.x.x[];
	  vr = u.y[] - dx*gu.y.x[];
	}	  
	if (hn > dry ){
	  ul = u.x[-1] + (1. - dx*gh.x[-1]/hn)*dx*gu.x.x[-1];
	  vl = u.y[-1] + (1. - dx*gh.x[-1]/hn)*dx*gu.y.x[-1];
	}
	else{
	  ul = u.x[-1] + dx*gu.x.x[-1];
	  vl = u.y[-1] + dx*gu.y.x[-1];
	}
	
	/**
	#### Riemann solver
	
	We can now call one of the approximate Riemann solvers to get
	the fluxes. */
	
	double fh, fu, fv;
	hllc (hCNl, hCNr, ul, ur, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? vl : vr)*fh;
	
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
	
	double sl = G/2.*(hi + hCNr)*(zi  - zCN);
	double sr = G/2.*(hCNl + hn)*(zn - zCN);
	  
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
    foreach_dimension()
      u.x.l = l;
    w.l = l;
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

  layer = qmalloc (nl, double);
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
    //    s.restriction = restriction_volume_average;
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

/**
# "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (*ref*). */

#define radiation(ref) (sqrt (G*max(h[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

#include "elevation.h"
#include "gauges.h"

/**
## Link to the homepage

* [Homepage](Readme)
*/
