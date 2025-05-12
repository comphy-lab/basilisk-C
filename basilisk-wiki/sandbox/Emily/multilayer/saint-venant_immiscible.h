/**
# A multilayer immiscible solver for the Saint-Venant equations

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

scalar zb[], hall[], eta[];  //Define hall[] as sum of hhl[]
vector u[]; //Define u[] as 'surface velocity' i.e. velocity at top layer

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
#include "multilayer_immiscible.h"

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
  vector * uol = (vector *) &output[nl];
  
  // new fields in ho[], uo[]
  foreach() {
    eta[] = zb[];
    hall[] = 0.;
    for (int l = 0; l < nl; l++) {
      scalar hi = input[l];
      scalar ho = output[l];
      scalar dh = updates[l];
      double hold = hi[];
      ho[] = hold + dt*dh[];
      vector uo = vector(output[nl + dimension*l]);
      if (ho[] > dry) {
	vector ui = vector(input[nl + dimension*l]);
	vector dhu = vector(updates[nl + dimension*l]);
	foreach_dimension()
	  uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
      }
      else { // dry
        foreach_dimension()
	  uo.x[] = 0.;
      }  
      hall[]+=ho[];
    }
    eta[]+=hall[];
  }
  /**
     In the case of [multiple layers](multilayer.h#viscous-friction-between-layers) we add the
     viscous friction between layers.
     Things are actually a little more subtle than this. If a layer is dry then you would want to use the values from the layer on the other side of it..... I will have to think a bit harder about how to deal with that - perhaps I can fill up hcol in a clever way?*/
  // This will not be correct if some of the layers are dry - it will currently break if any of them are dry hence check... What needs to happen in the case of some dry layers is fairly subtle...
  if (nl > 1)
    foreach(){
      bool iswet=true;
      for (int l = 0; l < nl; l++){
	if (iswet) {
	  scalar ho = output[l];
	  if (ho[] < dry) {
	    iswet = false;
	    break;
	  }
	}
      }
      if (iswet)
	vertical_viscosity (point, output, uol, dt);
    }
// fixme: on trees eta is defined as eta = zb + h and not zb +
  // ho in the refine_eta() and restriction_eta() functions below
  //scalar * list = list_concat ( {h,eta}, (scalar *) uol);
  //list = list_concat ( (scalar *) hol, list);
  //boundary (list);
  //free (list);
  boundary(all);
}

/**
When using an adaptive discretisation (i.e. a tree)., we need
to make sure that $\eta$ is maintained as $z_b + h$ whenever cells are
refined or restrictioned. */
#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + hall[];
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[] + hall[];
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
  
  vector u = vector(evolving[nl]);
  /**
  *Fh* and *Fq* will contain the fluxes for $h$ and $h\mathbf{u}$
  respectively and *S* is necessary to store the asymmetric topographic
  source term. */
  
  face vector Fh[], S[];
  tensor Fq[];
  
  /**
     The gradients are stored in locally-allocated fields. First-order
  reconstruction is used for the gradient fields. */
 
  vector geta[];
  tensor gu[];
  for (scalar s in {geta,gu}) {
    s.gradient = zero;
#if TREE
    s.prolongation = refine_linear;
#endif
  }
  gradients ({eta, u}, {geta, gu});
  /** We calculate and store gh for all the layers so we can access this for calculating the source term.*/
  vector * ghl = NULL;
  
  for (int l = 0; l < nl;l++) {
    scalar h0 = evolving[l];
    vector gh0 = new vector;
    for (scalar s in {gh0}) {
      s.gradient = zero;
#if TREE
      s.prolongation = refine_linear;
#endif
    }
    gradients ({h0}, {gh0});
    ghl = vectors_append (ghl, gh0);
  }
  /**
     We go through each layer. */
  
  for (int l = 0; l < nl; l++) {
    
    /**
    We recover the velocity field for the current layer and compute
    its gradient (for the first layer the gradient has already been
    computed above). 
    Question: How do we know where in gh and gu to store everything? or is it getting written over each time? It seems to be written over each time which means that you only have the gradients for that level - It would be good to have the gradients for each of the levels for reconstructing h_m in the
    */
    scalar h = evolving[l];
    vector u = vector (evolving[nl + dimension*l]);
    if (l > 0){
      for (scalar s in {gu}) {
	s.gradient = zero;
#if TREE
	s.prolongation = refine_linear;
#endif
      }
      gradients ((scalar *) {u}, (vector *) {gu});
    }
      /**
    The faces which are "wet" on at least one side are traversed. */

    foreach_face (reduction (min:dtmax)) {
      
      double hi = h[], hn = h[-1];
      if (hi > dry || hn > dry) {
	double halli=0., halln=0.;
	double sum_ghx=0.,sum_ghxn=0.;
	for (int m = l; m < nl; m++) {
	  vector gh = vector(ghl[m]);
	  scalar h = evolving[m];
	  sum_ghx+=gh.x[];
	  sum_ghxn+=gh.x[-1];
	  halli+=h[];
	  halln+=h[-1];
	}
	vector gh = vector(ghl[l]);
	  
	/**
	#### Left/right state reconstruction
      
	The gradients computed above are used to reconstruct the left
	and right states of the primary fields $h_l$, $\mathbf{u}_l$,
	$z_b$. The "interface" topography $z_{lr}$ is reconstructed
	using the hydrostatic reconstruction of [Audusse et al,
	2004](/src/references.bib#audusse2004) */
	double dx = Delta/2.;
	double zi = eta[] - halli;
	double zl = zi - dx*(geta.x[] - sum_ghx); 
	double zn = eta[-1] - halln;
	double zr = zn + dx*(geta.x[-1] - sum_ghxn);
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
	NEED to update this for nl levels...
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
	/**
	### Need to loop through all the different layers    
	 */
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	  /** 
	      The first version calculates the fluxes directly. I am now implementing a version where they are calculated and stored... We will see which is best...*/ 	
	for (int m = l+1; m < nl; m++) {
	  /**
	      This version calculates the fluxes directly. I am now implementing a version where they are calculated and stored... We will see which is best...
	      Not sure whether this will work here - I am calculating ghm[] and ghm[-1] directly rather than calculating the whole gradient vector each time. I may have the wrong Delta for the [-1] case - also not sure whether [-2] will work - it depends on whether the 5x5 stencil has been allocated or just the 3x3 ...
	  vector ghm0,ghmm1;
	  scalar hm = scalar(evolving[m]);
	  if (hm.gradient)
	    foreach_dimension() {
            #if EMBED
	      if (!fs.x[] || !fs.x[1]) {
		ghm0.x = 0.;
		ghmm1.x = 0.;
	      }
	      else {
            #endif
		ghm0.x[] = hm.gradient (hm[-1], hm[], hm[1])/Delta;
		ghmm1.x[] = hm.gradient (hm[-2], hm[-1], hm[])/Delta;
	      }
	    }
	  else // centered
	    foreach_dimension() {
            #if EMBED
	      if (!fs.x[] || !fs.x[1]) {
		ghm0.x = 0.;
		ghmm1.x = 0.;
	      }
	      else {
            #endif
		ghm0.x[] = (hm[1] - hm[-1])/(2.*Delta);
		ghmm1.x[] = (hm[] - hm[-2])/(2.*Delta);
	      }
	    }
	}
	}
          #if TREE
	  if (is_prolongation(cell)) {
	    hmi = coarse(hm);
	  }
	  if (is_prolongation(neighbor(-1))) {
	    hmn = coarse(hm,-1);
	  }
          #endif
	   */ 
	  scalar hhm = evolving[m];
	  vector gh = ghl[m];
	  double hhmi = hhm[];
	  double hhmn = hhm[-1];
	  double hhml = hhmi - dx*gh.x[]; 
	  double hhmr = hhmn + dx*gh.x[-1];
	  sl +=G/2.*(hl + hi)*min(rhol[m]/rhol[l],1.)*(hhmi - hhml);
	  sr +=G/2.*(hr + hn)*min(rhol[m]/rhol[l],1.)*(hhmn - hhmr);
	}
	/**
	#### Flux update */
      
	Fh.x[]   = fm.x[]*fh;
	Fq.x.x[] = fm.x[]*(fu - sl);
	S.x[]    = fm.x[]*(fu - sr);
	Fq.y.x[] = fm.x[]*fv;
      }
      else // dry
	Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;}

    boundary_flux ({Fh, S, Fq});

    /**
    #### Updates for evolving quantities
  
    We store the divergence of the fluxes in the update fields. Note that
    these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */
  
    scalar dh = updates[l];  
    vector dhu = vector(updates[nl + dimension*l]);
    foreach() {
      dh[] = - (Fh.x[1,0] - Fh.x[] + Fh.y[0,1] - Fh.y[])/(cm[]*Delta);
      foreach_dimension()
	dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);

      /**
      ALSO REMOVED FOR MY IMMISCIBLE VERSION
      For [multiple layers](multilayer.h#fluxes-between-layers) we
      need to store the divergence in each layer. */
 
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
      ALSO REMOVED FOR MY IMMISCIBLE VERSION
  For [multiple layers](multilayer.h#fluxes-between-layers) we need to
  add fluxes between layers. */
  delete ( (scalar*) ghl );
  return dtmax;
}

/**
## Initialisation and cleanup

We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0)
{
  assert (hhl == NULL && ul == NULL && wl == NULL);
  assert (nl > 0);
  rhol = malloc(nl * sizeof(double));
  mul = malloc(nl * sizeof(double));
  for (int l = 0; l < nl; l++) {
    scalar h = new scalar;
    scalar w = new scalar;
    vector u = new vector;
    hhl = list_append (hhl, h);
    if (l<nl-1)
      ul = vectors_append (ul, u);
    wl = list_append (wl, w);
    rhol[l]=1.;
    mul[l]=0.;
  }
  ul=vectors_append(ul,u);

  evolving = list_concat ( (scalar *) hhl, (scalar *) ul);
  foreach()
    for (scalar s in evolving)
      s[] = 0.;
  boundary (evolving);
  
  /**
  Do I want to set it so that by default hl is evenly spaced h at the start?*/
  /*
  layer = qmalloc (nl, double);
  for (int l = 0; l < nl; l++)
    layer[l] = 1./nl;
  */
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
  foreach(){
    hall[]=0;
    for (scalar hh in hhl)
      hall[]+=hh[];
    eta[] = zb[] + hall[];
  }
  boundary (all);
}

/**
At the end of the simulation, we free the memory allocated in the
*defaults* event. */

event cleanup (i = end, last) {
  free (evolving);
  free (hhl), hhl = NULL;
  free (ul), ul = NULL;
  free (wl), wl = NULL;
}

#include "elevation.h"
