/**
# A solver for the Two-layer equations
The Saint-Venant equations for two layers:

The [Saint-Venant equations](http://en.wikipedia.org/wiki/Shallow_water_equations)
can be written in integral form as the hyperbolic system of
conservation laws 
$$
\partial_t \int_{\Omega} \mathbf{q} d \Omega +
\int_{\partial \Omega} \mathbf{f} (
\mathbf{q}) \cdot \mathbf{n}d \partial
\Omega + \int_{\Omega} \mathbf{F} d \Omega= 0
$$
where $\Omega$ is a given subset of space, $\partial \Omega$ its boundary and
$\mathbf{n}$ the unit normal vector on this boundary. For
conservation of mass and momentum in the shallow-water context, $\Omega$ is a
subset of bidimensional space and $\mathbf{q}$ and
$\mathbf{f}$ are written 
$$
\mathbf{q} =  
h_1,	\;			 
h_1u_1_x,\;					 
h_1u_1_y,\;
h_2 ,\;
h_2u_2_x,\; 
h_2u_2_y
$$
\;\;\;\;\;\;
$$
\mathbf{f} (\mathbf{q}) = \left(\begin{array}{cc}
h_1u_1_x & h_1u_1_y\					\
h_1u_1_x^2 + \frac{1}{2} gh_1^2 & h_1u_1_xu_1_y\	\
h_1u_1_xu_1_y & h_1u_1_y^2 + \frac{1}{2} gh_1^2\	\
h_2u_2_x & h_2u_2_y\					\
h_2u_2_x^2 + \frac{1}{2} gh_2^2 & h_2u_2_xu_2_y\	\
h_2u_2_xu_2_y & h_2u_2_y^2 + \frac{1}{2} gh_2^2
\end{array}\right),
\;\;\;\;\;\;
\mathbf{F} = \left(\begin{array}{c}
0\								\
g h_1 \frac{\partial}{\partial x} \left( h_2 + z_b \right)\	\
g h_1 \frac{\partial}{\partial y} \left( h_2 + z_b \right)\	\
0\									\
g h_2 \frac{\partial}{\partial x} \left( z_b + \frac{\rho_1}{\rho_2} h_1 \right)\ \
g h_2 \frac{\partial}{\partial y} \left( z_b + \frac{\rho_1}{\rho_2} h_1 \right)
\end{array} \right),
$$
where $h_1$ the water depth, $h_2$ is the landslide thickness, $\mathbf{u}_1$ is the velocity vector of the water, $\mathbf{u}_2$ is the velocity vector of the landslide, and
$z_b$ the height of the topography. See also [Popinet, 
2011](/src/references.bib#popinet2011) for a more detailed
introduction.

## User variables and parameters

The primary fields are the water depth $h_1$, the landslide thickness $h_2$, the bathymetry $z_b$ and
the flow speeds of water and landslide respectively $\mathbf{u}_1$ and $\mathbf{u}_2$. $\eta$ is the water level i.e. $z_b +
h_1 + h_2$. Note that the order of the declarations is important as $z_b$
needs to be refined before $h_2$, before $h_1$, before $\eta$. */

scalar zb[], h2[], h1[], eta[],h[];
vector u2[], u1[];
//double default_sea_level=0.;

/**
The only physical parameter is the acceleration of gravity `G` and the densities of the two fluids. 
Cells are considered "dry" when the water depth is less than the `dry` parameter (this 
should not require tweaking). */

double G = 1.;
double dry = 1e-10;
double RHO1 = 1.;
double RHO2 = 2.;

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in predictor-corrector.h needs
to know which fields are updated. */

scalar * evolving = {h1, u1, h2, u2};

/**
We need to overload the default *advance* function of the
predictor-corrector scheme, because the evolving variables ($h_1$, $\mathbf{u}_1$, $h_2$ and
$\mathbf{u}_2$) are not the conserved variables $h_1$, $h_1\mathbf{u}_1$, $h_2$ and
$h_2\mathbf{u}_2$. */

static void advance_two_layer (scalar * output, scalar * input, 
			       scalar * updates, double dt)
{
  //  fprintf(stderr,"advance_two_layer\n");
  // recover scalar and vector fields from lists
  scalar h1i = input[0], h1o = output[0], dh1 = updates[0];
  scalar h2i = input[3], h2o = output[3], dh2 = updates[3];
  vector u1i = { input[1], input[2] }, u1o = { output[1], output[2] };
  vector dhu1 = { updates[1], updates[2] }, u2i = { input[4], input[5] };
  vector u2o = { output[4], output[5] }, dhu2 = { updates[4], updates[5] };
  // new fields in ho[], uo[]
  foreach() {
    double h1old = h1i[];
    double h2old = h2i[];
    h1o[] = h1old + dt*dh1[];
    h2o[] = h2old + dt*dh2[];
    eta[] = h1o[] + h2o[] + zb[];
    if (h2o[] > dry){
      //      fprintf(stderr,"H2 wet");
      foreach_dimension()
	u2o.x[] = (h2old*u2i.x[] + dt*dhu2.x[])/h2o[];}
    else{
      //      fprintf(stderr,"H2 dry");
      foreach_dimension()
	u2o.x[] = 0.;}
    if (h1o[] > dry){
      //      fprintf(stderr," H1 wet\n");
      foreach_dimension()
	u1o.x[] = (h1old*u1i.x[] + dt*dhu1.x[])/h1o[];}
    else{
      //      fprintf(stderr," H1 dry\n");
      foreach_dimension()
	u1o.x[] = 0.;}
  }
  //  fprintf(stderr,"boundary\n");
  //h1o.prolongation = refine_linear;
  boundary ({h2o, h1o, eta, u2o, u1o});
  //  fprintf(stderr,"done\n");
}

/**
When using an adaptive discretisation (i.e. a quadtree)., we need
to make sure that $\eta$ is maintained as $z_b + h_1 + h_2$ whenever cells are
refined or coarsened. */

/**
#if QUADTREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h1[] + h2[];
}

static void coarsen_eta (Point point, scalar eta)
{
  eta[] = zb[] + h1[] + h2[];
}
#endif
*/

/**
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "riemann.h"

double update_two_layer (scalar * evolving, scalar * updates, double dtmax)
{
  
/**
We first recover the currently evolving fields (as set by the
predictor-corrector scheme). */

  scalar h1 = evolving[0];
  vector u1 = { evolving[1], evolving[2] };
  scalar h2 = evolving[3];
  vector u2 = { evolving[4], evolving[5] };
   
/**
`Fh1`, `Fh2` `Fq1` and `Fq2` will contain the fluxes for $h_1$, $h_2$, $h_1\mathbf{u}_1$ and $h_2\mathbf{u}_2$
respectively and `S1` and `S2` are necessary to store the asymmetric topographic
source terms. */

  vector Fh1[], Fh2[], S1[], S2[];
  tensor Fq1[], Fq2[];
   
/**
The gradients are stored in locally-allocated fields. First-order
reconstruction is used for the gradient fields. */

  vector gh1[], gh2[], geta[];
  tensor gu1[], gu2[];
  for (scalar s in {gh1, gh2, geta, gu1, gu2}){
    s.gradient = zero;
#if QUADTREE   //Is this where the problem is?  Do I need to refine some other way?
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h1, h2, eta, u1, u2}, {gh1, gh2, geta, gu1, gu2});
  //fprintf(stdout,"%% updates\n");
/**
The faces which are "wet" on at least one side are traversed. 
First we see whether "wet" in bottom fluid, if so look for lake at rest solution
$h_2+z_b=C_2$ $h_1=C_1$
If $h_2$ is dry look for lake at rest solution
$h_2=0$ $h_1+z_b=C_0$
*/

  foreach_face (reduction (min:dtmax)) {
    // First the bottom layer
    double hi = h2[], hn = h2[-1,0];
    if (hi > dry || hn > dry) {
      //      fprintf(stderr,"H2 wet");
      
/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
      double dx = Delta/2.;
      double zi = eta[] - hi - h1[];
      double zl = zi - dx*(geta.x[] - gh2.x[]- gh1.x[]);
      double zn = eta[-1,0] - hn - h1[-1,0];
      double zr = zn + dx*(geta.x[-1,0] - gh2.x[-1,0] - gh1.x[-1,0]);
      double zlr = max(zl, zr);
      
      double hl = hi - dx*gh2.x[];
      double up = u2.x[] - dx*gu2.x.x[];
      double hp = max(0., hl + zl - zlr);
      
      double hr = hn + dx*gh2.x[-1,0];
      double um = u2.x[-1,0] + dx*gu2.x.x[-1,0];
      double hm = max(0., hr + zr - zlr);
      
/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

      double fh, fu, fv;
      kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
      fv = (fh > 0. ? u2.y[-1,0] + dx*gu2.y.x[-1,0] : u2.y[] - dx*gu2.y.x[])*fh;
      
/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
      if (is_prolongation(cell)) {
	hi = coarse(h2,0,0,0);
	zi = coarse(zb,0,0,0);
      }
      if (is_prolongation(neighbor(-1,0))) {
	hn = coarse(h2,-1,0,0);
	zn = coarse(zb,-1,0,0);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl+dx*RHO1/RHO2*gh1.x[]));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr-dx*RHO1/RHO2*gh1.x[-1,0]));
      
/**
#### Flux update */

      Fh2.x[]   = fm.x[]*fh;
      Fq2.x.x[] = fm.x[]*(fu - sl);
      S2.x[]    = fm.x[]*(fu - sr);
      Fq2.y.x[] = fm.x[]*fv;
      //fprintf(stdout,"2 %g %g %g %g %g\n",x,y,fh,fu,fv);
/**
#### Here we must calculate fluxes for $h_1$  */
      hi = h1[], hn = h1[-1,0];
      if (hi > dry || hn > dry) {
	//	fprintf(stderr," H1 wet\n");
/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
	double dx = Delta/2.;
	double zi = eta[] - hi;
	double zl = zi - dx*(geta.x[] - gh1.x[]);
	double zn = eta[-1,0] - hn;
	double zr = zn + dx*(geta.x[-1,0] - gh1.x[-1,0]);
	double zlr = max(zl, zr);
	
	double hl = hi - dx*gh1.x[];
	double up = u1.x[] - dx*gu1.x.x[];
	double hp = max(0., hl + zl - zlr);
	
	double hr = hn + dx*gh1.x[-1,0];
	double um = u1.x[-1,0] + dx*gu1.x.x[-1,0];
	double hm = max(0., hr + zr - zlr);
	
/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

	kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? u1.y[-1,0] + dx*gu1.y.x[-1,0] : u1.y[] - dx*gu1.y.x[])*fh;
	
/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
	if (is_prolongation(cell)) {
	  hi = coarse(h1,0,0,0);
	  zi = coarse(zb,0,0,0)+coarse(h2,0,0,0);
	}
	if (is_prolongation(neighbor(-1,0))) {
	  hn = coarse(h1,-1,0,0);
	  zn = coarse(zb,-1,0,0)+coarse(h2,-1,0,0);
	}
#endif
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	
/**
#### Flux update */

	Fh1.x[]   = fm.x[]*fh;
	Fq1.x.x[] = fm.x[]*(fu - sl);
	S1.x[]    = fm.x[]*(fu - sr);
	Fq1.y.x[] = fm.x[]*fv;
	//fprintf(stdout,"1 %g %g %g %g %g\n",x,y,fh,fu,fv);
      }
      else { // dry
	Fh1.x[] = Fq1.x.x[] = S1.x[] = Fq1.y.x[] = 0.;
	//	fprintf(stderr," H1 dry\n");
      }
    }   
    else {//h_2 is dry - Note that h_2 is not necessarily dry in the neighbouring cell...
      Fh2.x[] = Fq2.x.x[] = S2.x[] = Fq2.y.x[] = 0.;
      //      fprintf(stderr,"H2 dry");
/**
If $h_2=0$ then we must satisfy $h_1+z_b=C_0$
*/
      hi = h1[], hn = h1[-1,0];
      if (hi > dry || hn > dry) {
	//	fprintf(stderr," H1 wet\n");
/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
	double dx = Delta/2.;
	double zi = eta[] - hi;
	double zl = zi - dx*(geta.x[] - gh1.x[]);
	double zn = eta[-1,0] - hn;
	double zr = zn + dx*(geta.x[-1,0] - gh1.x[-1,0]);
	double zlr = max(zl, zr);
	
	double hl = hi - dx*gh1.x[];
	double up = u1.x[] - dx*gu1.x.x[];
	double hp = max(0., hl + zl - zlr);
	
	double hr = hn + dx*gh1.x[-1,0];
	double um = u1.x[-1,0] + dx*gu1.x.x[-1,0];
	double hm = max(0., hr + zr - zlr);
	
/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

	double fh, fu, fv;
	kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? u1.y[-1,0] + dx*gu1.y.x[-1,0] : u1.y[] - dx*gu1.y.x[])*fh;
	
/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
	if (is_prolongation(cell)) {
	  hi = coarse(h1,0,0,0);
	  zi = coarse(zb,0,0,0);
	}
	if (is_prolongation(neighbor(-1,0))) {
	  hn = coarse(h1,-1,0,0);
	  zn = coarse(zb,-1,0,0);
	}
#endif
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	
/**
#### Flux update */

	Fh1.x[]   = fm.x[]*fh;
	Fq1.x.x[] = fm.x[]*(fu - sl);
	S1.x[]    = fm.x[]*(fu - sr);
	Fq1.y.x[] = fm.x[]*fv;
	//fprintf(stdout,"1 %g %g %g %g %g\n",x,y,fh,fu,fv);
      }
      else {// dry
	Fh1.x[] = Fq1.x.x[] = S1.x[] = Fq1.y.x[] = 0.;
	//	fprintf(stderr," H1 dry\n");
      }
    }
  }
  boundary_flux ({Fh1, Fh2, S1, S2, Fq1, Fq2});


  


/**
#### Updates for evolving quantities

We store the divergence of the fluxes in the update fields. Note that
these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */

  scalar dh1 = updates[0], dh2 = updates[3];
  vector dhu1 = { updates[1], updates[2] }, dhu2 = {updates[4], updates[5]};
  
  foreach() {
    dh1[] = (Fh1.x[] + Fh1.y[] - Fh1.x[1,0] - Fh1.y[0,1])/(cm[]*Delta);
    dh2[] = (Fh2.x[] + Fh2.y[] - Fh2.x[1,0] - Fh2.y[0,1])/(cm[]*Delta);
    foreach_dimension(){
      dhu1.x[] = (Fq1.x.x[] + Fq1.x.y[] - S1.x[1,0] - Fq1.x.y[0,1])/(cm[]*Delta);
      dhu2.x[] = (Fq2.x.x[] + Fq2.x.y[] - S2.x[1,0] - Fq2.x.y[0,1])/(cm[]*Delta);
    }

    /**
    We also need to add the metric terms. They can be written (see
    eq. (8) of [Popinet, 2011](references.bib#popinet2011)) 
    $$
    S_g = h \left(\begin{array}{c}
    0\								\
    \frac{g}{2} h \partial_{\lambda} m_{\theta} + f_G u_y\	\
    \frac{g}{2} h \partial_{\theta} m_{\lambda} - f_G u_x
    \end{array}\right)
    $$
    with
    $$
    f_G = u_y \partial_{\lambda} m_{\theta} - u_x \partial_{\theta} m_{\lambda}
    $$

#### CHECK THIS - NOT SURE OF IT ####
    */

    double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
    double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
    double fG1 = u1.y[]*dmdl - u1.x[]*dmdt;
    double fG2 = u2.y[]*dmdl - u2.x[]*dmdt;
    dhu1.x[] += h1[]*(G*h1[]/2.*dmdl + fG1*u1.y[]);
    dhu1.y[] += h1[]*(G*h1[]/2.*dmdt - fG1*u1.x[]);
    dhu2.x[] += h2[]*(G*h2[]/2.*dmdl + fG2*u2.y[]);
    dhu2.y[] += h2[]*(G*h2[]/2.*dmdt - fG2*u2.x[]);

  }
  scalar Fh1x[],Fh1y[],Fh2x[],Fh2y[],S1x[],S1y[],S2x[],S2y[];
  foreach() {
    Fh1x[]=val(Fh1.x,0,0);
    Fh1y[]=val(Fh1.y,0,0);
    Fh2x[]=val(Fh2.x,0,0);
    Fh2y[]=val(Fh2.y,0,0);
    S1x[]=val(S1.x,0,0);
    S1y[]=val(S1.y,0,0);
    S2x[]=val(S2.x,0,0);
    S2y[]=val(S2.y,0,0);
  }
  //printf ("%% Flux_updates t=%g\n",t);
  //output_field ({h1, h2, zb, Fh1x,Fh1y,Fh2x,Fh2y,S1x,S1y,S2x,S2y,dh1,dh2}, stdout, n = 1 << 10, linear = false);
  return dtmax;
}
/**
We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0)
{
  
/**
We overload the default 'advance' and 'update' function of the predictor-corrector
scheme and setup the refinement and coarsening methods on quadtrees. */

  advance = advance_two_layer; 
  update = update_two_layer;  
#if QUADTREE
  for (scalar s in {h1,zb,u2,u1,eta}) {
    s.refine = s.prolongation = refine_linear;
    s.coarsen = coarsen_volume_average;
  }
  //eta.refine  =refine_eta;
  //eta.prolongation = refine_linear;
  //eta.coarsen = coarsen_eta;
#endif
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  foreach()
    eta[] = zb[] + h1[] + h2[];
  boundary (all);
}


/**
## Conservation of water surface elevation 

When using the default adaptive reconstruction of variables, the
Saint-Venant solver will conserve the water depth when cells are
refined or coarsened. However, this will not necessarily ensure that
the "lake-at-rest" condition (i.e. a constant water surface elevation)
is also preserved. In what follows, we redefine the `refine()` and
`coarsen()` methods of the water depth $h$ so that the water surface
elevation $\eta$ is conserved. 

We start with the reconstruction of fine "wet" cells: */
#if QUADTREE
static void refine_elevation1 (Point point, scalar h1)
{
  struct { double x, y; } g; // gradient of eta
  if (h1[] >= dry) {
    double eta1 = zb[] + h1[] + h2[];   // water surface elevation 
    foreach_dimension()
      g.x = gradient (zb[-1,0] + h1[-1,0] + h2[-1,0], eta1, zb[1,0] + h1[1,0]+ h2[1,0])/4.;
    // reconstruct water depth h from eta and zb
    foreach_child()
      h1[] = max(0, eta1 + g.x*child.x + g.y*child.y - zb[]-h2[]);
  }
  else {
/**
The "dry" case is a bit more complicated. We look in a 3x3
neighborhood of the coarse parent cell and compute a depth-weighted
average of the "wet" surface elevation $\eta$. We need to do this
because we cannot assume a priori that the surrounding wet cells are
necessarily close to e.g. $\eta = 0$. */
    double v = 0., eta1 = 0.; // water surface elevation
    // 3x3 neighbourhood
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h1[i,j] >= dry) {
	  eta1 += h1[i,j]*(zb[i,j] + h2[i,j] + h1[i,j]);
	  v += h1[i,j];
	}
    if (v > 0.)
      eta1 /= v; // volume-averaged eta of neighbouring wet cells
    else
      /**
If none of the surrounding cells is wet, we assume a default sealevel
at zero. 
XXXXXXXXXXXX
Is this best?  What if zb is negative?  Should we set eta to that?*/
      eta1 = 0.;
/**
We then reconstruct the water depth in each child using $\eta$ (of the
parent cell i.e. a first-order interpolation in contrast to the wet
case above) and $z_b$ of the child cells. */
    // reconstruct water depth h from eta and zb
    foreach_child()
      h1[] = max(0, eta1 - zb[] - h2[]);
  }
}

/**
Cell coarsening is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */
static void coarsen_elevation1 (Point point, scalar h1)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h1[] > dry) {
      eta += h1[]*(zb[] + h1[] + h2[]);
      v += h1[];
    }  
/**
... and use this in combination with $z_b$ (of the coarse cell) to
compute the water depth $h$.  */
  if (v > 0.)
    h1[] = max(0., eta/v - zb[] - h2[]);
  else // dry cell
    h1[] = 0.;
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h1.refine  = refine_elevation1;
  h1.prolongation  = refine_linear;
  h1.coarsen = coarsen_elevation1;
}
#else // Cartesian
void conserve_elevation (void) {}
#endif

/**
## "Radiation" boundary conditions

This can be used to implement open boundary conditions at low
[Froude numbers](http://en.wikipedia.org/wiki/Froude_number). The idea
is to set the velocity normal to the boundary so that the water level
relaxes towards its desired value (`ref`). */

#define radiation1(ref) (sqrt (G*max(h1[],0.)) - sqrt(G*max((ref) - zb[], 0.)))
#define radiation2(ref) (sqrt (G*max(h2[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

/**
## Tide gauges

An array of `Gauge` structures passed to `output_gauges()` will create
a file (called `name`) for each gauge. Each time `output_gauges()` is
called a line will be appended to the file. The line contains the time
and the value of each scalar in `list` in the (wet) cell containing
`(x,y)`. The `desc` field can be filled with a longer description of
the gauge. */

  typedef struct {
    char * name;
    double x, y;
    char * desc;
    FILE * fp;
  } Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  for (Gauge * g = gauges; g->name; g++) {
    if (!g->fp) {
      g->fp = fopen (g->name, "w");
      if (g->desc)
	fprintf (g->fp, "%s\n", g->desc);
    }
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    Point point = locate (xp, yp);
    if (point.level >= 0 && h1[] > dry) {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " %g", s[]);
    }
    else {
      fprintf (g->fp, "%g", t);
      for (scalar s in list)
	fprintf (g->fp, " NaN");
    }
    fputc ('\n', g->fp);
    fflush (g->fp);
  }
}
