/**
# A solver for the Two-layer equations

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
\mathbf{q} = \left(\begin{array}{c}
h\\
h u_x\\
h u_y\\
h_s\\
h_s u_{sx}\\
h_s u_{sy}\\
\end{array}\right),
\;\;\;\;\;\;
\mathbf{f} (\mathbf{q}) = \left(\begin{array}{cc}
h   \ u_x & h \ u_y\\
h   \ u_x^2 + \frac{1}{2} g h^2 & h \ u_x u_y\\
h   \ u_x u_y & h u_y^2 + \frac{1}{2} g h^2\\
h_s \ u_{sx} & h_s \ u_{sy}\\
h_s \ u_{sx}^2 + \frac{1}{2} g h_s^2 & h_s \ u_{sx} u_{sy}\\
h_s \ u_{sx} u_{sy} & h_s \ u_{sy}^2 + \frac{1}{2} g h_s^2
\end{array}\right),
\;\;\;\;\;\;
\mathbf{F} = \left(\begin{array}{c}
0\\
g h \frac{\partial}{\partial x} \left( h_s + z_b \right)\\
g h \frac{\partial}{\partial y} \left( h_s + z_b \right)\\
0\\
g h_s \frac{\partial}{\partial x} \left( z_b + \frac{\rho}{\rho_s} h \right)\\
g h_s \frac{\partial}{\partial y} \left( z_b + \frac{\rho}{\rho_s} h \right)
\end{array}\right),
$$

where $h$ the water depth, $h_s$ is the landslide thickness, $\mathbf{u}$ is the velocity vector of the water, $\mathbf{u_s}$ is the velocity vector of the landslide, and
$z_b$ the height of the topography. See also [Popinet, 
2011](/src/references.bib#popinet2011) for a more detailed
introduction.

## User variables and parameters

The primary fields are the water depth $h$, the landslide thickness $h_s$, the bathymetry $z_b$ and
the flow speeds of water and landslide respectively $\mathbf{u}$ and $\mathbf{u_s}$. $\eta$ is the water level i.e. $z_b +
h + h_s$. Note that the order of the declarations is important as $z_b$
needs to be refined before $h_s$, before $h$, before $\eta$. */

scalar zb[], hs[], h[], eta[];
vector us[], u[];
double default_sea_level=0.;

/**
The only physical parameter is the acceleration of gravity `G` and the densities of the two fluids. 
Cells are considered "dry" when the water depth is less than the `dry` parameter (this 
should not require tweaking). */

double G = 1.;
double dry = 1e-10;
double RHOratio = 0.5;


/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in predictor-corrector.h needs
to know which fields are updated. */

scalar * evolving = NULL;

/**
We need to overload the default *advance* function of the
predictor-corrector scheme, because the evolving variables ($h$, $\mathbf{u}, $h_s$ and
$\mathbf{u_s}$) are not the conserved variables $h$, $h\mathbf{u}, $h_s$ and
$h_s\mathbf{u_s}$. */

trace
static void advance_two_layer (scalar * output, scalar * input, 
			       scalar * updates, double dt)
{
  //fprintf(stderr,"advance_two_layer\n");
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  scalar hsi = input[1], hso = output[1], dhs = updates[1];
  vector uo = vector(output[2]), uso = vector(output[2+dimension]);

/*
  Order within the scalar lists e.g. input, output, updates, evolving, [h, hs, u, us]
  Note that this differs from what I previously did...  But this should all occur within this so as long as it is self consistent it should be fine
  */

  // new fields in ho[], uo[]
  foreach() {
    double hold = hi[];
    double hsold = hsi[];
    ho[] = hold + dt*dh[];
    hso[] = hsold + dt*dhs[];
    eta[] = ho[] + hso[] + zb[];
    if (hso[] > dry){
      //      fprintf(stderr,"HS wet");
      vector usi = vector(input[2+dimension]);
      vector dhus = vector(updates[2+dimension]);
      foreach_dimension()
	uso.x[] = (hsold*usi.x[] + dt*dhus.x[])/hso[];}
    else{
      //fprintf(stderr,"HS dry");
      //hs[]=0.;  // Sometimes hs is negative - will this stop that or will it cause issues?
      foreach_dimension()
	uso.x[] = 0.;}
    if (ho[] > dry){
      //fprintf(stderr," H wet\n");
      vector ui = vector(input[2]);
      vector dhu = vector(updates[2]);
      foreach_dimension()
	uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];}
    else{
      //fprintf(stderr," H dry\n");
      //      h[]=0.;  // Sometime h is negative - will this stop that or will it cause issues?
      foreach_dimension()
	uo.x[] = 0.;}
  }
  //  fprintf(stderr,"boundary\n");
  //ho.prolongation = refine_linear;
  scalar * list = list_copy({hso, ho, eta, uso, uo});
  boundary (list);
  free (list);
  //  fprintf(stderr,"done\n");
}

/**
When using an adaptive discretisation (i.e. a quadtree)., we need
to make sure that $\eta$ is maintained as $z_b + h + h_s$ whenever cells are
refined or coarsened. */

#if TREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h[] + hs[];
}

static void restriction_eta (Point point, scalar eta)
{
  eta[] = zb[] + h[] + hs[];
}
#endif

/**
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "riemann.h"

trace
double update_two_layer (scalar * evolving, scalar * updates, double dtmax)
{
  
/**
We first recover the currently evolving fields (as set by the
predictor-corrector scheme). 
NOTE change in order as mentioned above*/

  scalar h = evolving[0], dh = updates[0];
  scalar hs = evolving[1], dhs = updates[1];
  vector u = vector(evolving[2]);
  vector us = vector(evolving[2+dimension]);
   
/**
`Fh`, `Fhs` `Fq` and `Fqs` will contain the fluxes for $h$, $h_s$, $h\mathbf{u}$ and $h_s\mathbf{u_s}$
respectively and `S1` and `S2` are necessary to store the asymmetric topographic
source terms. */

  face vector Fh[], Fhs[], S1[], S2[];
  tensor Fq[], Fqs[];
   
/**
The gradients are stored in locally-allocated fields. First-order
reconstruction is used for the gradient fields. */

  vector gh[], ghs[], geta[];
  tensor gu[], gus[];
  for (scalar s in {gh, ghs, geta, gu, gus}){
    s.gradient = zero;
#if TREE   //Is this where the problem is?  Do I need to refine some other way?
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h, hs, eta, u, us}, {gh, ghs, geta, gu, gus});
  //fprintf(stdout,"%% updates\n");
/**
The faces which are "wet" on at least one side are traversed. 
First we see whether "wet" in bottom fluid - if so look for lake at rest solution
$h_s+z_b=C_0$ $h=C$
If hs is dry look for lake at rest solution
$h_s=0$ $h+z_b=C_0$
*/

  int hswet;
  foreach_face (reduction (min:dtmax)) {
    // First the bottom layer
    double hi = hs[], hn = hs[-1];
    if (hi > dry || hn > dry) {
      //      fprintf(stderr,"HS wet");
	   hswet=1;
      
/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
      double dx = Delta/2.;
      double zi = eta[] - hi - h[];
      double zl = zi - dx*(geta.x[] - ghs.x[]- gh.x[]);
      double zn = eta[-1,0] - hn - h[-1,0];
      double zr = zn + dx*(geta.x[-1,0] - ghs.x[-1,0] - gh.x[-1,0]);
      double zlr = max(zl, zr);
      
      double hl = hi - dx*ghs.x[];
      double up = us.x[] - dx*gus.x.x[];
      double hp = max(0., hl + zl - zlr);
      
      double hr = hn + dx*ghs.x[-1,0];
      double um = us.x[-1,0] + dx*gus.x.x[-1,0];
      double hm = max(0., hr + zr - zlr);
      
/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

      double fh, fu, fv;
      kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
      fv = (fh > 0. ? us.y[-1,0] + dx*gus.y.x[-1,0] : us.y[] - dx*gus.y.x[])*fh;
      
/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if TREE
      if (is_prolongation(cell)) {
	hi = coarse(hs);
	zi = coarse(zb);
      }
      if (is_prolongation(neighbor(-1,0))) {
	hn = coarse(hs,-1);
	zn = coarse(zb,-1);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl+dx*RHOratio*gh.x[]));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr-dx*RHOratio*gh.x[-1,0]));
      
/**
#### Flux update */

      Fhs.x[]   = fm.x[]*fh;
      Fqs.x.x[] = fm.x[]*(fu - sl);
      S2.x[]    = fm.x[]*(fu - sr);
      Fqs.y.x[] = fm.x[]*fv;
      //fprintf(stdout,"2 %g %g %g %g %g\n",x,y,fh,fu,fv);

    }   
    else {//h_s is dry - Note that h_s is not necessarily dry in the neighbouring cell...
      Fhs.x[] = Fqs.x.x[] = S2.x[] = Fqs.y.x[] = 0.;
      hswet=0;
      //fprintf(stderr,"HS dry");
    }
/**
Now we must calculate fluxes for $h$.
If $h_s>dry$ then we must satisfy $h_s+z_b=C_0$ $h=C$;
If $h_s=0$ then we must satisfy $h+z_b=C_0$
*/
      hi = h[], hn = h[-1];
      if (hi > dry || hn > dry) {
	//	fprintf(stderr," H wet\n");
/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
	double dx = Delta/2.;
	double zi = eta[] - hi;
	double zl = zi - dx*(geta.x[] - gh.x[]);
	double zn = eta[-1,0] - hn;
	double zr = zn + dx*(geta.x[-1,0] - gh.x[-1,0]);
	double zlr = max(zl, zr);
	
	double hl = hi - dx*gh.x[];
	double up = u.x[] - dx*gu.x.x[];
	double hp = max(0., hl + zl - zlr);
	
	double hr = hn + dx*gh.x[-1,0];
	double um = u.x[-1,0] + dx*gu.x.x[-1,0];
	double hm = max(0., hr + zr - zlr);
	
/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

	double fh, fu, fv;
	kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
	fv = (fh > 0. ? u.y[-1,0] + dx*gu.y.x[-1,0] : u.y[] - dx*gu.y.x[])*fh;
	
/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if TREE
	if (is_prolongation(cell)) {
	  hi = coarse(h);
	  zi = hswet?coarse(zb) + coarse(hs):coarse(zb);
	}
	if (is_prolongation(neighbor(-1,0))) {
	  hn = coarse(h,-1);
	  zn = hswet?coarse(zb,-1) + coarse(hs,-1):coarse(zb,-1);
	}
#endif
	double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
	double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));
	
/**
#### Flux update */

	Fh.x[]   = fm.x[]*fh;
	Fq.x.x[] = fm.x[]*(fu - sl);
	S1.x[]    = fm.x[]*(fu - sr);
	Fq.y.x[] = fm.x[]*fv;
	//fprintf(stdout,"1 %g %g %g %g %g\n",x,y,fh,fu,fv);
      }
      else {// dry
	Fh.x[] = Fq.x.x[] = S1.x[] = Fq.y.x[] = 0.;
	//fprintf(stderr," H dry\n");
      }
    }
  boundary_flux ({Fh, Fhs, S1, S2, Fq, Fqs});

/**
#### Updates for evolving quantities

We store the divergence of the fluxes in the update fields. Note that
these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */

  //scalar dh = updates[0], dhs = updates[1];
  vector dhu = vector( updates[2]), dhus = vector(updates[2+dimension]);
  
  foreach() {
    dh[] = (Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/(cm[]*Delta);
    dhs[] = (Fhs.x[] + Fhs.y[] - Fhs.x[1,0] - Fhs.y[0,1])/(cm[]*Delta);
    foreach_dimension(){
      dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S1.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);
      dhus.x[] = (Fqs.x.x[] + Fqs.x.y[] - S2.x[1,0] - Fqs.x.y[0,1])/(cm[]*Delta);
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
    double fG1 = u.y[]*dmdl - u.x[]*dmdt;
    double fG2 = us.y[]*dmdl - us.x[]*dmdt;
    dhu.x[] += h[]*(G*h[]/2.*dmdl + fG1*u.y[]);
    dhu.y[] += h[]*(G*h[]/2.*dmdt - fG1*u.x[]);
    dhus.x[] += hs[]*(G*hs[]/2.*dmdl + fG2*us.y[]);
    dhus.y[] += hs[]*(G*hs[]/2.*dmdt - fG2*us.x[]);

  }

  return dtmax;
}


// XXXXXXXXXXXXXXXXXXXXUP TO HERE XXXXXXXXXXXXXXXXXXXXXXXXXXXX need to add in initialisation page 7 of new sain venant...
/**
We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0){
  evolving = list_copy({h, hs, u, us});
  foreach()
    for (scalar s in evolving)
      s[]=0.;
  boundary(evolving);
/**
We overload the default 'advance' and 'update' function of the predictor-corrector
scheme and setup the refinement and coarsening methods on quadtrees. */

  advance = advance_two_layer; 
  update = update_two_layer;  
#if TREE
  for (scalar s in {zb,hs,eta,h,us,u}) {
    s.refine = s.prolongation = refine_linear;
    s.restriction = restriction_volume_average;
  }
  eta.refine  =refine_eta;
  eta.restriction = restriction_eta;
  eta.prolongation = refine_linear;
#endif
}

/**
The event below will happen after all the other initial events to take
into account user-defined field initialisations. */

event init (i = 0)
{
  foreach()
    eta[] = zb[] + h[] + hs[];
  boundary (all);
}

event cleanup (i = end, last) {
  free (evolving);
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
#if TREE
static void refine_elevation1 (Point point, scalar h)
{
  //struct { double x, y; } g; // gradient of eta
  if (h[] >= dry) {
    //fprintf(stderr,"Before: zb=%g, hs=%g, h=%g\n",zb[],hs[],h[]);
    double eta1 = hs[] >= dry ? zb[] + h[] + hs[]:zb[] + h[];   // water surface elevation  
    //fprintf(stderr,"After: eta=%g, zb=%g, hs=%g, h=%g\n",eta1,zb[],hs[],h[]);
        coord g; // gradient of eta
    if (gradient)
      foreach_dimension()
    	g.x = gradient (zb[-1] + h[-1] + hs[-1], eta1, zb[1] + h[1]+ hs[1])/4.;
    else
      foreach_dimension()
	g.x = (zb[1] - zb[-1])/(2.*Delta);
    	//g.x = (zb[1] + hs[1] + h[1] - zb[-1] - hs[-1] - h[-1])/(2.*Delta);
    // reconstruct water depth h from eta and zb
    foreach_child() {
      double etac = eta1;
      foreach_dimension()
    	etac += g.x*child.x;
      h[] = max(0, etac - zb[] - hs[]);
    }
  }
  else {
/**
The "dry" case is a bit more complicated. We look in a 3x3
neighborhood of the coarse parent cell and compute a depth-weighted
average of the "wet" surface elevation $\eta$. We need to do this
because we cannot assume a priori that the surrounding wet cells are
necessarily close to e.g. $\eta = 0$. */
    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood - 1 symbolises this rather than 5x5
    foreach_neighbor(1)
      if (h[] >= dry) {
	eta += hs[] > dry ? h[]*(zb[] + h[] + hs[]) : h[]*(zb[] + h[]);
	v += h[];
	}
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else
      /**
If none of the surrounding cells is wet, we assume a default sealevel
at zero. */
      //eta = 0.;
      eta = default_sea_level;
      //eta = min( default_sea_level, zb[]);
    /**
We then reconstruct the water depth in each child using $\eta$ (of the
parent cell i.e. a first-order interpolation in contrast to the wet
case above) and $z_b$ of the child cells. */
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0., eta - zb[] - hs[]);
  }
}

/**
Cell coarsening is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */
static void restriction_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += hs[] > dry ? h[]*(zb[] + hs[] + h[]) :  h[]*(zb[] + h[]);
      v += h[];
    }
/**
... and use this in combination with $z_b$ (of the coarse cell) to
compute the water depth $h$.  */
  if (v > 0.)
    h[] = max(0., eta/v - zb[] - hs[]);
  else // dry cell
    h[] = 0.;
}

/**
We also need to define a consistent prolongation function. For cells
which are entirely surrounded by wet cells, we can use the standard
linear refinement function, otherwise we use straight injection from
the parent cell. */

static void prolongation_elevation (Point point, scalar h)
{
  bool wet = true;
  foreach_neighbor(1)
    if (h[] <= dry)
      wet = false, break;
  if (wet)
    refine_linear (point, h);
  else {
    //fprintf(stderr,"prolongation: x=%g, y=%g, h[]=%g, zb[]=%g\n",x,y,h[],val(zb));
    double hc = h[], zc = zb[];
    foreach_child() {
      h[] = hc;
      zb[] = zc;
    }
  }
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h.refine  = refine_elevation1;
  h.prolongation  = prolongation_elevation;
  h.restriction = restriction_elevation;
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

#define radiation1(ref) (sqrt (G*max(h[],0.)) - sqrt(G*max((ref) - zb[], 0.)))
#define radiation2(ref) (sqrt (G*max(hs[],0.)) - sqrt(G*max((ref) - zb[], 0.)))

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
  scalar * list1 = list_append (NULL, h);
  for (scalar s in list)
    list1 = list_append (list1, s);

  for (Gauge * g = gauges; g->name; g++);
  for (Gauge * g = gauges; g->name; g++) {
    if (!g->fp) {
      g->fp = fopen (g->name, "w");
      if (g->desc)
	fprintf (g->fp, "%s\n", g->desc);
    }
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    Point point = locate (xp, yp);
    if (point.level >= 0 && h[] > dry) {
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
