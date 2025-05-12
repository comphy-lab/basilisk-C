/**
# A solver for the Saint-Venant equations

The
[Saint-Venant equations](http://en.wikipedia.org/wiki/Shallow_water_equations)
can be written in integral form as the hyperbolic system ofwind_pool0.
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

scalar zb[], pa[], h[], eta[];
vector u[], ts[];
double Pa2m = 0.00009916;

/**
The only physical parameter is the acceleration of gravity `G`. Cells are 
considered "dry" when the water depth is less than the `dry` parameter (this 
should not require tweaking). */

double G = 1.;
double dry = 1e-10;

/**
## Time-integration

### Setup

Time integration will be done with a generic
[predictor-corrector](predictor-corrector.h) scheme. */

#include "predictor-corrector.h"

/**
The generic time-integration scheme in predictor-corrector.h needs
to know which fields are updated. */

scalar * evolving = {h, u};

/**
We need to overload the default *advance* function of the
predictor-corrector scheme, because the evolving variables ($h$ and
$\mathbf{u}$) are not the conserved variables $h$ and
$h\mathbf{u}$. */

static void advance_saint_venant (scalar * output, scalar * input, 
				  scalar * updates, double dt)
{
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  vector
    ui = { input[1], input[2] },
    uo = { output[1], output[2] },
    dhu = { updates[1], updates[2] };

  // new fields in ho[], uo[]
  foreach() {
    double hold = hi[];
    ho[] = hold + dt*dh[];
    eta[] = ho[] + zb[];
    if (ho[] > dry)
      foreach_dimension()
	uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
    else
      foreach_dimension()
	uo.x[] = 0.;
  }
  boundary ({ho, eta, uo});
}

/**
When using an adaptive discretisation (i.e. a quadtree)., we need
to make sure that $\eta$ is maintained as $z_b + h$ whenever cells are
refined or coarsened. */

#if QUADTREE
static void refine_eta (Point point, scalar eta)
{
  foreach_child()
    eta[] = zb[] + h[];
}

static void coarsen_eta (Point point, scalar eta)
{
  eta[] = zb[] + h[];
}
#endif

/**
We use the main time loop (in the predictor-corrector scheme) to setup
the initial defaults. */

event defaults (i = 0)
{

/**
We overload the default 'advance' function of the predictor-corrector
scheme and setup the refinement and coarsening methods on quadtrees. */

  advance = advance_saint_venant;  
#if QUADTREE
  eta.refine  = refine_eta;
  eta.coarsen = coarsen_eta;
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
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "riemann.h"

double update (scalar * evolving, scalar * updates, double dtmax)
{

/**
We first recover the currently evolving fields (as set by the
predictor-corrector scheme). */

  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };

/**
The gradients are stored in locally-allocated fields. First-order
reconstruction is used for the gradient fields. */

  vector gh[], geta[], gpa[];
  tensor gu[];
  for (scalar s in {gh, geta, gpa, gu})
    s.gradient = zero;
  gradients ({h, eta, pa, u}, {gh, geta, gpa, gu});

/**
`Fh` and `Fq` will contain the fluxes for $h$ and $h\mathbf{u}$
respectively and `S` is necessary to store the asymmetric topographic
source term. */

  vector Fh[], S[];
  tensor Fq[];

/**
The faces which are "wet" on at least one side are traversed. */

  foreach_face (reduction (min:dtmax)) {
    double hi = h[], hn = h[-1,0];
    if (hi > dry || hn > dry) {

/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
      double dx = Delta/2.;
      double zi = eta[] - hi + Pa2m * pa[];
      double zl = zi - dx*(geta.x[] - gh.x[] + Pa2m *  gpa.x[]); 
      double zn = eta[-1,0] - hn + Pa2m *  pa[-1,0];
      double zr = zn + dx*(geta.x[-1,0] - gh.x[-1,0] + Pa2m *  gpa.x[-1,0]);
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
      kurganov (hm, hp, um, up, Delta, &fh, &fu, &dtmax);
      fv = (fh > 0. ? u.y[-1,0] + dx*gu.y.x[-1,0] : u.y[] - dx*gu.y.x[])*fh;

/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
      if (!(cell.flags & (fghost|active))) {
	hi = coarse(h,0,0);
	zi = coarse(zb,0,0) + Pa2m * coarse(pa,0,0);
      }
      if (!(neighbor(-1,0).flags & (fghost|active))) {
	hn = coarse(h,-1,0);
	zn = coarse(zb,-1,0) + Pa2m * coarse(pa,-1,0);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

/**
#### Flux update */

      Fh.x[]   = fh;
      Fq.x.x[] = fu - sl;
      S.x[]    = fu - sr;
      Fq.y.x[] = fv;
    }
    else // dry
      Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
  }
  
  boundary_normal ({Fh, S, Fq});

/**
#### Updates for evolving quantities

We store the divergence of the fluxes in the update fields. Note that
these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */

  scalar dh = updates[0];
  vector dhu = { updates[1], updates[2] };
  
  foreach() {
    dh[] = (Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/Delta;
    foreach_dimension()
      dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/Delta+( h[]>dry ? ts.x[]: 0. );
      //dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/Delta;
  }

  return dtmax;
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
static void refine_elevation (Point point, scalar h)
{
  // reconstruction of fine cells using elevation (rather than water depth)
  // (default refinement conserves mass but not lake-at-rest)
  if (h[] >= dry) {
    double eta = zb[] + h[];   // water surface elevation  
    struct { double x, y; } g; // gradient of eta
    foreach_dimension()
      g.x = gradient (zb[-1,0] + h[-1,0], eta, zb[1,0] + h[1,0])/4.;
    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta + g.x*child.x + g.y*child.y - zb[]);
  }
  else {

/**
The "dry" case is a bit more complicated. We look in a 3x3
neighborhood of the coarse parent cell and compute a depth-weighted
average of the "wet" surface elevation $\eta$. We need to do this
because we cannot assume a priori that the surrounding wet cells are
necessarily close to e.g. $\eta = 0$. */

    double v = 0., eta = 0.; // water surface elevation
    // 3x3 neighbourhood
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
	if (h[i,j] >= dry) {
	  eta += h[i,j]*(zb[i,j] + h[i,j]);
	  v += h[i,j];
	}
    if (v > 0.)
      eta /= v; // volume-averaged eta of neighbouring wet cells
    else

/**
If none of the surrounding cells is wet, we assume a default sealevel
at zero. 
I have just changed this next line EML 16/1/2014*/

      eta = min( 0., zb[]);

/**
We then reconstruct the water depth in each child using $\eta$ (of the
parent cell i.e. a first-order interpolation in contrast to the wet
case above) and $z_b$ of the child cells. */

    // reconstruct water depth h from eta and zb
    foreach_child()
      h[] = max(0, eta - zb[]);
  }
}

/**
Cell coarsening is simpler. We first compute the depth-weighted
average of $\eta$ over all the children... */

static void coarsen_elevation (Point point, scalar h)
{
  double eta = 0., v = 0.;
  foreach_child()
    if (h[] > dry) {
      eta += h[]*(zb[] + h[]);
      v += h[];
    }

/**
... and use this in combination with $z_b$ (of the coarse cell) to
compute the water depth $h$.  */
    
  if (v > 0.)
    h[] = max(0., eta/v - zb[]);
  else // dry cell
    h[] = 0.;
}

/**
Finally we define a function which will be called by the user to apply
these reconstructions.  */

void conserve_elevation (void)
{
  h.refine  = refine_elevation;
  h.coarsen = coarsen_elevation;
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

#define radiation(ref) (sqrt (G*h[]) - sqrt(G*max((ref) - zb[], 0.)))

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
    Point point = locate (g->x, g->y);
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
  }
}
