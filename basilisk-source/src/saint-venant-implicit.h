/**
# A semi-implicit Saint-Venant solver

<div class="message">
Note that the [multilayer solver](layered/hydro.h) provides the same
functionality and should be prefered for most applications.</div>

We solve the [Saint-Venant equations](saint-venant.h) semi-implicitly
to lift the timestep restriction due to gravity waves. This is just a
particular application of the "all Mach" semi-implicit solver. We will
use the Bell-Collela-Glaz advection scheme to transport mass and
momentum. */

#include "all-mach.h"
#include "tracer.h"

/**
In addition to the momentum $\mathbf{q}=h\mathbf{u}$ defined by the
all Mach solver, we will need the fluid depth *h* (i.e. the density
$\rho$) and the topography *zb*. Both the momentum $\mathbf{q}$ and
mass $h$ are advected tracers. The acceleration of gravity is
*G*. *dry* is the fluid depth below which a cell is considered
"dry". */

scalar zb[], h[], * tracers = {q,h};
double G = 1.;
double dry = 1e-4;

/**
We need fields to store the (varying) state fields $\rho c^2$,
$\alpha$ and the acceleration $\mathbf{a}$. */

scalar rhoc2v[];
face vector alphav[], av[];

event defaults (i = 0) {
  rho = h;
  alpha = alphav;
  a = av;
  rhoc2 = rhoc2v;

  /**
  We set limiting for *q*, *h*. */
  
  for (scalar s in {q,h})
    s.gradient = minmod2;

  /**
  As well as default values. */

  foreach()
    h[] = 1.;

  /**
  On trees, we ensure that limiting is also applied to prolongation. */
  
  #if TREE
  for (scalar s in {q,h})
    s.prolongation = refine_linear;
  #endif

  /**
  We setup the default display. */

  display ("squares (color = 'h > 0 ? zb + h : nodata', spread = -1);");
}

/**
After user initialisation we set the initial pressure and apply
boundary conditions. The boundary conditions for $\rho=h$ and
$\mathbf{q}$ are already applied by the [all Mach
solver](all-mach.h). */

event init (i = 0) {
  foreach() {
    p[] = G*sq(h[])/2.;
    dimensional (h[] == Delta);
  }
}

/**
By default the timestep is only limited by the CFL condition on the
advection velocity field. We add the option to define an "acoustic CFL
number" which takes into account the speed of gravity waves. */

double CFLa = HUGE;

event stability (i++) {
  if (CFLa < HUGE)
    foreach (reduction (min:dtmax))
      if (h[] > dry) {
	double dt = CFLa*Delta/sqrt(G*h[]);
	if (dt < dtmax)
	  dtmax = dt;
      }
}

/**
The properties of the fluid are given by the "Saint-Venant equation of
state" i.e. $p = gh^2/2$, $\rho c^2 = gh^2$. */

event properties (i++) {
  foreach() {
    rhoc2v[] = G*sq(max(h[],dry));
    ps[] = rhoc2v[]/2.;
  }

  /**
  The specific volume $\alpha=1/\rho$ is constructed by averaging,
  taking into account dry states. */
  
  foreach_face() {
    if ((h[] > dry && h[-1] > dry) ||
	(h[] > dry && h[] + zb[] >= zb[-1]) ||
	(h[-1] > dry && h[-1] + zb[-1] >= zb[]))
      alphav.x[] = 2./(max(h[],dry) + max(h[-1],dry));
    else
      alphav.x[] = 0.;
  }
}

/**
## Topographic source term

The acceleration due to the topography is $- g\nabla z_b$. On regular
Cartesian grids we can simply do */

#if !TREE
event acceleration (i++) {
  foreach_face()
    if (alpha.x[])
      av.x[] -= G*(zb[] - zb[-1])/Delta;
}
#else // TREE

/**
On trees things are a bit more complicated due to the necessity to
verify the "lake-at-rest" condition. For a lake, the pressure gradient
must balance the topographic source term i.e.
$$
\frac{1}{h}\nabla p = a = -g\nabla z_b
$$
with $a$ the acceleration. Using the equation of state and introducing
$$
\eta = z_b + h
$$
the elevation of the free surface (constant for a lake at rest), we get
$$
\frac{1}{2h}\nabla gh^2 = -g\nabla (\eta - h) = g\nabla h
$$
This identity is obviously verified mathematically, however it is not
necessarily verified by the discrete gradient operator. In the case of
Cartesian meshes it is simple to show that the naive discrete gradient
operator we used above verifies this identity. For tree meshes
this is not generally the case due to the prolongation operator used
to fill ghost cells at refinement boundaries. Rather than trying to
redefine the prolongation operator, we discretise the topographic
source term as
$$
a = -g\nabla z_b = -g \nabla\eta + \frac{g}{2h}\nabla h^2
$$
This is strictly identical mathematically, however if we now use the
same discrete operator to estimate $\nabla p$ and $\nabla h^2$, the
*discrete* lake-at-rest condition becomes 
$$ \nabla\eta = 0 $$
which is much easier to verify.

To do so, we define the following restriction and prolongation
operators for $\eta$. The main idea is to avoid using the elevation of
dry cells. Restriction is done by averaging the elevation of wet
cells. */

void eta_restriction (Point point, scalar s)
{
  double sum = 0., n = 0.;
  foreach_child()
    if (h[] > dry) {
      sum += s[];
      n++;
    }
  s[] = n > 0 ? sum/n : nodata;
}

/**
Prolongation uses linear interpolation with a cell-centered gradient
computed from simple finite-differences of wet cells. */

void eta_prolongation (Point point, scalar s)
{
  coord g;
  foreach_dimension() {
    if (h[] > dry) {
      if (h[1] > dry && h[-1] > dry)
	g.x = (s[1] - s[-1])/2.;
      else if (h[1] > dry)
	g.x = s[1] - s[];
      else if (h[-1] > dry)
	g.x = s[] - s[-1];
      else
	g.x = 0.;
    }
    else
      g.x = 0.;
  }
  double sc = s[];
  foreach_child() {
    s[] = sc;
    foreach_dimension()
      s[] += child.x*g.x/4.;
  }
}

/**
One can check that the prolongation values constructed using these
functions verify $\eta = constant$ for a lake-at-rest. */

event acceleration (i++) {
  
  /**
  To compute the acceleration due to topography, we first compute
  $\eta$ and $h^2$. */
  
  scalar eta[], h2[];
  foreach() {
    h2[] = sq(h[]);
    eta[] = zb[] + h[];
  }

  /**
  We then make sure that $\eta$ uses our restriction/prolongation
  functions. $h^2$ uses the same prolongation functions as $p$ by
  default. */

  eta.prolongation = eta_prolongation;
  eta.restriction = eta_restriction;
    
  /**
  We then compute the acceleration as
  $$
  a = -g \nabla\eta + g \frac{\alpha}{2}\nabla h^2  
  $$
  */
  
  foreach_face()
    if (alpha.x[])
      av.x[] += G*(eta[-1] - eta[] + alpha.x[]/2.*(h2[] - h2[-1]))/Delta;
}

#endif // TREE

/**
The utility functions in *elevation.h* need to know which gradient we
are using. */

double (* gradient)  (double, double, double) = minmod2;

#include "elevation.h"
#include "gauges.h"
