/**
# Axisymmetric rising Newtonian bubble in a Newtonian fluid */

#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
// Smering out the interfaces for high density and viscosity ratios.
#define FILTERED 1
#include "mytwo-phase-yield.h"
#include "tension.h"
// Reduced gravity, where gravity is treat as an interfacial force
#define REDUCED 0
#if REDUCED
#include "reduced.h"
#endif // REDUCED
#include "view.h"

/**
## Mesh

We define here the min, max and initial level of refinement along
witht he adaptation criteria for the velocity. */

#define lmin (9) // l=9 is 3pt/D
#define lmax (12) // l=11 is 12pt/D
#define cmax (1.e-2)

/**
## Code */

int main()
{ 
  /**
  The domain is $16\times 16$. */

  L0 = 16.;
  size (L0);
  origin (-L0/2., 0.); // y=r can not be negative

  /**
  We set the numerical parameters: *TOLERANCE* of the poisson and
  viscous solvers, maximum timestep and initial uniform grid size. */

  TOLERANCE = 1.e-4;

  DT = 1.e-2;
  
  N = 1 << (lmin);
  init_grid (N);

  /**
  We set the physical parameters. */

#if REDUCED
  G.x = -10.;
  Z.x = -(L0)/2.;
#endif // REDUCED

  rho1 = 1.;    mu1 = 2.8e-2;
  rho2 = 1.e-3; mu2 = 2.8e-5;
  f.sigma = 8.e-2;

  T1 = 0.;

  run ();
}

/**
## Boundary conditions 

Using (axi.h)[axi.h], the cooridinate system (x,y) corresponds to
(z,r). Therefore, the *bottom* boundary *y=r=0* is the axis of
symmetry and the left and right boundaries *x=z=-L0/2* and *x=z=L0/2*
are the outflow and inflow boundaries. */

//r=0 (default slip)
u.n[bottom]  = dirichlet (0.);
u.t[bottom]  = neumann (0.);

//r=top (default slip)
u.n[top]  = dirichlet (0.);
u.t[top]  = neumann (0.);

//z=left (no-slip)
u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);

//z=right (no-slip)
u.n[right] = dirichlet (0.);
u.t[right] = dirichlet (0.);

/**
We also make sure there is no flow through the top and bottom
boundary, otherwise the compatibility condition for the Poisson
equation can be violated. */

uf.n[bottom] = dirichlet (0.);
uf.n[top]    = dirichlet (0.);

/**
Note that we can only use the following boundary conditions with
reduced gravity. Indeed, when impose the gravity through the
acceleration term, the boundary condition for *uf* should be
proportional to $a*dt$ along the x-direction. */

#if REDUCED
uf.n[left]   = dirichlet (0.);
uf.n[right]  = dirichlet (0.);
#endif // REDUCED

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We define the initial shape of the bubble, characterized by the VOF
  tracer *f* defined in (two-phase.h)[two-phase.h]. We also refine the
  grid locally around the bubble. */

  astats s;
  int ic = 0;
  do {
    ic++;
    fraction (f, sq(x) + sq(y) - sq(0.1));
    boundary ((scalar *) {f});
    s = adapt_wavelet ({f}, (double[]) {1.e-2},
		       maxlevel = (lmax), minlevel = (2));
  } while ((s.nf || s.nc) && ic < 100);

  fraction (f, sq(x) + sq(y) - sq(0.1));
  
  /**
  We set the initial velocity of the fluid. */

  foreach() {
    foreach_dimension()
      u.x[] = 0.;
  }
  boundary ((scalar *) {u});
}

/**
## Properties 

The properties of both phase are updated in the event *properties*
described in (two-phase.h)[two-phase.h].

We also add the acceleration of gravity, acting the *-x=-z*
direction. */

#if !REDUCED
event acceleration (i++) 
{
  face vector av = a;
    foreach_face(x)
      av.x[] += -10.;
}
#endif // !REDUCED

/**
## Mesh adaptation */

event adapt (i++)
{
  adapt_wavelet ((scalar *) {f, u, yuyz}, 
		 (double[]) {1.e-2, (cmax), (cmax), 1.e-2}, 
		 maxlevel = (lmax), minlevel = (2));
}

/**
## End of simulation */

event end_run (t = 0.5) 
{
  return 1;
}

/**
## Post-processing and results */

event logfile (i += 10) 
{
  double vb = 0., bb = 0., ub = 0.;
  foreach(reduction(+:vb) reduction(+:bb) reduction(+:ub)) {
    double dv = (1. - clamp(f[],0,1))*dv()/cm[]; // cm is proportinal to y=r
    vb += dv; // Volume
    bb += x*dv; // Barycenter
    ub += u.x[]*dv; // Volume averaged velocity
  }
  bb /= vb;
  ub /= vb;
    
  fprintf (ferr,
	   "%d %g %g %g %g %g\n", 
	   i, dt, t, vb, bb, ub);
  fflush (ferr);
}

event shape (t = end)
{
  output_facets (f, stdout);  
}

event movie (i += 10)
{
  clear ();

  view (fov = 20,
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 400, height = 400);
    
  draw_vof ("f", lw = 5);
  squares ("u.x", linear = true, map = cool_warm);
  mirror ({0,1}) {
    cells();
    squares (color = "level", min = (lmin), max = (lmax));
  }
  save ("bubble.mp4");

  draw_vof ("f", lw = 5);
  squares ("p", linear = false, map = cool_warm);
  mirror ({0,1}) {
    cells();
    squares (color = "yuyz", min = (0), max = (1));
  }
  save ("yield.mp4");
}

/**
## Plots

![*u.x* velocity field and mesh levels](bubble-axi/bubble.mp4)

![VOF interface and yield criteria (0 yielded, 1 unyielded)](bubble-axi/yield.mp4)

~~~gnuplot Time evolution of the bubble speed
set key top right spacing 1.1
set xlabel "t"
set ylabel "u"
plot "log" u 3:6 w l lw 2 lc rgb "blue" t "Newtonian"
~~~

~~~gnuplot Time evolution of the bubble volume
set ylabel "V"
plot "log" u 3:4 w l lw 2 lc rgb "blue" t "Newtonian"
~~~

~~~gnuplot Bubble shape at the final time
reset
plot 'out' u 1:2 w l lc rgb "blue" t "Newtonian"
~~~
*/
