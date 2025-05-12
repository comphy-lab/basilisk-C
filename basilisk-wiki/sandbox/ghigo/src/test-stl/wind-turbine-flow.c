/**
# Incompressible flow past a fixed windturbine */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "distance.h"
#include "lambda2.h"
#include "view.h"
#include "../myperfs.h"

/**
## Reference solution */

#define Re   (1000.)      // Reynolds number, u*l/nu
#define l    (100.)       // Characteristic size of the windturbine
#define uref (1.)         // Reference velocity
#define tref ((l)/(uref)) // Reference time, tref=l/u

/**
## Importing the geometry from an *stl* file

This function computes the solid and face fractions given a pointer to
an STL file, an adaptation criteria (maximum relative error on
distance) and a minimum and maximum level. */

#define scmin (1.e-14) // Minimun volume and face fraction,
          
void fraction_from_stl (scalar c, face vector f, FILE * fp, double eps,
			int minlevel, int maxlevel)
{
  /**
  We read the STL file and compute the bounding box of the model. */
  
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  
  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */

  scalar d[];
  distance (d, p);
#if TREE
  while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel, minlevel).nf);
#endif // TREE

  /**
  We also compute the volume fraction from the distance field. We
  first construct a vertex field interpolated from the centered field
  and then call the appropriate VOF functions. We have to be carefull
  with the orientation of the normals (*-* sign here). */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
  		     smin = (scmin), cmin = (scmin));
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
We define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 3pt/l)
#define lmax (9) // Max mesh refinement level (l=9 is 50pt/l)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $1024\times 1024\times 1024$. */
  
  L0 = 1024.;
  size (L0);
  origin (-(L0)/2., -(L0)/2., 0.); // Centered in x and y, z is the altitude
  
  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
  NITERMAX = 200;
  
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);

  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions in the y-direction. */

u.n[bottom] = dirichlet ((uref));
u.t[bottom] = dirichlet (0);
u.r[bottom] = dirichlet (0);
p[bottom]   = neumann (0);
pf[bottom]  = neumann (0);

u.n[top] = neumann (0);
u.t[top] = neumann (0);
u.r[top] = neumann (0);
p[top]   = dirichlet (0);
pf[top]  = dirichlet (0);

/**
The boundary z=0 is the ground with a no-slip boundary condition. */

u.n[back] = dirichlet (0);
u.t[back] = dirichlet (0);
u.r[back] = dirichlet (0);
p[back]   = neumann (0);
pf[back]  = neumann (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[bottom] = (uref);
uf.n[back]   = 0;
uf.n[front]  = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(l)/(Re)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */

#if SLOPELIMITER
  for (scalar s in {u}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
  cut-cells. */
#endif // TREE

  /**
  We initialize the embedded boundary. We read the *stl* file that
  contains the network geometry. */

  FILE * fp = fopen ("../data/wind-turbine/windTurbineGeom.stl", "r");
  fraction_from_stl (cs, fs, fp, 1.e-4, (lmin), max (10, (lmax)));
  fclose (fp);
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. Since we do not have an analytic expression for the
  level-set function, we perform this operation only once. */
  
  adapt_wavelet ({cs}, (double[]) {1.e-30},
		 maxlevel = (lmax), minlevel = (1));
  fractions_cleanup (cs, fs,
  		     smin = (scmin), cmin = (scmin));
#endif // TREE
  
  /**
  We also define the volume fraction at the previous timestep
  *csm1=cs*. */

  csm1 = cs;
    
  /**
  We define the no-slip boundary conditions for the velocity. */
   
  u.n[embed] = dirichlet (0);
  u.t[embed] = dirichlet (0);
  u.r[embed] = dirichlet (0);
  p[embed]   = neumann (0);

  uf.n[embed] = dirichlet (0);
  uf.t[embed] = dirichlet (0);
  uf.r[embed] = dirichlet (0);
  pf[embed]   = neumann (0);

  /**
  We initialize the velocity. */

  foreach()
    u.y[] = cs[]*(uref);
  boundary ((scalar *) {u});
  
  /**
  We finally initialize the reference velocity field. */
  
  foreach()
    un[] = u.y[];
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We do not reset the embedded fractions to avoid interpolation errors
  on the geometry as we do not have an analytic expression for the
  level-set function. */

  fractions_cleanup (cs, fs,
  		     smin = (scmin), cmin = (scmin));
}
#endif // TREE

/**
## Outputs */

event logfile (i++; t <= 100.*(tref))
{
  /**
  We look for a stationary solution. */

  double du = change (u.y, un);

  stats nx = statsf (u.x), ny = statsf(u.y), nz = statsf (u.z), np = statsf(p);
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %ld %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   nx.sum/(nx.volume + SEPS), nx.min, nx.max,
	   ny.sum/(ny.volume + SEPS), ny.min, ny.max,
	   nz.sum/(nz.volume + SEPS), nz.min, nz.max,
	   np.sum/(np.volume + SEPS), np.min, np.max,
	   du,
	   grid->tn, perf.t
	   );
  fflush (stderr);
}

/**
## Snapshots */

event snapshots (t = end)
{
  view (fov = 16,
  	quat = {0.575, -0.191, -0.261, 0.752},
  	ty = -0.2,
  	bg = {0.3,0.4,0.6},
  	width = 800, height = 800);

  clear();
  draw_vof ("cs", "fs", lw = 0.5);
  save ("vof.png");
  
  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.x", n= {1,0,0}, map = cool_warm);
  save ("ux.png");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.y", n= {1,0,0}, map = cool_warm);
  save ("uy.png");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.y", n= {0,1,0}, map = cool_warm);
  save ("uy-2.png");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("p", n= {1,0,0});
  save ("p.png");

  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});

  draw_vof ("cs", "fs", lw = 0.5);
  cells (n= {1,0,0});
  squares ("u.y", n= {1,0,0}, map = cool_warm);
  isosurface ("l2", -0.002);
  save ("l2.png");
}

event animations (i += 10)
{
  view (fov = 16,
  	quat = {0.575, -0.191, -0.261, 0.752},
  	ty = -0.2,
  	bg = {0.3,0.4,0.6},
  	width = 800, height = 800);

  clear();
  cells (n= {1,0,0});
  save ("mesh.mp4");
  
  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.x", n= {1,0,0}, map = cool_warm);
  save ("ux.mp4");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.y", n= {1,0,0}, map = cool_warm);
  save ("uy.mp4");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("u.y", n= {0,1,0}, map = cool_warm);
  save ("uy-2.mp4");

  draw_vof ("cs", "fs", lw = 0.5);
  squares ("p", n= {1,0,0});
  save ("p.mp4");

  scalar l2[];
  lambda2 (u, l2);
  boundary ({l2});

  draw_vof ("cs", "fs", lw = 0.5);
  cells (n= {1,0,0});
  squares ("u.y", n= {1,0,0}, map = cool_warm);
  isosurface ("l2", -0.002);
  save ("l2.mp4");
}

/**
## Results

![$\lambda_2$ criteria](wind-turbine-flow/l2.mp4)(loop)

*/
