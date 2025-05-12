/**
# Incompressible pulsatile flow in a 3D network model */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "distance.h"
#include "view.h"
#include "../myperfs.h"

/**
## Reference solution */

#define l    (10.)  // Length of the network in the x-direction
#define diam (1.)   // Approx. characterstic vessel diameter
#define Re   (10.)  // Reynolds number
#define T0   (0.8)  // Oscillation period
#define uref (1.)   // Maximum inlet velocity
#define tref (min ((diam)/(uref), (T0))) // Reference time

/**
## Importing the geometry from an *stl* file

This function computes the solid and face fractions given a pointer to
an STL file, an adaptation criteria (maximum relative error on
distance) and a minimum and maximum level. */

#define scmin (1.e-3) // Minimun volume and face fraction,
                      // necessary to avoid spurious velocities

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
  and then call the appropriate VOF functions. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
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
We define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 3pt/d)
#define lmax (8) // Max mesh refinement level (l=8 is 25pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $10\times 10\times 10$. */
  
  L0 = (l);
  size (L0);
  origin (0., -L0/2., -L0/2.);
  
  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
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

We use inlet boundary conditions. */

#define wave(t,T) (max (0., sin (2.*M_PI*(t)/(T))))
#define profile(y,z) (sq (diam) - (sq (y) + sq (z)))

u.n[left] = dirichlet (cs[] ? cs[]*(uref)*(wave (t, (T0))) : 0.);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann   (0);
pf[left]  = neumann   (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
u.r[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left] = fs.n[] ? (uref)*(wave (t, (T0))) : 0.;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(diam)/(Re)*fm.x[];
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

  FILE * fp = fopen ("../data/network/network_v2-binary.stl", "r");
  fraction_from_stl (cs, fs, fp, 5e-4, (lmin), (lmax));
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
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-6,(cmax),(cmax),(cmax)},
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

event logfile (i++)
{
  
  stats nx = statsf (u.x), ny = statsf(u.y), nz = statsf (u.z), np = statsf(p);
  fprintf (stderr, "%g %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   Re,
	   i, t/(T0), dt/(T0),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   nx.sum/(nx.volume + SEPS), nx.min, nx.max,
	   ny.sum/(ny.volume + SEPS), ny.min, ny.max,
	   nz.sum/(nz.volume + SEPS), nz.min, nz.max,
	   np.sum/(np.volume + SEPS), np.min, np.max
	   );
  fflush (stderr);
}

/**
## Snapshots */

event snapshots (t = {1.*(T0), 1.25*(T0), 1.5*(T0), 2.*(T0)})
{
  char name[80];
  
  view (fov = 20, camera = "front",
	tx = -0.5, ty = 1.e-12,
	bg = {0.3,0.4,0.6},
	width = 800, height = 800);

  clear();
  cells ();
  sprintf (name, "mesh-tsT0-%.1f.png", t/(T0));
  save (name);
  
  draw_vof ("cs", "fs", lw = 0.5);
  sprintf (name, "vof-tsT0-%.1f.png", t/(T0));
  save (name);

  squares ("u.x", min = -2.1, max = 2.1, map = cool_warm);
  sprintf (name, "ux-tsT0-%.1f.png", t/(T0));
  save (name);

  squares ("u.y", min = -2., max = 2., map = cool_warm);
  sprintf (name, "uy-tsT0-%.1f.png", t/(T0));
  save (name);
  
  squares ("p", min = -200., max = 200.);
  sprintf (name, "p-tsT0-%.1f.png", t/(T0));
  save (name);
}

event animations (i += 10)
{
  view (fov = 20, camera = "front",
	tx = -0.5, ty = 1.e-12,
	bg = {0.3,0.4,0.6},
	width = 800, height = 800);

  clear();
  cells ();
  save ("mesh.mp4");
  
  draw_vof ("cs", "fs", lw = 0.5);
  save ("vof.mp4");

  squares ("u.x", min = -2.1, max = 2.1, map = cool_warm);
  save ("ux.mp4");

  squares ("u.y", min = -2., max = 2., map = cool_warm);
  save ("uy.mp4");

  squares ("p", min = -200., max = 200.);
  save ("p.mp4");
}

/**
## Results

![Time evolution of the mesh](network-pulsed/mesh.mp4)(loop)

![Time evolution of the embedded fractions](network-pulsed/vof.mp4)(loop)

![Time evolution of the velocity $u_x$](network-pulsed/ux.mp4)(loop)

![Time evolution of the pressure $p$](network-pulsed/p.mp4)(loop)

~~~gnuplot Time evolution of the average, minimun and maximum velocity $u_x$
set terminal svg font ",16"
set key top right spacing 1.1
set xtics 0,0.25,10
set xlabel 't/T'
set ylabel 'u_x'
set xrange [0:2]
set yrange [*:*]
plot 'log' u 3:15 w l lc rgb "black" t "Average", \
     ''    u 3:16 w l lc rgb "blue"  t "Minimun", \
     ''    u 3:17 w l lc rgb "red"   t "Maximum"
~~~

~~~gnuplot Time evolution of the average, minimun and maximum velocity $u_y$
set ylabel 'u_y'
plot 'log' u 3:18 w l lc rgb "black" t "Average", \
     ''    u 3:19 w l lc rgb "blue"  t "Minimun", \
     ''    u 3:20 w l lc rgb "red"   t "Maximum"
~~~

~~~gnuplot Time evolution of the average, minimun and maximum velocity $u_z$
set ylabel 'u_z'
plot 'log' u 3:21 w l lc rgb "black" t "Average", \
     ''    u 3:22 w l lc rgb "blue"  t "Minimun", \
     ''    u 3:23 w l lc rgb "red"   t "Maximum"
~~~

~~~gnuplot Time evolution of the average, minimun and maximum pressure $p$
set ylabel 'p'
plot 'log' u 3:24 w l lc rgb "black" t "Average", \
     ''    u 3:25 w l lc rgb "blue"  t "Minimun", \
     ''    u 3:26 w l lc rgb "red"   t "Maximum"
~~~
*/
