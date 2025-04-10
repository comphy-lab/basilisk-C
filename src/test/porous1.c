/**
# Stokes flow through a complex porous medium, randomly refined

This is similar to [porous.c]() but tougher since the mesh refinement
is now completely arbitrary (i.e. independent from the solution). This
further tests robustness of the treatment of arbitrary embedded
boundaries, with arbitrary levels of refinement. 

~~~gnuplot Randomly refined mesh and embedded boundary
set size ratio -1
plot 'inter' w l t ''
~~~
*/

#include "embed.h"
#include "navier-stokes/centered.h"

int maxlevel = 5;

/**
The porous medium is the same as in [porous.c](). */

void porous (scalar cs, face vector fs)
{
  int ns = 800;
  double xc[ns], yc[ns], R[ns];
  srand (0);
  for (int i = 0; i < ns; i++)
    xc[i] = 0.5*noise(), yc[i] = 0.5*noise(), R[i] = 0.01 + 0.02*fabs(noise());

  /**
  Once we have defined the random centers and radii, we can compute
  the levelset function $\phi$ representing the embedded boundary. */

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;

    /**
    Since the medium is periodic, we need to take into account all
    the disk images using periodic symmetries. */
      
    for (double xp = -L0; xp <= L0; xp += L0)
      for (double yp = -L0; yp <= L0; yp += L0)
	for (int i = 0; i < ns; i++)
	  phi[] = intersection (phi[], (sq(x + xp - xc[i]) +
					sq(y + yp - yc[i]) - sq(R[i])));
    phi[] = -phi[];
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}

/**
The domain is the periodic unit square centered on the origin. */

int main()
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */
  
  size (1. [0]);
  
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);

  /**
  We turn off the advection term. The choice of the maximum timestep
  and of the tolerance on the Poisson and viscous solves is not
  trivial. This was adjusted by trial and error to minimize (possibly)
  splitting errors and optimize convergence speed. */
  
  stokes = true;
  DT = 2e-5 [0];
#if 1
  TOLERANCE = HUGE;
  NITERMIN = 2;
#else
  TOLERANCE = 1e-3;
  NITERMIN = 2;
#endif
  N = 1 << maxlevel;

  run();
}

scalar un[];

event init (t = 0)
{

  /**
  The mesh is refined randomly along the embedded boundary. */
  
  porous (cs, fs);
  refine (cs[] && cs[] < 1. &&
	  level < 6 && (int)(100.*rand()/(double)RAND_MAX) == 0);
  porous (cs, fs);
  refine (cs[] && cs[] < 1. &&
	  level < 7 && (int)(200.*rand()/(double)RAND_MAX) == 0);
  porous (cs, fs);
  refine (cs[] && cs[] < 1. &&
	  level < 8 && (int)(400.*rand()/(double)RAND_MAX) == 0);
  
  /**
  We define the porous embedded geometry. */

  porous (cs, fs);
  dump();

  FILE * fp = fopen ("inter", "w");
  output_cells (fp);
  output_facets (cs, fp, fs);
  foreach()
    fprintf (fp, "cs %g %g %g\n", x, y, cs[]);
  foreach_face()
    fprintf (fp, "fs %g %g %g\n", x, y, fs.x[]);
  fclose (fp);
  
  /**
  The gravity vector is aligned with the channel and viscosity is
  unity. */
  
   a[] = {1.,0.};
  mu = fm;

  /**
  The boundary condition is zero velocity on the embedded boundary. */

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  
  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.x[];
}

/**
We check for convergence of the solution. */

event logfile (i++; i <= 400)
{
  double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
  fprintf (stderr, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	   maxlevel, i,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
}

/**
~~~gnuplot Convergence history
reset
set xlabel 'Iterations'
set logscale y
set ytics format '%.0e'
set key bottom left
plot '../porous1.ref' u 2:9 w l t '', '' u 2:10 w l t '', \
    '' u 2:11 w l t '', '' u 2:12 w l t '', '' u 2:13 w l t '', \
    'log' u 2:9 w p t 'du', '' u 2:10 w p t 'resp', \
    '' u 2:11 w p t 'resu', '' u 2:12 w p t 'u.x.sum', '' u 2:13 w p t 'p.max'
~~~

## See also

* [Stokes flow through a complex porous medium](porous.c)
*/
