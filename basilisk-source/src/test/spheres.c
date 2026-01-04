/**
# Stokes flow past a periodic array of spheres

This is the 3D equivalent of [flow past a periodic array of
cylinders](cylinders.c).

We compare the numerical results with the solution given by the
multipole expansion of [Zick and Homsy, 1982](#zick1982). 

Note that we do not use an adaptive mesh since the 3D gaps are much
wider than for the 2D case. */

#include "grid/multigrid3D.h"
#include "embed.h"
#include "navier-stokes/centered.h"

/**
This is adapted from Table 2 of [Zick and Homsy, 1982](#zick1982), where
the first column is the volume fraction $\Phi$ of the spheres and the
second column is the drag coefficient $K$ such that the force exerted
on each sphere in the array is:
$$
F = 6\pi\mu a K U
$$
with $a$ the sphere radius, $\mu$ the dynamic vicosity and $U$ the
average fluid velocity. */

static double zick[7][2] = {
  {0.027,   2.008},
  {0.064,   2.810},
  {0.125,   4.292},
  {0.216,   7.442},
  {0.343,  15.4},
  {0.45,   28.1},
  {0.5236, 42.1}
};

/**
We can vary the maximum level of refinement, *nc* is the index of the
case in the table above, the radius of the cylinder will be computed
using the volume fraction $\Phi$. */

int maxlevel = 5, nc;
double radius;

int main()
{

  /**
  The domain is the periodic unit cube, centered on the origin. */
  
  size (1.);
  origin (-L0/2., -L0/2., -L0/2.);
  foreach_dimension()
    periodic (right);

  /**
  We turn off the advection term. The choice of the maximum timestep
  and of the tolerance on the Poisson and viscous solves is not
  trivial. This was adjusted by trial and error to minimize (possibly)
  splitting errors and optimize convergence speed. */
   
  stokes = true;
  DT = 2e-2;
  TOLERANCE = HUGE;
  NITERMIN = 10;

  /**
  We do the 7 cases computed by Zick & Homsy. The radius is computed
  from the volume fraction. */
  
  for (nc = 0; nc < 7; nc++) {
    maxlevel = 5;
    N = 1 << maxlevel;
    radius = pow (3.*zick[nc][0]/(4.*pi), 1./3.);
    run();
  }
}

/**
We need an extra field to track convergence. */

scalar un[];

event init (t = 0)
{

  /**
  We initialize the embedded geometry. */

  solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(radius));

  /**
  And set acceleration and viscosity to unity. */
  
  a[] = {1.,0.,0.};
  mu = fm;

  /**
  The boundary condition is zero velocity on the embedded boundary. */

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  u.r[embed] = dirichlet(0);
  
  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.x[];
}

/**
We check for a stationary solution. */

event logfile (i++; i <= 500)
{
  double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
  fprintf (fout, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	   maxlevel, i,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
  fflush (fout);
  
  if (i > 1 && du < 1e-3) {

    /**
    $K$ is computed using formula 4.2 of Zick an Homsy, although the
    $1 - \Phi$ factor is a bit mysterious. */
    
    stats s = statsf(u.x);
    double Phi = 4./3.*pi*cube(radius)/cube(L0);
    double V = s.sum/s.volume;
    double dp = 1., mu = 1.;
    double K = dp*sq(L0)/(6.*pi*mu*radius*V*(1. - Phi));
    fprintf (stderr,
	     "%d %g %g %g %g %g\n",
	     maxlevel, radius, Phi, V, K,
	     zick[nc][1]);

    /**
    We stop. */
    
    return 1; /* stop */
  }
}

/**
The drag coefficient closely matches the results of Zick & Homsy.

~~~gnuplot Drag coefficient as a function of volume fraction
set xlabel 'Volume fraction'
set ylabel 'K'
set logscale y
set grid
set key top left
plot 'log' u 3:6 ps 1 lw 2 t 'Zick and Homsy, 1982',	     \
     '' u 3:5 ps 1 pt 6 lw 2 t '5 levels',		     \
     'spheres.6' u 3:5 ps 1 pt 8 lw 2 t '6 levels'
~~~

This can be further quantified by plotting the relative error. Better
than second-order convergence with spatial resolution is obtained.

~~~gnuplot Relative error on the drag coefficient
set ylabel 'Relative error'
plot 'log' u 3:(abs($6-$5)/$5) w lp t '5 levels',	\
     'spheres.6' u 3:(abs($6-$5)/$5) w lp t '6 levels'
~~~

## References

~~~bib
@Article{zick1982,
  Title                    = {{Stokes flow through periodic arrays of spheres}},
  Author                   = {Zick, A.A. and Homsy, G.M.},
  Journal                  = {Journal of Fluid Mechanics},
  Year                     = {1982},
  Number                   = {1},
  Pages                    = {13--26},
  Volume                   = {115}
}

@article{sangani1982,
  title={Slow flow through a periodic array of spheres},
  author={Sangani, AS and Acrivos, A},
  journal={International Journal of Multiphase Flow},
  volume={8},
  number={4},
  pages={343--360},
  year={1982},
  publisher={Elsevier}
}
~~~

## See also

* [Stokes flow through a complex 3D porous medium](/src/examples/porous3D.c)
*/
