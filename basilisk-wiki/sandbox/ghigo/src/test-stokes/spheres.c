/**
# Stokes flow past a periodic array of spheres

This is the 3D equivalent of [flow past a periodic array of
cylinders](cylinders.c).

We compare the numerical results with the solution given by the
multipole expansion of [Zick and Homsy, 1982](#zick1982).

Note that we do not use an adaptive mesh since the 3D gaps are much
wider than for the 2D case. */

#include "grid/multigrid3D.h"
#include "../myembed.h"
#if CENTERED == 1
#include "../mycentered2.h"
#else
#include "../mycentered.h"
#endif // CENTERED

/**
## Reference solution

This is adapted from Table 2 of [Zick and Homsy, 1982](#zick1982),
where the first column is the volume fraction $\Phi$ of the spheres
and the second column is the drag coefficient $K$ such that the force
exerted on each sphere in the array is:
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
We also define the shape of the domain. The radius of the cylinder
will be computed using the volume fraction $\Phi$. */

double radius;

void p_shape (scalar c, face vector f)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq (x) + sq (y) + sq (z) - sq (radius);
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
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
We vary the maximum level of refinement and the index *nc* of the case
in the table above. */

#if LVL
int lmax = (LVL);
#else // lmax = 5
int lmax = 5; // level of refinement (from 12pt/d to 32pt/d)
#endif // LVL
int nc;

int main()
{
  /**
  The domain is $1\times 1\times 1$ and periodic. */
  
  origin (-L0/2., -L0/2., -L0/2.);

  periodic (left);
  periodic (bottom);
  periodic (back);

  /**
  We set the maximum timestep. */

  DT = 2.e-2;

  /**
  We set the tolerance of the Poisson solver. We reduce the tolerance
  of the viscous Poisson solver as the average value of the velocity
  varies between $10^{-1}$ and $10^{-3}$, depending on the volume
  fraction $\Phi$. */

  stokes       = true;
  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-6;

  /**
  We run the 7 cases computed by Zick & Homsy. The radius is computed
  from the volume fraction. */
  
  for (nc = 0; nc < 7; nc++) {

    radius = pow (3.*zick[nc][0]/(4.*pi), 1./3.);

    /**
    We initialize the grid. */

    N = 1 << (lmax);
    init_grid (N);

    run();
  }
}

/**
## Boundary conditions */

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[];
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
  The gravity vector is aligned with the $x$-direction. */
  
  const face vector g[] = {1., 0., 0.};
  a = g;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p})
    s.third = false;
#else
  for (scalar s in {u, p})
    s.third = true;
#endif // ORDER2
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
  cut-cells. */
#endif // TREE
  
  /**
  We initialize the embedded boundary. */
    
  p_shape (cs, fs);
  
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

  /**
  We finally initialize the reference velocity field. */
  
  foreach()
    un[] = u.x[];
}

/**
## Embedded boundaries */

/**
## Outputs 

We check for a stationary solution. */

event logfile (i++; i <= 500)
{
  double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
  fprintf (stdout, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	   lmax, i,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
  fflush (stdout);
  
  if (i > 1 && du < 1e-3) {

    /**
    $K$ is computed using formula 4.2 of Zick an Homsy, although the
    $1 - \Phi$ factor is a bit mysterious. */

    coord Fp, Fmu;
    embed_force (p, u, mu, &Fp, &Fmu);

    stats s = statsf(u.x);
    double Phi = 4./3.*pi*cube(radius)/cube(L0);
    double V = s.sum/s.volume;
    double dp = 1., mu = 1.;
    double K = dp*sq(L0)/(6.*pi*mu*radius*V*(1. - Phi));
    fprintf (stderr,
	     "%d %g %g %g %g %g %g\n",
	     lmax, radius, Phi, V, K,
	     (Fp.x + Fmu.x)/(6.*pi*mu*radius*V*sq (1 - Phi)), // Why sq (1 - phi) ??
	     zick[nc][1]
	     );
    fflush (stderr);

    /**
    We stop. */
    
    return 1; /* stop */
  }
}

/**
## Results

#### Drag coefficient for *l=5*

The drag coefficient closely matches the results of Zick & Homsy.

~~~gnuplot Drag coefficient as a function of volume fraction
set terminal svg font ",16"
set key bottom right spacing 1.1
set grid
set xtics 0,0.1,2
set xlabel 'volume fraction \phi'
set ylabel 'C_D'
set xrange [0:0.6]
set yrange [*:100]
set logscale y
plot 'log' u 3:7 w p pt 7 ps 1.5 lc rgb "black" t 'Zick and Homsy, 1982',	     \
     ''    u 3:5 w p pt 6 ps 1.5 lc rgb "blue"  t 'Basilisk, l=5',		     \
     ''    u 3:6 w p pt 2 ps 1.5 lc rgb "red"   t 'Basilisk, using embed{\_}force, l=5'
~~~

#### Relative error

We now plot the relative error for different levels and observe better
than second-order convergence. Note that this is only true for the
value of the drag evaluated with the norm of the velocity only, not
using *embed_force*.

~~~gnuplot Relative error on the drag coefficient
set key top right
set ylabel 'err'
set yrange [*:1.e0]
plot 'log'             u 3:(abs($7-$5)/$7) w lp pt 7 ps 0.7 lc rgb "black" t 'l=5',	\
     '../spheres1/log' u 3:(abs($7-$5)/$7) w lp pt 5 ps 0.7 lc rgb "blue"  t 'l=6',	\
     '../spheres2/log' u 3:(abs($7-$5)/$7) w lp pt 3 ps 0.7 lc rgb "red"   t 'l=7'
~~~

~~~gnuplot Relative error on the drag coefficient using *embed_force*
plot 'log'             u 3:(abs($7-$6)/$7) w lp pt 7 ps 0.7 lc rgb "black" t 'l=5',	\
     '../spheres1/log' u 3:(abs($7-$6)/$7) w lp pt 5 ps 0.7 lc rgb "blue"  t 'l=6',	\
     '../spheres2/log' u 3:(abs($7-$6)/$7) w lp pt 3 ps 0.7 lc rgb "red"   t 'l=7'
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
