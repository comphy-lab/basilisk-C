/**
# Settling cylinder impacting the bottom wall of a channel

This test case is based on the numerical work of [Koblitz et al.,
2017](#Koblitz2017). We investigate the settling in the horizontal
direction of a cylinder until it impacts the far right wall.

This test case is governed by the density ratio $r=\rho_s/\rho$ and
the Reynolds number $Re = \frac{u_{\mathrm{ref}}d}{\nu}$, where
$u_{\mathrm{ref}} = \sqrt{\frac{\pi d}{2}(\frac{\rho_s}{\rho} - 1)g}$
is the "steady" settling velocity.

We solve here the 2D Navier-Stokes equations and describe the cylinder
using an [embedded boundary](/src/embed.h). */

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.25)    // Diameter of the cylinder
#define rr   (1.25)    // Density ratio of the cylinder

#define nu   (0.1)
#define grav (981.)

#define uref (sqrt (M_PI*(d)/2.*((rr) - 1.)*(grav))) // Characteristic speed of the cylinder
#define tref ((d)/(uref)) // Characteristic time

/**
We also define the shape of the domain. */

#define circle(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (circle ((x - p.x), (y - p.y)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = (rr); // Ratio of solid and fluid density
const double p_v = (p_volume_cylinder ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_cylinder ((d), (rr))),
		    (p_moment_inertia_cylinder ((d), (rr)))}; // Particle moment of interia
const coord  p_g = {(grav), 0.};

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (6)  // Min mesh refinement level (l=6 is 3pt/d)
#define lmax (10) // Max mesh refinement level (l=10 is 43pt/d, l=13 is 341pt/d)
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $24d\times 24d$. */

  L0 = 24.*(d);
  size (L0);
  origin (0., -(L0)/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-3*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  
  /**
  We initialize the grid. */

#if !TREE
  N = 1 << (lmax);
#else
  N = 1 << (lmin);
#endif // !TREE
  init_grid (N);
    
  run();
}

/**
## Boundary conditions

We use no-slip boundary conditions. */

u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);

u.n[right] = dirichlet (0);
u.t[right] = dirichlet (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[bottom] = 0;
uf.n[top]    = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (nu)*fm.x[];
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
  for (scalar s in {u, p, pf}) {
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
  As we are computing an equilibrium solution when the particle
  reaches its settling velocity, we remove the Neumann pressure
  boundary condition which is responsible for instabilities. */

#if PNEUMANN // Non-zero Neumann bc for pressure
  for (scalar s in {p, pf}) {
    s.neumann_zero = false;
  }
#else // Zero Neumann bc for pressure
  for (scalar s in {p, pf}) {
    s.neumann_zero = true;
  }
#endif // PNEUMANN

  /**
  We initialize the embedded boundary.

  We first define the particle's initial position. */

  p_p.x = 8.*(d);
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, p_p);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
                        maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, p_p);
}

/**
## Embedded boundaries */


/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
                 maxlevel = (lmax), minlevel = (1));

  /**
  We do not need here to reset the embedded fractions to avoid
  interpolation errors on the geometry as the is already done when
  moving the embedded boundaries. It might be necessary to do this
  however if surface forces are computed around the embedded
  boundaries. */
}
#endif // TREE

/**
## Outputs */

event logfile (i++; t < 40.*(tref))
{
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
           i, t/(tref), dt/(tref),
           mgp.i, mgp.nrelax, mgp.minlevel,
           mgu.i, mgu.nrelax, mgu.minlevel,
           mgp.resb, mgp.resa,
           mgu.resb, mgu.resa,
           p_p.x/(d),    p_p.y/(d),
           p_u.x/(uref), p_u.y/(uref),
           p_w.x*(tref), p_w.y*(tref)
           );
  fflush (stderr);
}

event snapshots (t = {0., 10.*(0.025474), 20.*(0.025474), 30.*(0.025474), 35.*(0.025474)})
{
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  
  view (fov = 7, camera = "front",
	tx = -0.5, ty = 0.,
	bg = {1,1,1},
	width = 600, height = 200);

  draw_vof ("cs", "fs", filled = -1);
  squares ("omega", min = -30, max = 30, map = cool_warm);
  char name1[80];
  sprintf (name1, "vorticity-t-%.0f.png", t/(tref));
  save (name1);
}

/* event animation (i += 10) */
/* { */
/*   scalar omega[]; */
/*   vorticity (u, omega); */
/*   boundary ({omega}); */
  
/*   view (fov = 7, camera = "front", */
/* 	tx = -0.5, ty = 0., */
/* 	bg = {1,1,1}, */
/* 	width = 600, height = 200); */

/*   draw_vof ("cs", "fs", filled = -1); */
/*   squares ("omega", min = -30, max = 30, map = cool_warm); */
/*   save ("vorticity.mp4"); */
/* } */

/**
## Results

![Vorticity](cylinder-impacting-wall/vorticity.mp4)(loop)

![Vorticity at $t/t_{\mathrm{ref}}=0$.](cylinder-impacting-wall/vorticity-t-0.png)

![Vorticity at $t/t_{\mathrm{ref}}=10$.](cylinder-impacting-wall/vorticity-t-10.png)

![Vorticity at $t/t_{\mathrm{ref}}=20$.](cylinder-impacting-wall/vorticity-t-20.png)

![Vorticity at $t/t_{\mathrm{ref}}=30$.](cylinder-impacting-wall/vorticity-t-30.png)

![Vorticity at $t/t_{\mathrm{ref}}=35$.](cylinder-impacting-wall/vorticity-t-35.png)

~~~gnuplot Time evolution of the particles' trajectory
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 'x/d'
set ylabel 'y/d'
set xrange [0:24]
set yrange [-4:4]
plot 'log' u 14:15 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the particles' horizontal position
set key bottom right
set xlabel 't/tref'
set ylabel 'x/d'
set xrange [0:40]
set yrange [8:25]
plot '../data/Koblitz2017/Koblitz2017-fig6a.csv' u 1:(24 - $2) w p ps 0.5 pt 4 lc rgb "black" t "fig. 6a, Koblitz et al. 2017", \
     'log' u 2:14 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the particles' vertical position
set key top right
set ylabel 'y/d'
set yrange [-1:1]
plot 'log' u 2:15 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the particles' horizontal velocity
set ylabel 'u/uref'
set yrange [0:1]
plot '../data/Koblitz2017/Koblitz2017-fig6b.csv' u 1:(-$2) w p ps 0.5 pt 4 lc rgb "black" t "fig. 6b, Koblitz et al. 2017", \
     'log' u 2:16 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the particles' vertical velocity
set ylabel 'v/uref'
set yrange [-0.1:0.1]
plot 'log' u 2:17 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the particles' rotational velocity
set ylabel 'w*tref'
set yrange [-0.1:0.1]
plot 'log' u 2:18 w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the number of cells between the particle and the wall
set ylabel 'cells'
set yrange [0.1:*]
set logscale y

# Reference data
d    = 0.25;
L0   = 24.*d;
lmax = 10;

plot 'log' u 2:((L0 - ($14*d + d/2.))/L0*2**(lmax)) w l lw 2 lc rgb "blue" t "basilisk"
~~~

## References

~~~bib
@article{Koblitz2017,
  title={Direct numerical simulation of particulate flows with an overset grid method},
  author={Koblitz, AR and Lovett, Sean and Nikiforakis, Nikolaos and Henshaw, William D},
  journal={Journal of computational physics},
  volume={343},
  pages={414--431},
  year={2017},
  publisher={Elsevier}
}
*/
