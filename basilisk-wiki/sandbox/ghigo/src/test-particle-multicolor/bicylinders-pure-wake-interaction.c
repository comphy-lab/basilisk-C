/**
# Pure wake interaction of two settling cylinders

This test case is based on the numerical work of [Uhlmann,
2005](#Uhlmann2005), then reproduced by [Koblitz et al.,
2017](#Koblitz2017). We investigate the settling in the horizontal
direction of two cylinders of identical diameter $d$ but of different
density ratios, initialized with a horizontal and vertical offset
where the lighter particle is ahead of the heavier one. The heavier
trailing particle therefore passes the lighter leading particle,
subjecting it to perturbations in its wake.

This test case is governed by the density ratio $r=\rho_s/\rho$ and
the Reynolds number $Re = \frac{u_{\mathrm{ref}}d}{\nu}$, where
$u_{\mathrm{ref}} = \sqrt{\frac{\pi d}{2}(\frac{\rho_s}{\rho} - 1)g}$
is the "steady" settling velocity.

We solve here the 2D Navier-Stokes equations and describe the
cylinders using an [embedded boundary](/src/embed.h). */

#define p_n (4)

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-multicolor.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.2)     // Diameter of the cylinders
#define h    (10.*(d)) // Width of the channel
#define l    (50.*(d)) // Length of the channel
#define r1   (1.5)     // Density ratio of heavy cylinder
#define r2   (1.25)    // Density ratio of light cylinder

#define nu   (0.0008)
#define grav (9.81)

#define uref (sqrt (M_PI*(d)/2.*((r1) - 1.)*(grav))) // Characteristic speed of the heavy cylinder
#define tref ((d)/(uref)) // Characteristic time

/**
We also define the shape of the domain. The first two particles are in
fact the channel, while the other particles are the cylinders. */

#define circle(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(y,w)   ((y) - (w)) // + over, - under

// Channel
double p0_phi (double xc, double yc)
{
  return intersection (-(wall (yc,  (h)/2.)), (wall (yc, -(h)/2.)));
}
double p1_phi (double xc, double yc)
{
  return intersection (-(wall (xc, (l))), (wall (xc, 0.)));
}
// Cylinders
double p2_phi (double xc, double yc)
{
  return circle (xc, yc);
}

/**
We finally define the gravity acceleration. */

const coord  p_g = {(grav), 0.}; // Gravity

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (7)  // Min mesh refinement level (l=7 is 2pt/d)
#define lmax (13) // Max mesh refinement level (l=13 is 152pt/d)
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $50d\times 50d$ with an additional buffer of $2d$ on
  each side. */

  L0 = (l) + 4.*(d);
  size (L0);
  origin (-2.*(d), -L0/2.);

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
## Boundary conditions */

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
  We initialize the embedded boundary.

  We define the shape of the four particles. */

  pl[0].phi = p0_phi;
  pl[1].phi = p1_phi;
  pl[2].phi = p2_phi;
  pl[3].phi = p2_phi;

  /**
  Next we define the particles' physical parameters. The channel is
  defined using the default parameters. */

  // Ratio of solid and fluid density
  pl[2].r = (r1);
  pl[3].r = (r2);
  // Particle volume
  pl[2].v = (p_volume_cylinder ((d)));
  pl[3].v = (p_volume_cylinder ((d)));
  // Particle moment of inertia
  foreach_dimension() {
    pl[2].i.x = (p_moment_inertia_cylinder ((d), pl[2].r));
    pl[3].i.x = (p_moment_inertia_cylinder ((d), pl[3].r));
  }

  /**
  We then define the particles' initial position. */

  // Heavy partilce
  pl[2].c.x = 4.*(d);
  pl[2].c.y = -0.65*(d);
  // Light partilce
  pl[3].c.x = 6.*(d);
  pl[3].c.y = 0.65*(d);
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, pl);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, pl);
}

/**
## Embedded boundaries 

We verify here that the velocity and pressure gradient boundary
conditions for the channel particle are correctly set. */

event check (i++)
{
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {

      // Velocity
      bool dirichlet;
      double ub;

      for (int p_i = 0; p_i <= 1; p_i++) {
	scalar col = pl[p_i].col;
	if (col[] > 0. && col[] < 1.) {
	  ub = u.x.boundary[embed] (point, point, u.x, &dirichlet);
	  assert (dirichlet);
	  assert (ub == 0. && pl[0].u.x == 0.);
	  ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
	  assert (dirichlet);
	  assert (ub == 0. && pl[0].u.y == 0.);
	}
      }

      // Pressure
      bool neumann;
      double pb;
      
      pb = p.boundary[embed] (point, point, p, &neumann);
      assert (!neumann);
      assert (pb == 0.);
      
      // Pressure gradient
      double gb;
      
      gb = g.x.boundary[embed] (point, point, g.x, &dirichlet);
      assert (dirichlet);
      assert (gb == 0.);
      gb = g.y.boundary[embed] (point, point, g.y, &dirichlet);
      assert (dirichlet);
      assert (gb == 0.);
    }
  }
}

/**
## Adaptive mesh refinement */

#if ADAPTOMEGA
scalar omega[];
#endif // ADAPTOMEGA

#if TREE
event adapt (i++)
{
#if ADAPTOMEGA
  vorticity (u, omega);
  boundary ({omega});
  adapt_wavelet ({cs,u,omega}, (double[]) {1.e-2,(cmax),(cmax),1.e-2},
  		 maxlevel = (lmax), minlevel = (1));
#else
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));
#endif // ADAPTOMEGA

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

event logfile (i++; t < 60.*(tref))
{
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   pl[2].c.x/(d),    pl[2].c.y/(d),
	   pl[2].u.x/(uref), pl[2].u.y/(uref),
	   pl[2].w.x*(tref), pl[2].w.y*(tref),
	   pl[3].c.x/(d),    pl[3].c.y/(d),
	   pl[3].u.x/(uref), pl[3].u.y/(uref),
	   pl[3].w.x*(tref), pl[3].w.y*(tref)
	   );
  fflush (stderr);

  /**
  We stop the simulation when the heavy particle is at a distance $2d$
  from the far right boundary. */
  
  if (pl[2].c.x - ((l) - 2.*(d)) >= 0.)
    return 1; // Stop
}

event snapshots (t = {0., 0.8, 3.2, 5.6, 8.})
{
#if !ADAPTOMEGA
  scalar omega[];
#endif // !ADAPTOMEGA
  vorticity (u, omega);
  boundary ({omega});
  
  view (fov = 4, camera = "front",
	tx = -0.5 + 2.*(d)/(L0), ty = 0.,
	bg = {1,1,1},
	width = 1000, height = 200);

  draw_vof ("cs", "fs", filled = -1);
  squares ("omega", min = -30, max = 30, map = cool_warm);
  char name1[80];
  sprintf (name1, "vorticity-t-%.1f.png", t/(tref));
  save (name1);
}

/* event animation (i += 10) */
/* { */
/*   scalar omega[]; */
/*   vorticity (u, omega); */
/*   boundary ({omega}); */
  
/*   view (fov = 4, camera = "front", */
/* 	tx = -0.5 + 2.*(d)/(L0), ty = 0., */
/* 	bg = {1,1,1}, */
/* 	width = 1000, height = 200); */

/*   draw_vof ("cs", "fs", filled = -1); */
/*   squares ("omega", min = -30, max = 30, map = cool_warm); */
/*   save ("vorticity.mp4"); */
/* } */

/**
## Results

![Vorticity](bicylinders-pure-wake-interaction/vorticity.mp4)(loop)

![Vorticity at $t/t_{\mathrm{ref}}=0$.](bicylinders-pure-wake-interaction/vorticity-t-0.0.png)

![Vorticity at $t/t_{\mathrm{ref}}=5$ ($t=0.8$).](bicylinders-pure-wake-interaction/vorticity-t-5.0.png)

![Vorticity at $t/t_{\mathrm{ref}}=19.9$ ($t=3.2$).](bicylinders-pure-wake-interaction/vorticity-t-19.9.png)

![Vorticity at $t/t_{\mathrm{ref}}=34.8$ ($t=5.6$).](bicylinders-pure-wake-interaction/vorticity-t-34.8.png)

![Vorticity at $t/t_{\mathrm{ref}}=49.7$ ($t=8$).](bicylinders-pure-wake-interaction/vorticity-t-49.7.png)

~~~gnuplot Time evolution of the particles' trajectory
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 'x/d'
set ylabel 'y/d'
set xrange [0:50]
set yrange [-5:5]

# Reference values
d = 0.2;
r1 = 1.5;
grav = 9.81;
uref = sqrt (pi*(d)/2.*((r1) - 1.)*(grav));
tref = d/uref;

plot '../data/Uhlmann2005/Uhlmann2005-fig9-heavy-particle.csv' u ($1/d):($2/d) w p ps 0.5 pt 4 lc rgb "black" t "fig. 9, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig9-light-particle.csv' u ($1/d):($2/d) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 9, Uhlmann 2005, light particle", \
     'log' u 14:15 w l lw 2 lc rgb "blue" t "pl[2]",			\
     ''    u 20:21 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

~~~gnuplot Time evolution of the particles' horizontal position
set xlabel 't/tref'
set ylabel 'x/d'
set xrange [0:60]
set yrange [0:100]
plot '../data/Uhlmann2005/Uhlmann2005-fig11a-heavy-particle.csv' u ($1/tref):($2/d) w p ps 0.5 pt 4 lc rgb "black" t "fig. 11a, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig11a-light-particle.csv' u ($1/tref):($2/d) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 11a, Uhlmann 2005, light particle", \
     '../data/Koblitz2017/Koblitz2017-fig13c-heavy-particle.csv' u 1:(-$2) w l lw 1 lc rgb "black" t "fig. 13c, Koblitz et al. 2017, heavy particle", \
     '../data/Koblitz2017/Koblitz2017-fig13c-light-particle.csv' u 1:(-$2) w l lw 1 lc rgb "brown" t "fig. 13c, Koblitz et al. 2017, light particle", \
     'log' u 2:14 w l lw 2 lc rgb "blue" t "pl[2]",			\
     ''    u 2:20 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

~~~gnuplot Time evolution of the particles' vertical position
set ylabel 'y/d'
set yrange [-5:6]
plot '../data/Uhlmann2005/Uhlmann2005-fig12a-heavy-particle.csv' u ($1/tref):($2/d) w p ps 0.5 pt 4 lc rgb "black" t "fig. 12a, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig12a-light-particle.csv' u ($1/tref):($2/d) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 12a, Uhlmann 2005, light particle", \
     '../data/Koblitz2017/Koblitz2017-fig13a-heavy-particle.csv' u 1:($2 - 5) w l lw 1 lc rgb "black" t "fig. 13a, Koblitz et al. 2017, heavy particle", \
     '../data/Koblitz2017/Koblitz2017-fig13a-light-particle.csv' u 1:($2 - 5) w l lw 1 lc rgb "brown" t "fig. 13a, Koblitz et al. 2017, light particle", \
     'log' u 2:15 w l lw 2 lc rgb "blue" t "pl[2]",		 \
     ''    u 2:21 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

~~~gnuplot Time evolution of the particles' horizontal velocity
set ylabel 'u/uref'
set yrange [0:2]
plot '../data/Uhlmann2005/Uhlmann2005-fig11b-heavy-particle.csv' u ($1/tref):($2/uref) w p ps 0.5 pt 4 lc rgb "black" t "fig. 11b, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig11b-light-particle.csv' u ($1/tref):($2/uref) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 11b, Uhlmann 2005, light particle", \
     '../data/Koblitz2017/Koblitz2017-fig13d-heavy-particle.csv' u 1:(-$2) w l lw 1 lc rgb "black" t "fig. 13d, Koblitz et al. 2017, heavy particle", \
     '../data/Koblitz2017/Koblitz2017-fig13d-light-particle.csv' u 1:(-$2) w l lw 1 lc rgb "brown" t "fig. 13d, Koblitz et al. 2017, light particle", \
     'log' u 2:16 w l lw 2 lc rgb "blue" t "pl[2]",		 \
     ''    u 2:22 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

~~~gnuplot Time evolution of the particles' vertical velocity
set ylabel 'v/uref'
set yrange [-0.5:0.8]
plot '../data/Uhlmann2005/Uhlmann2005-fig12b-heavy-particle.csv' u ($1/tref):($2/uref) w p ps 0.5 pt 4 lc rgb "black" t "fig. 12b, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig12b-light-particle.csv' u ($1/tref):($2/uref) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 12b, Uhlmann 2005, light particle", \
     '../data/Koblitz2017/Koblitz2017-fig13b-heavy-particle.csv' u 1:2 w l lw 1 lc rgb "black" t "fig. 13b, Koblitz et al. 2017, heavy particle", \
     '../data/Koblitz2017/Koblitz2017-fig13b-light-particle.csv' u 1:2 w l lw 1 lc rgb "brown" t "fig. 13b, Koblitz et al. 2017, light particle", \
     'log' u 2:17 w l lw 2 lc rgb "blue" t "pl[2]",		 \
     ''    u 2:23 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

~~~gnuplot Time evolution of the particles' rotational velocity
set ylabel 'w*tref'
set yrange [-0.1:0.3]
plot '../data/Uhlmann2005/Uhlmann2005-fig13b-heavy-particle.csv' u ($1/tref):($2*tref) w p ps 0.5 pt 4 lc rgb "black" t "fig. 13b, Uhlmann 2005, heavy particle", \
     '../data/Uhlmann2005/Uhlmann2005-fig13b-light-particle.csv' u ($1/tref):($2*tref) w p ps 0.5 pt 6 lc rgb "brown" t "fig. 13b, Uhlmann 2005, light particle", \
     '../data/Koblitz2017/Koblitz2017-fig13e-heavy-particle.csv' u 1:2 w l lw 1 lc rgb "black" t "fig. 13e, Koblitz et al. 2017, heavy particle", \
     '../data/Koblitz2017/Koblitz2017-fig13e-light-particle.csv' u 1:2 w l lw 1 lc rgb "brown" t "fig. 13e, Koblitz et al. 2017, light particle", \
     'log' u 2:18 w l lw 2 lc rgb "blue" t "pl[2]",		 \
     ''    u 2:24 w l lw 2 lc rgb "red"  t "pl[3]"
~~~

## References

~~~bib
@article{Uhlman2005,
  title={An immersed boundary method with direct forcing for the simulation of particulate flows},
  author={Uhlman, M.},
  journal={Journal of Computational Physics},
  volume={209},
  pages={448--476},
  year={2005}
}

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
