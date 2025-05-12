/**
# Sphere accelerating in a quiescent fluid for $Re=300$

This test case is the 3D counterpart of the test case
[cylinder-accelerating.c](). We study here the flow induced by a
uniformly accelerating sphere. Two dimensionless parameters govern
this test case:

* the Reynolds number $Re = u d/\nu$;

* the non-dimensional acceleration rate $\alpha = a d^3/nu^2$.

We solve here the Navier-Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define nu   (1.e-2) // Viscosity
#define Re_s (100.) // Starting particle Reynolds number
#define Re_t (300.) // Terminal particle Reynolds number
#define A    (500.) // Adimensional acceleration rate
#define u_s  ((Re_s)*(nu)/(d)) // Starting velocity
#define u_t  ((Re_t)*(nu)/(d)) // Terminal velocity

#define uref (u_s) // Reference velocity, uref
#define tref (min ((d)/(uref), sqrt ((d)/(A*sq (nu)/cube ((d)))))) // Reference time, tref=min (d/u, sqrt(D/a))

/**
We also define the shape of the domain, and make sure that the
embedded boundaries are compatible with periodic boundaries. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	for (double zp = -(L0); zp <= (L0); zp += (L0))
	  phi[] = intersection (phi[],
				(sphere ((x + xp - p.x),
					 (y + yp - p.y),
					 (z + zp - p.z))));
  }
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
We define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 2pt/d)
#define lmax (8) // Max mesh refinement level (l=8 is 16pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $16\times 16\times 16$ and tri-periodic. */

  L0 = 16.;
  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);

  foreach_dimension()
    periodic (left);
    
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
  We initialize the embedded boundary. */

  /**
  We first define the sphere's initial position. */

  foreach_dimension()
    p_p.x = (L0)/4.;
    
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

  /**
  We initialize the particle's speed and acceleration. */

  foreach_dimension() {
    p_u.x  = (u_s);
    p_au.x = (A)*sq (nu)/cube ((d));
  }
}

/**
## Embedded boundaries 

The cylinder's position is advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  foreach_dimension() {
    p_u.x += (p_au.x)*(dt);
    p_p.x += (p_u.x)*(dt);
  }
}

/**
We verify here that the velocity and pressure gradient boundary
conditions are correctly computed. */

event check (i++)
{
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {

      // Normal pointing from fluid to solid
      coord b, n;
      embed_geometry (point, &b, &n);

      // Velocity
      bool dirichlet;
      double ub;

      ub = u.x.boundary[embed] (point, point, u.x, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      p_u.x == 0.);
      ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      p_u.y == 0.);
      ub = u.z.boundary[embed] (point, point, u.y, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      p_u.z == 0.);

      // Pressure
      bool neumann;
      double pb;

      pb = p.boundary[embed] (point, point, p, &neumann);
      assert (!neumann);
      assert (pb +
	      rho[]/(cs[] + SEPS)*(p_au.x*n.x + p_au.y*n.y+ p_au.z*n.z) == 0.);

      // Pressure gradient
      double gb;
      
      gb = g.x.boundary[embed] (point, point, g.x, &dirichlet);
      assert (dirichlet);
      assert (gb - p_au.x == 0.);
      gb = g.y.boundary[embed] (point, point, g.y, &dirichlet);
      assert (dirichlet);
      assert (gb - p_au.y == 0.);
      gb = g.z.boundary[embed] (point, point, g.z, &dirichlet);
      assert (dirichlet);
      assert (gb - p_au.z == 0.);
    }
  }
}

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax),(cmax)},
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

event logfile (i++; t <= 6.*(tref))
{
  double Re = p_u.x*(d)/(nu);
  
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CDx = (Fp.x + Fmu.x)/(0.5*sq ((uref))*pi*sq ((d)/2.));
  double CDy = (Fp.y + Fmu.y)/(0.5*sq ((uref))*pi*sq ((d)/2.));
  double CDz = (Fp.z + Fmu.z)/(0.5*sq ((uref))*pi*sq ((d)/2.));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   Re,
	   p_p.x/(d), p_p.y/(d), p_p.z/(d),
	   p_u.x/(u_s), p_u.y/(u_s), p_u.z/(u_s),
	   CDx, CDy, CDz);
  fflush (stderr);
}

/**
## Results

We first plot the time evolution of the Reynolds number $Re$

~~~gnuplot Time evolution of the Reynolds number $Re$
set terminal svg font ",16"
set key top right spacing 1.1
set xtics 0,2,100
set xlabel 't*u/d'
set ylabel 'Re'
set xrange[0:6]
plot 'log' u 2:14 w l lw 2 lc rgb "black" notitle
~~~

We then plot the time evolution of the position and velocity of the
cylinder.

~~~gnuplot Time evolution of the position of the cylinder
set ylabel 'p/d'
plot 'log' u 2:15 w l lw 2 lc rgb "black" t 'x', \
     ''    u 2:16 w l lw 2 lc rgb "blue"  t 'y', \
     ''    u 2:17 w l lw 2 lc rgb "red"   t 'z'
~~~

~~~gnuplot Time evolution of the velocity of the cylinder
set key bottom right
set ylabel 'u/u_s'
plot 'log' u 2:18 w l lw 2 lc rgb "black" t 'u_x', \
     ''    u 2:19 w l lw 2 lc rgb "blue"  t 'u_y', \
     ''    u 2:20 w l lw 2 lc rgb "red"   t 'u_z'
~~~

We finally plot the time evolution of the drag and lift coefficients
$C_D$ and $C_L$.

~~~gnuplot Time evolution of the drag coeffficients
set key top right
set ylabel 'C_{D}'
set yrange[-5:5]
plot 'log' u 2:21 w l lw 2 lc rgb "black" t "C_{D,x}", \
     ''    u 2:22 w l lw 2 lc rgb "blue"  t "C_{D,y}", \
     ''    u 2:23 w l lw 2 lc rgb "red"   t "C_{D,z}"
~~~
*/
