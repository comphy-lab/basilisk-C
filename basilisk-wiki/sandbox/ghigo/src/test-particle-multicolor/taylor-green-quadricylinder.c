/**
# Four buoyant moving cylinders in Taylor--Green vortices

[Taylor--Green
vortices](http://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex)
are one of the few exact non-trivial solutions of the incompressible
Euler equations. In this test case, we use this solution as initial
condition and check whether the numerical scheme can respect the
balance between non-linear advection terms and pressure
gradients. Numerical diffusion will in particular introduce
dissipation. This dissipation can be quantified and is a useful
measure of the accuracy of the numerical scheme.

Furthermore, we add four cylinders who therefore should move with the
fluid. We set the gravity acceleration to 0.

We solve here the Euler equations and add the cylinders using an
[embedded boundary](/src/embed.h). */

#define p_n (4)

#define PARTICLE_PERIOD_X (1)
#define PARTICLE_PERIOD_Y (1)

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-multicolor.h"
#include "view.h"

/**
## Reference solution */

#define d0   (0.1)
#define d1   (0.05)

#define uref (1.) // Reference velocity, uref
#define tref (1.) // Reference time, tref=d/u

/**
We also define the shape of the cylinders. */

#define circle(x,y,r) (sq (x) + sq (y) - sq (r))

double p0_phi (double xc, double yc)
{
  return circle (xc, yc, (d0)/2.);
}
double p1_phi (double xc, double yc)
{
  return circle (xc, yc, (d1)/2.);
}

/**
We finally define the gravity acceleration. */

const coord  p_g = {0., 0.}; // Gravity, zero

/**
## Setup

We define the mesh adaptation parameters. */

#define lmin (6) // Min mesh refinement level (l=7 is 3pt/d)
#define lmax (9) // Max mesh refinement level (l=10 is 25pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is unity, centered on the origin and periodic in all
  directions. */

  L0 = 1.;
  size (L0);
  origin (-L0/2., -L0/2.);

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
  
  N = 1 << (lmax);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions */

/**
## Properties */

/**
## Initial conditions */

event init (i = 0)
{
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
  We initialize the embedded boundary.

  We define the cylinders' shape. */

  pl[0].phi = p0_phi;
  pl[1].phi = p0_phi;
  pl[2].phi = p1_phi;
  pl[3].phi = p1_phi;

  /**
  Next we define the cylinders' physical parameters. */

  // Ratio of solid and fluid density
  pl[0].r = pl[1].r = 2.;
  pl[2].r = pl[3].r = 1.;
  // Particle volume
  pl[0].v = pl[1].v = (p_volume_cylinder ((d0)));
  pl[2].v = pl[3].v = (p_volume_cylinder ((d1)));
  // Particle moment of inertia
  foreach_dimension() {
    pl[0].i.x = pl[1].i.x = (p_moment_inertia_cylinder ((d0), pl[0].r));
    pl[2].i.x = pl[3].i.x = (p_moment_inertia_cylinder ((d1), pl[2].r));
  }
  
  /**
  We then define the cylinders' initial position. */

  for (int i = 0; i < (p_n); i++)
    foreach_dimension() 
      pl[i].c.x = -3./8. + 1./4.*((double) i);

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
  
  /**
  We initialize the cylinders' velocity with the initial Taylor--Green
  solution. */

  for (int i = 0; i < (p_n); i++) {
    pl[i].u.x = - cos(2.*pi*pl[i].c.x)*sin(2.*pi*pl[i].c.y);
    pl[i].u.y =   sin(2.*pi*pl[i].c.x)*cos(2.*pi*pl[i].c.y);
  }
  
  /**
  We also initialize the velocity, pressure and pressure gradient with
  the initial Taylor--Green solution for velocity and pressure. */
  
  foreach() {
    
    u.x[] = cs[] ? - cos(2.*pi*x)*sin(2.*pi*y) : 0.;
    u.y[] = cs[] ?   sin(2.*pi*x)*cos(2.*pi*y) : 0.;
    p[]   = cs[] ? - (cos(4.*pi*x) + cos(4.*pi*y))/4. : 0.;
    pf[]  = p[];
  }
  boundary ((scalar *) {u, p, pf});

  foreach()
    foreach_dimension()
      g.x[] = center_gradient (p);
  boundary ((scalar *) {g});
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

scalar divv[];

event end_timestep (i++)
{
  foreach () {
    divv[] = 0.;
    foreach_dimension()
      divv[] += (uf.x[1] - uf.x[])*pow (Delta, dimension - 1);
    if (cs[] > 0. && cs[] < 1.) {
      coord b, n;
      double area = embed_geometry (point, &b, &n);
      foreach_dimension() {
	bool dirichlet = true;
	double ufb = area*(uf.x.boundary[embed] (point, point, uf.x, &dirichlet));
	assert (dirichlet);
	divv[] += ufb*n.x*pow (Delta, dimension - 1);  
      }
    }	
  }
}

event logfile (i++; t < 2.*(tref))
{
  /**
  We log the time-evolution of the error, the maximum divergence and of the total
  kinetic energy. */
  
  scalar e[], ke[];
  foreach() {
    double u0 = - cos(2.*pi*x)*sin(2.*pi*y);
    double v0 =   sin(2.*pi*x)*cos(2.*pi*y);
    e[]   = cs[] ? norm(u) - sqrt(sq(u0) + sq(v0)) : nodata;
    ke[]  = cs[] ? sq(u.x[]) + sq(u.y[]) : nodata;
  }
  boundary ((scalar *) {e, divv, ke});

  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   n.avg, n.rms, n.max,
	   normf(divv).max, statsf(ke).sum,
	   pl[0].c.x, pl[0].c.y,
	   pl[0].u.x, pl[0].u.y,
	   pl[0].w.x, pl[0].w.y,
	   pl[1].c.x, pl[1].c.y,
	   pl[1].u.x, pl[1].u.y,
	   pl[1].w.x, pl[1].w.y,
	   pl[2].c.x, pl[2].c.y,
	   pl[2].u.x, pl[2].u.y,
	   pl[2].w.x, pl[2].w.y,
	   pl[3].c.x, pl[3].c.y,
	   pl[3].u.x, pl[3].u.y,
	   pl[3].w.x, pl[3].w.y
	   );
  fflush (stderr);
}

event snapshot (i += 10)
{
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  squares ("omega");
  save ("vorticity.mp4");
}

/**
## Results

![Vorticity](taylor-green-quadricylinder/vorticity.mp4)(loop)

We plot the time evolution of the error. We observe small variations
of the velocity.

~~~gnuplot Time evolution of the error
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,1,10
set ytics format "%.0e" 1.e-18,1.e-4,1.e4
set xlabel 't'
set ylabel '||error||'
set xrange [0:2]
set yrange [1.e-8:1.e3]
set logscale y
plot 'log' u 2:4 w l lw 2 lc rgb "black" t '||error||_{1}', \
     ''    u 2:5 w l lw 2 lc rgb "blue"  t '||error||_{2}', \
     ''    u 2:6 w l lw 2 lc rgb "red"   t '||error||_{inf}'
~~~

~~~gnuplot Time evolution of the maximum divergence of the centered velocity field
set ylabel 'div'
plot 'log' u 2:7 w l lw 2 lc rgb "black" notitle
~~~

~~~gnuplot Time evolution of the kinetic energy
set ylabel 'ke'
plot 'log' u 2:8 w l lw 2 lc rgb "black" notitle
~~~

~~~gnuplot Time evolution of the particles' trajectories
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xtics -100,1,100
set ytics -100,1,100
set xlabel 'x'
set ylabel 'y'
set xrange [-1:1]
set yrange [-1:1]
plot 'log' u 9:10  w l lw 2 lc rgb "black"     t 'pl[0]', \
     ''    u 15:16 w l lw 2 lc rgb "blue"      t 'pl[1]', \
     ''    u 21:22 w l lw 2 lc rgb "red"       t 'pl[2]', \
     ''    u 27:28 w l lw 2 lc rgb "sea-green" t 'pl[3]'
~~~

~~~gnuplot Time evolution of the particles' translation velocities
set xlabel 't'
set ylabel 'p_u'
set xrange [0:2]
set yrange [*:*]
plot 'log' u 2:11 w l lw 2 lc rgb "black"     t 'p_u[0].x', \
     ''    u 2:12 w l lw 2 lc rgb "blue"      t 'p_u[0].y', \
     ''    u 2:17 w l lw 2 lc rgb "red"       t 'p_u[1].x', \
     ''    u 2:18 w l lw 2 lc rgb "sea-green" t 'p_u[1].y', \
     ''    u 2:23 w l lw 2 lc rgb "orange"    t 'p_u[2].x', \
     ''    u 2:24 w l lw 2 lc rgb "magenta"   t 'p_u[2].y', \
     ''    u 2:29 w l lw 2 lc rgb "brown"     t 'p_u[3].x', \
     ''    u 2:30 w l lw 2 lc rgb "gray"      t 'p_u[3].y'
~~~

~~~gnuplot Time evolution of the particles' rotational velocities
set ytics -100,0.01,100
set ylabel 'p_w'
plot 'log' u 2:13 w l lw 2 lc rgb "black"     t 'p_w[0].x', \
     ''    u 2:14 w l lw 2 lc rgb "blue"      t 'p_w[0].y', \
     ''    u 2:19 w l lw 2 lc rgb "red"       t 'p_w[1].x', \
     ''    u 2:20 w l lw 2 lc rgb "sea-green" t 'p_w[1].y', \
     ''    u 2:25 w l lw 2 lc rgb "orange"    t 'p_w[2].x', \
     ''    u 2:26 w l lw 2 lc rgb "magenta"   t 'p_w[2].y', \
     ''    u 2:31 w l lw 2 lc rgb "brown"     t 'p_w[3].x', \
     ''    u 2:32 w l lw 2 lc rgb "gray"      t 'p_w[3].y'
~~~
*/
