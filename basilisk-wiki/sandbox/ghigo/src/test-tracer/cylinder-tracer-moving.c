/**
# Flow past a moving cylinder at different Reynolds number and Peclet numbers $Pe$

This test case is the moving embedded boundaries counterpart of the
test case [cylinder-tracer.c]().

We also compute the advection-diffusion of a passive scalar for
different Peclet numbers $Pe$. Depending on the test case, we impose
either Dirichlet or Neumann boundary conditions for the passive scalar
on the embedded boundaries.

We solve here the 2D Navier-Stokes equations and add the cylinder
using an [embedded boundary](/src/embed.h). */

#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "../mytracer.h"
#include "../mydiffusion.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define uref (1.) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

#if RE
#define Re ((double) (RE))
#else // Re = 40
#define Re (40.) // Particle Reynolds number Re = ud/nu
#endif // RE

#if PE
#define Pe ((double) (PE))
#else // Pe = 10
#define Pe (10.) // Peclet number Pe = ud/df
#endif // PE

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (cylinder ((x - p.x), (y - p.y)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
## Setup

We define the tracer field *f* and the tracer diffusion coefficient
*df*. */

scalar f[];
face vector df[];
scalar * tracers = {f};

/**
We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
We define the mesh adaptation parameters. */

#define lmin (6) // Min mesh refinement level (l=6 is 2pt/d)
#define lmax (9) // Max mesh refinement level (l=9 is 16pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $32\times 32$. */

  L0 = 32.*(d);
  size (L0);
  origin (-L0/2., -L0/2.);

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

event properties (i++)
{
  foreach_face() {
    muv.x[] = (uref)*(d)/(Re)*fm.x[];
    df.x[]  = (uref)*(d)/(Pe)*fm.x[];
  }
  boundary ((scalar *) {muv, df});
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
  We first define the cylinder's position. */

  p_p.x = 4.*(d);
  
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
  We set the tracer field to 1 on the embedded boundary. */

#if NEUMANN
  f[embed] = neumann (1);
#else
  f[embed] = dirichlet (1);
#endif // NEUMANN
  
  /**
  We initialize the particle's velocity. */

  p_u.x  = -(uref);

  /**
  We finally initialize the reference velocity field. */
  
  foreach()
    un[] = u.x[];
}

/**
## Embedded boundaries 

The cylinder's position is advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  p_p.x -= (uref)*(dt);
}

/**
## Tracer diffusion */

event tracer_diffusion (i++)
{
  diffusion (f, dt, df, theta = cm);
}

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u,f}, (double[]) {1.e-2,(cmax),(cmax),1.e-2},
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

event logfile (i++; t <= 11.*(tref))
{
  double du = change (u.x, un);
  
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(d));

  fprintf (stderr, "%g %g %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g\n",
	   (Re), (Pe),
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CL,
	   statsf(f).min, statsf(f).max,
	   du);
  fflush (stderr);
}

event fields (t = {1.*(tref), 5.*(tref), 10.*(tref)})
{
  /**
  We compute the gradient of *f*. */

  vector gf[];
  foreach()
    foreach_dimension()
      gf.x[] = center_gradient (f);
  boundary ((scalar *) {gf});
  
  double step = ((L0)/(1 << (lmax)));

  // Centerline y=0
  char namex[80];
  sprintf (namex, "fields-x-%.0f-pid-%d.dat", t/(tref), pid());
  FILE * fpx = fopen (namex, "w");

  for (double x = -(L0)/2. + step/2.; x < (L0)/2. - step/2.; x += step)
    fprintf (fpx, "%g %g %g\n", (x - p_p.x)/(d),
	     ((x - p_p.x) <= -(d)/2. || (x - p_p.x) >= (d)/2. ? interpolate (f,    x, p_p.y) : 0.),
	     ((x - p_p.x) <= -(d)/2. || (x - p_p.x) >= (d)/2. ? interpolate (gf.x, x, p_p.y) : 0.));  
  fclose (fpx);

  // Centerline x=0
  char namey[80];
  sprintf (namey, "fields-y-%.0f-pid-%d.dat", t/(tref), pid());
  FILE * fpy = fopen (namey, "w");

  for (double y = -(L0)/2. + step/2.; y < (L0)/2. - step/2.; y += step)
    fprintf (fpy, "%g %g %g\n", (y - p_p.y)/(d),
	     ((y - p_p.y) <= -(d)/2. || (y - p_p.y) >= (d)/2. ? interpolate (f,    p_p.x, y) : 0.),
	     ((y - p_p.y) <= -(d)/2. || (y - p_p.y) >= (d)/2. ? interpolate (gf.y, p_p.x, y) : 0.));  
  fclose (fpy);
}

/**
## Animations */

event movie (i += 10)
{
  scalar omega[];
  vorticity (u, omega);
  boundary ((scalar *) {omega});
  
  view (fov = 4, camera = "front",
	tx = -(p_p.x + (d))/L0, ty = 1.e-12,
	bg = {1,1,1},
	width = 400, height = 200);

  travelling (5, 11, fov = 10, tx = 1.e-12, ty = 1.e-12);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("omega", linear = true, map = cool_warm);
  save ("vorticity.mp4", opt = "-r 12");

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("f", spread = -1, map = cool_warm);
  save ("scalar.mp4", opt = "-r 12");
}

/**
## Results

![Vorticity](cylinder-tracer-moving/vorticity.mp4)(loop)

![Tracer field](cylinder-tracer-moving/scalar.mp4)(loop)

#### Drag coefficient

We plot the time evolution of the drag coefficient $C_D$.

~~~gnuplot Drag coefficient $C_D$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 't/(d/u)'
set ylabel 'C_{D}'
set xrange [0:11]
set yrange [1.2:3]
plot '../cylinder-tracer/log' u 4:16 w l lw 2 lc rgb "blue" t "fixed", \
     'log'                    u 4:16 w l lw 2 lc rgb "red"  t "moving"
~~~

#### Scalar field *f*

We plot the distribution along the centerlines $y=0$ and $x=0$ of the
scalar field *f*.

~~~gnuplot Scalar field *f* along the centerline $y=0$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 'x/d'
set ylabel 'f'
set xrange [-5:15]
set yrange [-0.1:1.1]
set arrow from -0.5, graph 0 to -0.5, graph 1 nohead
set arrow from  0.5, graph 0 to 0.5, graph 1 nohead
plot '< cat ../cylinder-tracer/fields-x-1-pid-*'  u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-x-5-pid-*'  u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-x-10-pid-*' u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat fields-x-1-pid-*'  u 1:2 w l lw 2 lc rgb "blue"      t 't=1', \
     '< cat fields-x-5-pid-*'  u 1:2 w l lw 2 lc rgb "red"       t 't=5', \
     '< cat fields-x-10-pid-*' u 1:2 w l lw 2 lc rgb "sea-green" t 't=10'
~~~

~~~gnuplot Scalar field *f* along the centerline $x=0$
set xlabel 'y/d'
set xrange [-10:10]
plot '< cat ../cylinder-tracer/fields-y-1-pid-*'  u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-y-5-pid-*'  u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-y-10-pid-*' u 1:2 w l lw 2 lc rgb "black" notitle, \
     '< cat fields-y-1-pid-*'  u 1:2 w l lw 2 lc rgb "blue"      t 't=1', \
     '< cat fields-y-5-pid-*'  u 1:2 w l lw 2 lc rgb "red"       t 't=5', \
     '< cat fields-y-10-pid-*' u 1:2 w l lw 2 lc rgb "sea-green" t 't=10'
~~~

#### Gradient of the scalar field *gf*

We plot the distribution along the centerlines $y=0$ and $x=0$ of the
gradient of the scalar field *gf*.

~~~gnuplot Gradient of the scalar field *gf.x* along the centerline $y=0$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 'x/d'
set ylabel 'gf.x'
set xrange [-1:1]
set yrange [-10:10]
set arrow from -0.5, graph 0 to -0.5, graph 1 nohead
set arrow from  0.5, graph 0 to 0.5, graph 1 nohead
plot '< cat ../cylinder-tracer/fields-x-1-pid-*'  u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-x-5-pid-*'  u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-x-10-pid-*' u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat fields-x-1-pid-*'  u 1:3 w l lw 2 lc rgb "blue"      t 't=1', \
     '< cat fields-x-5-pid-*'  u 1:3 w l lw 2 lc rgb "red"       t 't=5', \
     '< cat fields-x-10-pid-*' u 1:3 w l lw 2 lc rgb "sea-green" t 't=10'
~~~

~~~gnuplot Gradient of the scalar field *gf.y* along the centerline $x=0$
set xlabel 'y/d'
set ylabel 'gf.y'
plot '< cat ../cylinder-tracer/fields-y-1-pid-*'  u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-y-5-pid-*'  u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat ../cylinder-tracer/fields-y-10-pid-*' u 1:3 w l lw 2 lc rgb "black" notitle, \
     '< cat fields-y-1-pid-*'  u 1:3 w l lw 2 lc rgb "blue"      t 't=1', \
     '< cat fields-y-5-pid-*'  u 1:3 w l lw 2 lc rgb "red"       t 't=5', \
     '< cat fields-y-10-pid-*' u 1:3 w l lw 2 lc rgb "sea-green" t 't=10'
~~~
*/
