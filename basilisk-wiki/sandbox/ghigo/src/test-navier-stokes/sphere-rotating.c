/**
# Sphere rotating in a free-stream flow for $Re=200$ and $\alpha =1$

This test case is the 3D counterpart of the test case
[cylinder-rotating.c](). We study here the flow induced by the
rotation of a sphere in a free-stream flow. Two dimensionless
parameters govern this test case:

* the Reynolds number $Re = u d/\nu$;

* the non-dimensional rotation rate $\alpha = \omega R/u$.

We solve here the Navier-Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h) for $Re=200$ and $\alpha = 1$. */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define uref (1.) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u
#define Re   (200.) // Particle Reynolds number
#define A    (1.) // Adimensional rotation rate

/**
We also define the shape of the domain. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sphere ((x - p.x), (y - p.y), (z - p.z)));
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
  The domain is $16\times 16\times 16$. */

  L0 = 16.;
  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);

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

u.n[left] = dirichlet ((uref));
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
u.r[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = (uref);
uf.n[bottom] = 0;
uf.n[top]    = 0;
uf.n[back]   = 0;
uf.n[front]  = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(d)/(Re)*fm.x[];
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
  We initialize the particle's rotation speed, in the
  counter-clockwise direction. */

  p_w.x = (A)*(uref)/((d)/2.);
  p_w.z = (A)*(uref)/((d)/2.);

  /**
  We initialize the velocity. */

  foreach()
    u.x[] = cs[]*(uref);    
  boundary ((scalar *) {u});  
}

/**
## Embedded boundaries 

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
      assert (ub +
	      p_w.z*(y + b.y*Delta - p_p.y) == 0.);

      ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      (p_w.z*(x + b.x*Delta - p_p.x) - p_w.x*(z + b.z*Delta - p_p.z)) == 0.);

      ub = u.z.boundary[embed] (point, point, u.z, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      p_w.x*(y + b.y*Delta - p_p.y) == 0.);

      // Pressure
      bool neumann;
      double pb;

      pb = p.boundary[embed] (point, point, p, &neumann);
      assert (!neumann);
      assert (pb +
	      (rho[])/(cs[]+ SEPS)*((-(sq (p_w.x) + sq (p_w.z))*(x + b.x*Delta - p_p.x) +
				     (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.x)*n.x +
				    (-(sq (p_w.x) + sq (p_w.z))*(y + b.y*Delta - p_p.y) +
				     (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.y)*n.y +
				    (-(sq (p_w.x) + sq (p_w.z))*(z + b.z*Delta - p_p.z) +
				     (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.z)*n.z
				    ) == 0.);

      // Pressure gradient
      double gb;
      
      gb = g.x.boundary[embed] (point, point, g.x, &dirichlet);
      assert (dirichlet);
      assert (gb - (-(sq (p_w.x) + sq (p_w.z))*(x + b.x*Delta - p_p.x) +
		    (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.x) == 0.);
      gb = g.y.boundary[embed] (point, point, g.y, &dirichlet);
      assert (dirichlet);
      assert (gb - (-(sq (p_w.x) + sq (p_w.z))*(y + b.y*Delta - p_p.y) +
		    (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.y) == 0.);
      gb = g.z.boundary[embed] (point, point, g.z, &dirichlet);
      assert (dirichlet);
      assert (gb - (-(sq (p_w.x) + sq (p_w.z))*(z + b.z*Delta - p_p.z) +
		    (p_w.x*(x + b.x*Delta - p_p.x) + p_w.z*(z + b.z*Delta - p_p.z))*p_w.z) == 0.);      
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

event logfile (i++; t <= 5.*(tref))
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD =  (Fp.x + Fmu.x)/(0.5*sq ((uref))*pi*sq ((d)/2.));
  double CLy = (Fp.y + Fmu.y)/(0.5*sq ((uref))*pi*sq ((d)/2.));
  double CLz = (Fp.z + Fmu.z)/(0.5*sq ((uref))*pi*sq ((d)/2.));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CLy, CLz);
  fflush (stderr);
}

/**
## Results

We plot the time evolution of the drag and lift coefficients
$C_D$, $C_{L,y}$ and , $C_{L,z}$.

~~~gnuplot Drag and lift coeffficient $C_D$ and $C_L$
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 't*u/R'
set ylabel 'C_{D,L}'
set xrange[0:10]
set yrange[-1:2]
plot 'log' u ($2/0.5):14 w l lw 2 lc rgb "black" t "C_D", \
     ''    u ($2/0.5):15 w l lw 2 lc rgb "blue"  t "C_{L,y}", \
     ''    u ($2/0.5):16 w l lw 2 lc rgb "red"   t "C_{L,z}"
~~~
*/
