/**
# Cylinder rotating in a free-stream flow for $Re=200$ and $\alpha =1$

We reproduce here the results of [Mittal and Kumar,
2003](#mittal2003). The authors studied the flow induced by the
rotation of a cylinder in a free-stream flow. Two dimensionless
parameters govern this test case:

* the Reynolds number $Re = u d/\nu$;

* the non-dimensional rotation rate $\alpha = \omega R/u$.

We solve here the Navier-Stokes equations and add the cylinder using
an [embedded boundary](/src/embed.h) for $Re=200$ and $\alpha = 1$. */

#include "grid/quadtree.h"
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

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (6) // Min mesh refinement level (l=6 is 2pt/d)
#define lmax (10) // Max mesh refinement level (l=10 is 32pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $32\times 32$. */

  L0 = 32.;
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
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions. */

u.n[left] = dirichlet ((uref));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = (uref);
uf.n[bottom] = 0;
uf.n[top]    = 0;

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
  counter-clockwise direction. In 2D, both component of the rotation
  represent $w_z$ and should be identical. */
  
  p_w.x = (A)*(uref)/((d)/2.);
  p_w.y = (A)*(uref)/((d)/2.);

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
	      p_w.x*(y + b.y*Delta - p_p.y) == 0.);
      ub = u.y.boundary[embed] (point, point, u.y, &dirichlet);
      assert (dirichlet);
      assert (ub -
	      p_w.y*(x + b.x*Delta - p_p.x) == 0.);

      // Pressure
      bool neumann;
      double pb;

      pb = p.boundary[embed] (point, point, p, &neumann);
      assert (!neumann);
      assert (pb +
	      (rho[])/(cs[] + SEPS)*(-sq (p_w.x)*(x + b.x*Delta - p_p.x)*n.x -
				     sq (p_w.y)*(y + b.y*Delta - p_p.y)*n.y) == 0.);

      // Pressure gradient
      double gb;
      
      gb = g.x.boundary[embed] (point, point, g.x, &dirichlet);
      assert (dirichlet);
      assert (gb - (-sq (p_w.x)*(x + b.x*Delta - p_p.x)) == 0.);
      gb = g.y.boundary[embed] (point, point, g.y, &dirichlet);
      assert (dirichlet);
      assert (gb - (-sq (p_w.y)*(y + b.y*Delta - p_p.y)) == 0.);
    }
  }
}

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

event logfile (i++; t <= 24.*(tref))
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(d));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CL);
  fflush (stderr);
}

/**
## Snapshots */

event snapshot (t = end)
{
  scalar omega[];
  vorticity (u, omega);
  
  char name2[80];

  /**
  We first plot the whole box. */
 
  view (fov = 20,
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh-level-%d.png", lmax);
  save (name2);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name2, "p-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega-level-%d.png", lmax);
  save (name2);

  /**
  We then zoom on the cylinder. */

  view (fov = 6,
	tx = -(p_p.x + 6.)/(L0), ty = -(p_p.x + 2.)/(L0),
	bg = {1,1,1},
	width = 800, height = 400);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh-zoom-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux-zoom-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy-zoom-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name2, "p-zoom-level-%d.png", lmax);
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega-zoom-level-%d.png", lmax);
  save (name2);
}

/**
## Results

We plot the time evolution of the drag and lift coefficients $C_D$ and
$C_L$. We compare the values with those obtained by [Chen et al.,
1993](#chen1993) and [Mittal and Kumar,2003](#mittal2003).

~~~gnuplot Time evolution of the drag coeffficient $C_D$
set terminal svg font ",16"
set key top right spacing 1.1
set xtics 0,5,50
set ytics -10,1,10
set grid y
set xlabel 't*u/R'
set ylabel 'C_{D}'
set xrange[0:50]
set yrange[0:5]
plot '../data/Chen1993/Chen1993-fig22.csv' u 1:2 w p pt 6 ps 0.75 lc rgb "black" t "Chen et al., 1993, C_D", \
     'log' u ($2/0.5):14 w l lc rgb "blue" t "Basilisk, l=10"
~~~

~~~gnuplot Time evolution of the lift coeffficient $C_L$
set ylabel 'C_{L}'
set yrange[-5:1]
plot 'log' u ($2/0.5):15 w l lc rgb "blue" t "Basilisk, l=10"
~~~

~~~gnuplot Phase diagram of the drag and lift coeffficients
set key bottom right
set xtics -50,0.5,50
set ytics -100,5,100
set xlabel 'C_{D}'
set ylabel 'C_{L}'
set xrange[-2:1.5]
set yrange[-28:2]
plot 'log' u 14:15 w l lc rgb "blue" t "Basilisk, l=10"
~~~

## References

~~~bib
@article{chen1993,
  title={Development of the wake behind a circular cylinder impulsively started into rotatory and rectilinear motion},
  author={Chen, Y.-M and Ou, Y.-R. and Pearlstein, A. J.},
  journal={Journal of Fluid Mechanics},
  volume={253},
  pages={449--484},
  year={1993},
  publisher={Cambridge University Press}
}

@article{mittal2003,
  title={Flow past a rotating cylinder},
  author={Mittal, S. and Kumar, B.},
  journal={Journal of Fluid Mechanics},
  volume={476},
  pages={303--334},
  year={2003},
  publisher={Cambridge University Press}
}
~~~
*/
