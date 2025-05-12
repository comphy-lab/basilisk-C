/**
# Drifting cylinder in a planar Poiseuille flow

This test case is based on the numerical work of [Feng et al.,
1994](#feng1994). We investigate the settling of a cylinder of
diameter $d$ in a long channel of aspect ration $W/d = 4$ under a
planar Poiseuille flow.

We study the effect of the density ration $r$ of the particle on the
settling dynamics of the particle. Due to the added-mass effect for
density ratios close to 1, we choose to simulate the cases where
$r=[1.2,1.5]$.

We solve here the 2D Navier-Stokes equations and describe the cylinder
using an [embedded boundary](/src/embed.h). */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle-color.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define h    (4.*(d)) // Width of the channel
#define Re   (120.) // Bulk Reynolds number Uh/nu
#define Fr   (43.56) // Froud number gh/U^2
#if DENSITY == 1 // r = 1.5, lagging
#define r    (1.5) // Ratio of solid to fluid density
#else // r = 1.2, lagging
#define r    (1.2) // Ratio of solid to fluid density
#endif // DENSITY
#define uref (1.) // Characteristic speed
#define tref ((d)/(uref)) // Characteristic time
#define grav ((Fr)*sq (uref)/(h)) // Gravity acceleration Fr*U^2/h

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(y,w)     ((y) - (w)) // + over, - under

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	phi[] = intersection (phi[],
			      intersection ((cylinder ((x + xp - p.x),
						       (y + yp - p.y))),
					    intersection (
							  -(wall (y,  (h)/2.)),
							   (wall (y, -(h)/2.)))));
  }
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = (r); // Ratio of solid and fluid density
const double p_v = (p_volume_cylinder ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_cylinder ((d), (r))),
		    (p_moment_inertia_cylinder ((d), (r)))}; // Particle moment of interia
const coord  p_g = {-(grav), 0.}; // Gravity

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
Since both the walls of the channel and the cylinder are described
using the same embedded volume fraction *cs*, we use the color field
*p_col* to color the cylinder only. */

void p_shape_col (scalar c, coord p)
{
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = HUGE;
    for (double xp = -(L0); xp <= (L0); xp += (L0))
      for (double yp = -(L0); yp <= (L0); yp += (L0))
	phi[] = intersection (phi[],
			      (cylinder ((x + xp - p.x),
					 (y + yp - p.y))));
  }
  boundary ({phi});
  fractions (phi, c);
}

/**
Finally, we define the mesh adaptation parameters and vary the maximum
level of refinement. */

#define lmin (8) // Min mesh refinement level (l=8 is 2pt/d)
#define lmax (13) // Max mesh refinement level (l=13 is 64pt/d)
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $128\times 128$. */

  L0 = 128.;
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

We use inlet boundary conditions with a parabolic velocity profile. */

#define profile(y) ((y) <= (h)/2. && (y) >= -(h)/2. ?		      \
		    4.*((y) + (h)/2.)/(h)*(1. - ((y) + (h)/2.)/(h)) : \
		    0.)

u.n[left] = dirichlet ((uref)*(profile (y)));
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

uf.n[left]   = (uref)*(profile (y));

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(h)/(Re)*fm.x[];
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
  As we are computing an equilibrium solution, we remove the Neumann
  pressure boundary condition which is responsible for
  instabilities. */

  for (scalar s in {p, pf}) {
    s.neumann_zero = true;
  }
  
  /**
  We initialize the embedded boundary. */

  /**
  We first define the particle's initial position. */

  p_p.x = (L0)/4.; // Due to gravity and r>1, the particle goes in the opposite direction of the flow
  p_p.y = -(h)/4.;
  
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
  We initialize the velocity. */

  foreach() 
    u.x[] = cs[] >= 1. ? (uref)*(profile (y)) : 0.;
  boundary ((scalar *) {u});
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

event logfile (i++; t < 60.*(tref))
{
  coord Fp, Fmu;
  p_shape_col (p_col, (p_p));
  boundary ((scalar *) {p_col});
  embed_color_force (p, u, mu, p_col, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((p_u.y) + SEPS)*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((p_u.y) + SEPS)*(d));

  double uslip  = ((p_u.x) - (uref)*(profile ((p_p.y))));
  double Reslip = fabs(uslip)*(d)/((uref)*(h)/(Re));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x/(d), p_p.y/(d),
	   p_u.x/(uref), p_u.y/(uref),
	   p_w.x, p_w.y,
	   CD, CL,
	   uslip/(uref), Reslip
	   );
  fflush (stderr);
}

/**
## Results

#### Results for $\rho/\rho_s = 1.2$

We compare our results to those present in fig. 21 of [Feng et al.,
1994](#feng1994).

According to [Feng et al., 1994](#feng1994), heavy particles ($r>1$)
lag in velocity and have a tendency to settle towards the center of
the channel due to an inertial lift and the wall effects. An overshoot
can be observed for large $r$. The stabilized position is close to the
channel center but not exactly $0$ due to the curvature of the
undisturbed velocity profile.

On the contrary, lighter particles lead in velocity and settle closer
to the wall. For higher slip velocities, the effect of the Poiseuille
velocity profile is reduced and the particle behaves as if it is
sedimenting towards the well, until the lubrification become strong
enough to stop the sedimentation.

~~~gnuplot Time evolution of the particle's trajectory
reset
set terminal svg font ",16"
set key bottom left spacing 1.1
set xlabel 'x/d'
set ylabel 'y/d'
set xrange [-64:64]
set yrange [-1.25:0.25]
set grid y
plot 'log' u 14:15 w l lc rgb "blue" t "Basilisk, l=13"  
~~~

~~~gnuplot Time evolution of the slip velocity
set key top right
set ylabel "u_{slip}/u_{m}"
set yrange [-3:3]
plot 'log' u 14:22 w l lc rgb "blue" t "Basilisk, l=13"  
~~~

~~~gnuplot Time evolution of the slip Reynolds number
set key bottom left
set ylabel "Re_{slip}"
set yrange [0:100]
plot 'log' u 14:23 w l lc rgb "blue" t "Basilisk, l=13"  
~~~

## References

~~~bib
@article{feng1994,
  title={Direct Simulation of Initial Value Problems for the Motion of Solid Bodies in a Newtonian Fluid. Part 2. Couette adn Poiseuille Flows.},
  author={Feng, J. and Hu, H. and Joseph, D.},
  journal={Journal of Fluid Mechanics},
  volume={227},
  pages={271--304},
  year={1994}
}
~~~
*/
