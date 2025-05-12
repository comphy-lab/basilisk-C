/**
# Settling sphere in a large container at $Re = 41$

This test case is based on the numerical work of [Uhlman,
2005](#uhlman2005) and the experimental work of [Mordant et al.,
2000](#mordant2000). We investigate the settling of a sphere of
diameter $D$ in a large domain, with parameters corresponding to case
1 in [Mordant et al., 2000](#mordant2000).

This test case is governed by the Reynolds number $Re =
\frac{UD}{\nu}$, where $U$ is the "steady" settling velocity, and the
Stokes number $St= 1/9 Re \rho_p/\rho$.

We solve here the 3D Navier-Stokes equations and describe the sphere
using an [embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "../myperfs.h"
#include "lambda2.h"
#include "view.h"

/**
## Reference solution */

#define d    (2./12.)
#define grav (9.81) // Gravity

#if RE == 1 // Re = 360
#define nu   (0.00104238) // Viscosity
#elif RE == 2 // Re = 280
#define nu   (0.00267626) // Viscosity
#else // Re = 41
#define nu   (0.005416368) // Viscosity
#endif // RE

#define uref (sqrt ((d)*(grav))) // Characteristic speed
#define tref (sqrt ((d)/(grav))) // Characteristic time

/**
We also define the shape of the domain. */

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
We finally define the particle parameters. */

#if RE == 2 // Re = 280
const double p_r = (7.71); // Ratio of solid and fluid density
const coord  p_i = {(p_moment_inertia_sphere ((d), 7.71)),
		    (p_moment_inertia_sphere ((d), 7.71)),
		    (p_moment_inertia_sphere ((d), 7.71))}; // Particle moment of interia
#else // RE = 41 and 360
const double p_r = (2.56); // Ratio of solid and fluid density
const coord  p_i = {(p_moment_inertia_sphere ((d), 2.56)),
		    (p_moment_inertia_sphere ((d), 2.56)),
		    (p_moment_inertia_sphere ((d), 2.56))}; // Particle moment of interia
#endif // RE
const double p_v = (p_volume_sphere ((d))); // Particle volume
const coord  p_g = {0., -(grav), 0.};

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters and vary the maximum level of
refinement. */

#define lmin (9) // Min mesh refinement level (l=9 is 2.5pt/d)
#if LMAX // 11, 12, 13, 14
#define lmax ((int) (LMAX))
#else // 12
#define lmax (12) // Max mesh refinement level (l=12 is 21pt/d)
#endif // LMAX
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $32^3$. It needs to be sufficiently big to allow for
  long settling times without the particle reaching the bottom. */

  L0 = 32.;
  size (L0);
  origin (-L0/2., 0., -L0/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-3*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);

  /**
  We initialize the grid. */
  
  N = 1 << (lmin); // todo: check the influence of coarser grid
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use no-slip boundary conditions. */

u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);

u.n[right] = dirichlet (0);
u.t[right] = dirichlet (0);
u.r[right] = dirichlet (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
u.r[bottom] = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);
u.r[top] = dirichlet (0);

u.n[back] = dirichlet (0);
u.t[back] = dirichlet (0);
u.r[back] = dirichlet (0);

u.n[front] = dirichlet (0);
u.t[front] = dirichlet (0);
u.r[front] = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[bottom] = 0;
uf.n[top]    = 0;
uf.n[back]   = 0;
uf.n[front]  = 0;

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
  If the simulation is not restarted, we define the initial mesh and
  the initial velocity. */
  
  if (!restore (file = "restart")) { // No restart
    
    /**
    We initialize the embedded boundary. */

    /**
    We first define the particle's initial position. */

    p_p.y = (L0) - 5.*(d);
  
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
  If the simulation is restarted, the proper restarting operations are
  performed in [myembed-moving.h](../myembed-moving.h). */
}

/**
## Embedded boundaries */

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
## Restarts and dumps

Every few characteristic time, we also dump the fluid data for
post-processing and restarting purposes. */

#if DUMP
event dump_data (t += 5.*(tref))
{
  // Dump fluid
  char name [80];
  sprintf (name, "dump-level-%d-t-%g", (lmax), t/(tref));
  dump (name);

  // Dump particle
  char p_name [80];
  sprintf (p_name, "p_dump-level-%d-t-%g", (lmax), t/(tref));
  particle pp = {p_p, p_u, p_w, p_au, p_aw};
  struct p_Dump pp_Dump = {p_name, &pp};
  p_dump (pp_Dump);
}
#endif // DUMP

event dump_end (t = end)
{
  // Dump fluid
  char name [80];
  sprintf (name, "dump-level-%d-final", (lmax));
  dump (name);

  // Dump particle
  char p_name [80];
  sprintf (p_name, "p_dump-level-%d-final", (lmax));
  particle pp = {p_p, p_u, p_w, p_au, p_aw};
  struct p_Dump pp_Dump = {p_name, &pp};
  p_dump (pp_Dump);
}

/**
## Profiling */

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs */

event coeffs (i++; t < 75.*(tref))
{
  char name1[80];
  sprintf (name1, "level-%d.dat", lmax);
  static FILE * fp = fopen (name1, "a"); // In case of restart
  fprintf (fp, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x, p_p.y, p_p.z,
	   p_u.x/(uref), p_u.y/(uref), p_u.z/(uref),
	   fabs(p_u.y)*(d)/(nu), 1./9.*fabs(p_u.y)*(d)/(nu)*(p_r)
	   );
  fflush (fp);

  double cell_wall = fabs (p_p.y - (d)/2.)/((L0)/(1 << (lmax)));
  if (cell_wall <= 1.)
    return 1; // stop
}

event logfile (t = end)
{
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x, p_p.y, p_p.z,
	   p_u.x/(uref), p_u.y/(uref), p_u.z/(uref),
	   fabs(p_u.y)*(d)/(nu), 1./9.*fabs(p_u.y)*(d)/(nu)*(p_r)
	   );
  fflush (stderr);
}

/**
## Snapshots
*/

event snapshots (t += 5.*(tref))
{
  scalar l2[];
  lambda2 (u, l2);

  scalar omega[];
  vorticity (u, omega); // Vorticity in xy plane

  char name2[80];

  clear();
  view (fov = 3, //theta = t/(4.*(tref)) - 0.5,
	tx = 0., ty = -(p_p.y + 4.*(d))/(L0),
	bg = {0.3,0.4,0.6},
	width = 800, height = 800);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.1, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em1-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);
  
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em2-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.001, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em3-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);
}

/**
## Animations
*/

event movie (i += 250)
{
  scalar l2[];
  lambda2 (u, l2);

  scalar omega[];
  vorticity (u, omega); // Vorticity in xy plane

  char name2[80];

  clear();
  view (fov = 3, theta = t/(4.*(tref)) - 0.5,
	tx = 0., ty = -(p_p.y + 4.*(d))/(L0),
	bg = {0.3,0.4,0.6},
	width = 800, height = 800);
  
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.1, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em1-level-%d.mp4", lmax);
  save (name2);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em2-level-%d.mp4", lmax);
  save (name2);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.001, color = "omega", map = cool_warm);
  sprintf (name2, "l2-1em3-level-%d.mp4", lmax);
  save (name2);
}

/**
## Results for $Ga = 49.14$

#### $\lambda_2$

![$\lambda_2 = -0.1$ isosurface for l=11](sphere-settling-large-domain/l2-1em1-level-12.mp4)(loop)

![$\lambda_2 = -0.01$ isosurface for l=11](sphere-settling-large-domain/l2-1em2-level-12.mp4)(loop)

![$\lambda_2 = -0.001$ isosurface for l=11](sphere-settling-large-domain/l2-1em3-level-12.mp4)(loop)

#### Settling velocity

We plot the time evolution of the ratio of the computed and analytic
settling velocities. We compare our results to the experimental
results of cases 1,2 and 4 from [Mordant et al., 2000](#mordant2000).

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
set terminal svg font ",16"
set key font ",10" top right spacing 0.7
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:40]
set yrange [-4:1]

#Exponential fit from Mordant et al., 2000
f0(x) = -0.0741/sqrt(0.0005*9.81)*(1. - exp(-3*x*sqrt(0.0005/9.81)/0.055));

plot f0(x) w l lw 2 lc rgb "black" t "Mordant et al., 2000, Ga=49.14", \
     "level-12.dat" u 2:18 w l lw 2 lc rgb "blue"  t "Basilisk, l=12"
~~~

## Results on supercomputer (cedar)

#### Settling velocity for $Ga = 49.14$

We plot the time evolution of the ratio of the computed and analytic
settling velocities. We compare our results to the experimental
results of case 1 from [Mordant et al., 2000](#mordant2000).

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",16"
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:40]
set yrange [-1.5:0]

#Exponential fit from Mordant et al., 2000
f0(x) = -0.0741/sqrt(0.0005*9.81)*(1. - exp(-3*x*sqrt(0.0005/9.81)/0.055));

plot f0(x) w l lw 2 lc rgb "black" t "Mordant et al., 2000, Ga=49.14", \
     "case-0/lmax-12/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-0/lmax-13/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-0/lmax-14/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,x}$
set key font ",16" top right spacing 1.1
set ylabel "u_{p,x}/u_{ref}"
set yrange [-1:1]
plot "case-0/lmax-12/level-12.dat" u 2:17 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-0/lmax-13/level-13.dat" u 2:17 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-0/lmax-14/level-14.dat" u 2:17 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,y}$
set ylabel "u_{p,z}/u_{ref}"
plot "case-0/lmax-12/level-12.dat" u 2:19 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-0/lmax-13/level-13.dat" u 2:19 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-0/lmax-14/level-14.dat" u 2:19 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,y}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
plot "case-0/lmax-12/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-0/lmax-13/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-0/lmax-14/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:*]
set yrange [*:*]
plot "case-0/lmax-12/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-0/lmax-13/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-0/lmax-14/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x"
set ylabel "z"
set zlabel "y"
set xyplane 0
set xrange [-5:2]
set yrange [-5:5]
set zrange [0:*]
splot "case-0/lmax-12/level-12.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
      "case-0/lmax-13/level-13.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
      "case-0/lmax-14/level-14.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

#### Settling velocity for $Ga = 255.35$

We plot the time evolution of the ratio of the computed and analytic
settling velocities. We compare our results to the experimental
results of case 2 from [Mordant et al., 2000](#mordant2000).

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",16"
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:40]
set yrange [-2.5:0]

#Exponential fit from Mordant et al., 2000
f1(x) = -0.218/sqrt(0.0015*9.81)*(1. - exp(-3*x*sqrt(0.0015/9.81)/0.142));

plot f1(x) w l lw 2 lc rgb "black" t "Mordant et al., 2000, Ga=255.35", \
     "case-1/lmax-12/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-1/lmax-13/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-1/lmax-14/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,x}$
set ylabel "u_{p,x}/u_{ref}"
set yrange [-1:1]
plot "case-1/lmax-12/level-12.dat" u 2:17 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-1/lmax-13/level-13.dat" u 2:17 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-1/lmax-14/level-14.dat" u 2:17 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,y}$
set ylabel "u_{p,z}/u_{ref}"
plot "case-1/lmax-12/level-12.dat" u 2:19 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-1/lmax-13/level-13.dat" u 2:19 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-1/lmax-14/level-14.dat" u 2:19 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,y}^2} $
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
plot "case-1/lmax-12/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-1/lmax-13/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-1/lmax-14/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-2:-1.6]
set yrange [*:*]
plot "case-1/lmax-12/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-1/lmax-13/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-1/lmax-14/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x"
set ylabel "z"
set zlabel "y"
set xyplane 0
set xrange [-5:2]
set yrange [-5:5]
set zrange [0:*]
splot "case-1/lmax-12/level-12.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
      "case-1/lmax-13/level-13.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
      "case-1/lmax-14/level-14.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

#### Settling velocity for $Ga = 206.27$

We plot the time evolution of the ratio of the computed and analytic
settling velocities. We compare our results to the experimental
results of case 3 from [Mordant et al., 2000](#mordant2000).

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",16"
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:40]
set yrange [-4.5:0]


#Exponential fit from Mordant et al., 2000
f1(x) = -0.316/sqrt(0.0008*9.81)*(1. - exp(-3*x*sqrt(0.0008/9.81)/0.108));

plot f1(x) w l lw 2 lc rgb "black" t "Mordant et al., 2000, Ga=206.27", \
     "case-2/lmax-12/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-2/lmax-13/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-2/lmax-14/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,x}$
set ylabel "u_{p,x}/u_{ref}"
set yrange [-1:1]
plot "case-2/lmax-12/level-12.dat" u 2:17 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-2/lmax-13/level-13.dat" u 2:17 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-2/lmax-14/level-14.dat" u 2:17 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,y}$
set ylabel "u_{p,z}/u_{ref}"
plot "case-2/lmax-12/level-12.dat" u 2:19 w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-2/lmax-13/level-13.dat" u 2:19 w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-2/lmax-14/level-14.dat" u 2:19 w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,y}^2} $
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
plot "case-2/lmax-12/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-2/lmax-13/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-2/lmax-14/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-4:-3]
set yrange [*:*]
plot "case-2/lmax-12/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
     "case-2/lmax-13/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
     "case-2/lmax-14/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x"
set ylabel "z"
set zlabel "y"
set xyplane 0
set xrange [-5:2]
set yrange [-5:5]
set zrange [0:*]
splot "case-2/lmax-12/level-12.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "blue"      t "Basilisk, l=12", \
      "case-2/lmax-13/level-13.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "red"       t "Basilisk, l=13", \
      "case-2/lmax-14/level-14.dat" u 14:16:(-($15-32)) w l lw 2 lc rgb "sea-green" t "Basilisk, l=14"
~~~

## References

~~~bib
@article{mordant2000,
  title={Velocity measurement of a settling sphere},
  author={Mordant, N. and Pinton, J.-F.},
  journal={The European Physical Journal B},
  volume={18},
  pages={343--352},
  year={2000}
}

@article{uhlman2005,
  title={An immersed boundary method with direct forcing for the simulation of particulate flows},
  author={Uhlman, M.},
  journal={Journal of Computational Physics},
  volume={209},
  pages={448--476},
  year={2005}
}
~~~
*/
