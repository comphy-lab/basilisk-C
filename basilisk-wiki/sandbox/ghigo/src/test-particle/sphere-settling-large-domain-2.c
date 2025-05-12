/**
# Settling sphere in a large container at different $Ga$

We investigate the settling of a sphere of diameter $d$ in a large
domain. This test case is governed by the Galileo number $Ga = \frac{u
d}{\nu}$, where $u = \sqrt{\lvert \frac{\rho_s}{\rho} - 1 \rvert d g}$
is the "steady" settling velocity.

More precisely, we reproduce here 4 test cases taken from [Uhlmann et
al., 2004](#Uhlmann2014) (the article was redacted for an unknow
reason), where the authors vary the Galileo number $Ga$ to reproduce 4
different particle trajectory regimes:

* $\frac{\rho_{p}}{\rho}=1.5,\, Ga = 144$ (case A): steady axisymmetric and vertical trajectory;
* $\frac{\rho_{p}}{\rho}=1.5,\, Ga = 178$ (case B): steady oblique trajectory;
* $\frac{\rho_{p}}{\rho}=1.5,\, Ga = 190$ (case C): low frequency oscillating oblique trajectory;
* $\frac{\rho_{p}}{\rho}=1.5,\, Ga = 250$ (case D): chaotic or intermittent oblique trajectory.

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

#define d       (1.) // Sphere diameter
#define grav    (2.) // Gravity
#define density (1.5)

#if GA // 
#define Ga   ((double) (GA))
#else // Ga = 178
#define Ga   (178) // Particle Reynolds number Re = ud/nu
#endif // GA

#define uref (sqrt (fabs ((density) - 1.)*(d)*(grav))) // Characteristic speed
#define tref ((d)/(uref)) // Characteristic time

/**
We also define the shape of the domain. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sphere ((x - p.x),
		     (y - p.y),
		     (z - p.z)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
We finally define the particle parameters. */

const double p_r = ((density)); // Ratio of solid and fluid density
const coord  p_i = {(p_moment_inertia_sphere ((d), (density))),
		    (p_moment_inertia_sphere ((d), (density))),
		    (p_moment_inertia_sphere ((d), (density)))}; // Particle moment of interia
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

#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d)
#if LMAX // 11, 12, 13, 14
#define lmax ((int) (LMAX))
#else // 12
#define lmax (12) // Max mesh refinement level (l=12 is 16pt/d)
#endif // LMAX
#if CMAX
#define cmax ((((double) (CMAX))*1.e-3)*(uref))
#else // 1.e-2
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field
#endif // CMAX

/**
We also vary the CFL number. */

#if MYCFL
#define mycfl (((double) MYCFL)/100.)
#else
#define mycfl (0.5)
#endif // MYCFL

int main ()
{
  /**
  The domain is $256^3$. It needs to be sufficiently big to allow for
  long settling times without the particle reaching the bottom. */

  L0 = 256.;
  size (L0);
  origin (-L0/2., 0., -L0/2.);

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
    muv.x[] = (uref)*(d)/(Ga)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## CFL condition */

event stability (i++)
{
  if (CFL > (mycfl))
    CFL = (mycfl);
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
  should also define the gradient of *u* at the full center of
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

event coeffs (i++; t < 200.*(tref))
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
	   p_p.x/(d), p_p.y/(d), p_p.z/(d),
	   p_u.x/(uref), p_u.y/(uref), p_u.z/(uref),
	   (Ga), fabs(p_u.y)*(d)/((uref)*(d)/(Ga))
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
	   p_p.x/(d), p_p.y/(d), p_p.z/(d),
	   p_u.x/(uref), p_u.y/(uref), p_u.z/(uref),
	   (Ga), fabs(p_u.y)*(d)/((uref)*(d)/(Ga))
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
  boundary ({l2});

  char name2[80];

  clear();
  view (fov = 7, theta = 0, relative = false,
	tx = -(p_p.x)/(L0), ty = -(p_p.y + 25.*(d))/(L0),
	bg = {1,1,1},
	width = 400, height = 800);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.5, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  sprintf (name2, "l2-5em1-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);
  
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.05, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  sprintf (name2, "l2-5em2-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.005, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  sprintf (name2, "l2-5em3-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);

  clear();
  view (fov = 2, theta = 0, relative = false,
	tx = -(p_p.x)/(L0), ty = -(p_p.y + 7.*(d))/(L0),
	bg = {1,1,1},
	width = 400, height = 800);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.5, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  cells (n = {0, 0, 1}, alpha = (p_p.z));
  sprintf (name2, "l2-zoom-5em1-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);
  
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.05, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  cells (n = {0, 0, 1}, alpha = (p_p.z));
  sprintf (name2, "l2-zoom-5em2-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);

  draw_vof ("cs", "fs");
  isosurface ("l2", -0.005, fc = {0.3,0.4,0.6}, lc = {0.3,0.4,0.6});
  cells (n = {0, 0, 1}, alpha = (p_p.z));
  sprintf (name2, "l2-zoom-5em3-level-%d-t-%g.png", lmax, t/(tref));
  save (name2);
}

/**
## Results on supercomputer (cedar) for *CFL=0.5* and *CMAX=1.e-2* */

/**
#### Settling velocity for $Ga = 100$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.2:-1.0]
set title "Ga=100, CFL=0.5, cmax=1e-2"
plot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-100/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 144$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.35:-1.0]
set title "Case A, Ga=144, CFL=0.5, cmax=1e-2"

# Steady settling values from Uhlamm et al, 2014
uAS = -1.292; # Case A, small domain
uAL = -1.285; # Case A, large domain
uAC15 = -1.2063; # Case AC-15, large domain
uAC36 = -1.2274; # Case AC-36, large domain

plot uAS w l lw 2 lc rgb "black" t "Ga=144 (AS), Uhlmann et al., 2014", \
     uAL w l lw 2 lc rgb "brown" t "Ga=144 (AL), Uhlmann et al., 2014", \
     uAC15 w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     uAC36 w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uAC15 - uAL)/uAL)*100. w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     abs ((uAC36 - uAL)/uAL)*100. w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-144/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 165$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.4:-1.0]
set title "Ga=165, CFL=0.5, cmax=1e-2"
plot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-165/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 178$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case B, Ga=178, CFL=0.5, cmax=1e-2"

# Steady settling values from Uhlamm et al, 2014
uBS = -1.363; # Case B, small domain
uBL = -1.356; # Case B, large domain
uBC15 = -1.2434; # Case BC-15, large domain
uBC48 = -1.3067; # Case BC-48, large domain

plot uBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     uBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     uBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     uBC48 w l lw 2 lc rgb "gray" t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uBC15 - uBL)/uBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((uBC48 - uBL)/uBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:0.5]

# Steady horizontal values from Uhlamm et al, 2014
vBS = 0.1270; # Case B, small domain
vBL = 0.1245; # Case B, large domain
vBC15 = 0.2090; # Case BC-15, large domain
vBC48 = 0.1110; # Case BC-48, large domain

plot vBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     vBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     vBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     vBC48 w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vBC15 - vBL)/vBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((vBC48 - vBL)/vBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]

# Steady angle values from Uhlamm et al, 2014
thetaBS = 5.323; # Case B, small domain
thetaBL = 5.225; # Case B, large domain

plot thetaBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     thetaBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-178/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 190$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case C, Ga=190, CFL=0.5, cmax=1e-2"

# Average settling values from Uhlamm et al, 2014
uCS = -1.383; # Case C, small domain
uCL = -1.376; # Case C, large domain
uCC48 = -1.3233; # Case C-48 

plot uCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     uCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     uCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uCC48 - uCL)/uCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]

# Average horizontal values from Uhlamm et al, 2014
vCS = 0.137; # Case C, small domain
vCL = 0.136; # Case C, large domain
vCC48 = 0.1201; # Case C-48 

plot vCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     vCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     vCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vCC48 - vCL)/vCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-190/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 250$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.75:0]
set title "Case D, Ga=250, CFL=0.5, cmax=1e-2"

# Average settling values from Uhlamm et al, 2014
uDL = -1.4604; # Case D, large domain
uDC36 = -1.3946; # Case DC-36, large domain

plot uDL w l lw 2 lc rgb "brown" t "Ga=250 (DL), Uhlmann et al., 2014", \
     uDC36 w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uDC36 - uDL)/uDL)*100. w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]
plot "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.7:-1.3]
set yrange [*:*]
plot "case-250/CFL-0.50/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.50/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.50/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
## Results on supercomputer (cedar) for *CFL=0.25* and *CMAX=1.e-2* */

/**
#### Settling velocity for $Ga = 100$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.2:-1.0]
set title "Ga=100, CFL=0.25, cmax=1e-2"
plot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-100/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 144$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.35:-1.0]
set title "Case A, Ga=144, CFL=0.25, cmax=1e-2"

# Steady settling values from Uhlamm et al, 2014
uAS = -1.292; # Case A, small domain
uAL = -1.285; # Case A, large domain
uAC15 = -1.2063; # Case AC-15, large domain
uAC36 = -1.2274; # Case AC-36, large domain

plot uAS w l lw 2 lc rgb "black" t "Ga=144 (AS), Uhlmann et al., 2014", \
     uAL w l lw 2 lc rgb "brown" t "Ga=144 (AL), Uhlmann et al., 2014", \
     uAC15 w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     uAC36 w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uAC15 - uAL)/uAL)*100. w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     abs ((uAC36 - uAL)/uAL)*100. w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-144/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 165$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.4:-1.0]
set title "Ga=165, CFL=0.25, cmax=1e-2"
plot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-165/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 178$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case B, Ga=178, CFL=0.25, cmax=1e-2"

# Steady settling values from Uhlamm et al, 2014
uBS = -1.363; # Case B, small domain
uBL = -1.356; # Case B, large domain
uBC15 = -1.2434; # Case BC-15, large domain
uBC48 = -1.3067; # Case BC-48, large domain

plot uBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     uBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     uBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     uBC48 w l lw 2 lc rgb "gray" t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uBC15 - uBL)/uBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((uBC48 - uBL)/uBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:0.5]

# Steady horizontal values from Uhlamm et al, 2014
vBS = 0.1270; # Case B, small domain
vBL = 0.1245; # Case B, large domain
vBC15 = 0.2090; # Case BC-15, large domain
vBC48 = 0.1110; # Case BC-48, large domain

plot vBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     vBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     vBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     vBC48 w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vBC15 - vBL)/vBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((vBC48 - vBL)/vBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]

# Steady angle values from Uhlamm et al, 2014
thetaBS = 5.323; # Case B, small domain
thetaBL = 5.225; # Case B, large domain

plot thetaBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     thetaBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-178/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 190$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case C, Ga=190, CFL=0.25, cmax=1e-2"

# Average settling values from Uhlamm et al, 2014
uCS = -1.383; # Case C, small domain
uCL = -1.376; # Case C, large domain
uCC48 = -1.3233; # Case C-48 

plot uCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     uCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     uCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uCC48 - uCL)/uCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]

# Average horizontal values from Uhlamm et al, 2014
vCS = 0.137; # Case C, small domain
vCL = 0.136; # Case C, large domain
vCC48 = 0.1201; # Case C-48 

plot vCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     vCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     vCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vCC48 - vCL)/vCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-190/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 250$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.75:0]
set title "Case D, Ga=250, CFL=0.25, cmax=1e-2"

# Average settling values from Uhlamm et al, 2014
uDL = -1.4604; # Case D, large domain
uDC36 = -1.3946; # Case DC-36, large domain

plot uDL w l lw 2 lc rgb "brown" t "Ga=250 (DL), Uhlmann et al., 2014", \
     uDC36 w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uDC36 - uDL)/uDL)*100. w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]
plot "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.7:-1.3]
set yrange [*:*]
plot "case-250/CFL-0.25/lmax-12/cmax-10e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-10e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-10e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
## Results on supercomputer (cedar) for *CFL=0.25* and *CMAX=1.e-3* */

/**
#### Settling velocity for $Ga = 100$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.2:-1.0]
set title "Ga=100, CFL=0.25, cmax=1e-3"
plot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-100/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-100/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-100/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 144$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.35:-1.0]
set title "Case A, Ga=144, CFL=0.25, cmax=1e-3"

# Steady settling values from Uhlamm et al, 2014
uAS = -1.292; # Case A, small domain
uAL = -1.285; # Case A, large domain
uAC15 = -1.2063; # Case AC-15, large domain
uAC36 = -1.2274; # Case AC-36, large domain

plot uAS w l lw 2 lc rgb "black" t "Ga=144 (AS), Uhlmann et al., 2014", \
     uAL w l lw 2 lc rgb "brown" t "Ga=144 (AL), Uhlmann et al., 2014", \
     uAC15 w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     uAC36 w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uAC15 - uAL)/uAL)*100. w l lw 2 lc rgb "purple" t "Ga=144 (AC-15), Uhlmann et al., 2014", \
     abs ((uAC36 - uAL)/uAL)*100. w l lw 2 lc rgb "gray"   t "Ga=144 (AC-36), Uhlmann et al., 2014", \
     "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs (($18 - uAL)/uAL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-144/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-144/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-144/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 165$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.4:-1.0]
set title "Ga=165, CFL=0.25, cmax=1e-3"
plot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [*:*]
plot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [*:-1]
set yrange [*:*]
plot "case-165/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-165/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-165/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 178$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case B, Ga=178, CFL=0.25, cmax=1e-3"

# Steady settling values from Uhlamm et al, 2014
uBS = -1.363; # Case B, small domain
uBL = -1.356; # Case B, large domain
uBC15 = -1.2434; # Case BC-15, large domain
uBC48 = -1.3067; # Case BC-48, large domain

plot uBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     uBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     uBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     uBC48 w l lw 2 lc rgb "gray" t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uBC15 - uBL)/uBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((uBC48 - uBL)/uBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs (($18 - uBL)/uBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:0.5]

# Steady horizontal values from Uhlamm et al, 2014
vBS = 0.1270; # Case B, small domain
vBL = 0.1245; # Case B, large domain
vBC15 = 0.2090; # Case BC-15, large domain
vBC48 = 0.1110; # Case BC-48, large domain

plot vBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     vBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     vBC15 w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     vBC48 w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vBC15 - vBL)/vBL)*100. w l lw 2 lc rgb "purple" t "Ga=178 (BC-15), Uhlmann et al., 2014", \
     abs ((vBC48 - vBL)/vBL)*100. w l lw 2 lc rgb "gray"   t "Ga=178 (BC-48), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vBL)/vBL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]

# Steady angle values from Uhlamm et al, 2014
thetaBS = 5.323; # Case B, small domain
thetaBL = 5.225; # Case B, large domain

plot thetaBS w l lw 2 lc rgb "black" t "Ga=178 (BS), Uhlmann et al., 2014", \
     thetaBL w l lw 2 lc rgb "brown" t "Ga=178 (BL), Uhlmann et al., 2014", \
     "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-178/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-178/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-178/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 190$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.5:0]
set title "Case C, Ga=190, CFL=0.25, cmax=1e-3"

# Average settling values from Uhlamm et al, 2014
uCS = -1.383; # Case C, small domain
uCL = -1.376; # Case C, large domain
uCC48 = -1.3233; # Case C-48 

plot uCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     uCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     uCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uCC48 - uCL)/uCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs (($18 - uCL)/uCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]

# Average horizontal values from Uhlamm et al, 2014
vCS = 0.137; # Case C, small domain
vCL = 0.136; # Case C, large domain
vCC48 = 0.1201; # Case C-48 

plot vCS w l lw 2 lc rgb "black" t "Ga=190 (CS), Uhlmann et al., 2014", \
     vCL w l lw 2 lc rgb "brown" t "Ga=190 (CL), Uhlmann et al., 2014", \
     vCC48 w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the horinzontal velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "err(sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref})"
set yrange [0:200]
plot abs ((vCC48 - vCL)/vCL)*100. w l lw 2 lc rgb "purple" t "Ga=190 (CC-48), Uhlmann et al., 2014", \
     "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs ((sqrt($17*$17 + $19*$19) - vCL)/vCL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.45:-1.25]
set yrange [*:*]
plot "case-190/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-190/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-190/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
#### Settling velocity for $Ga = 250$

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
reset
set terminal svg font ",14"
set key font ",14" top right spacing 1.1
set xlabel "t/t_{ref}"
set ylabel "u_{p,y}/u_{ref}"
set xrange [0:200]
set yrange [-1.75:0]
set title "Case D, Ga=250, CFL=0.25, cmax=1e-3"

# Average settling values from Uhlamm et al, 2014
uDL = -1.4604; # Case D, large domain
uDC36 = -1.3946; # Case DC-36, large domain

plot uDL w l lw 2 lc rgb "brown" t "Ga=250 (DL), Uhlmann et al., 2014", \
     uDC36 w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:18 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:18 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:18 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Error for of the settling velocity $u_{p,y}$
set ylabel "err(u_{p,y}/u_{ref})"
set yrange [0:10]
plot abs ((uDC36 - uDL)/uDL)*100. w l lw 2 lc rgb "purple" t "Ga=250 (DC-36), Uhlmann et al., 2014", \
     "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(abs (($18 - uDL)/uDL)*100) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the velocity $\sqrt{u_{p,x}^2 + u_{p,z}^2}$
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set yrange [0:*]
plot "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Angle (in degrees) of the particle motion with respect to the vertical
set ylabel "angle (deg)"
set yrange [0:*]
plot "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:(180./pi*atan (sqrt($17*$17 + $19*$19)/abs($18))) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the time step
set ylabel "dt/t_{ref}"
set yrange [1.e-6:*]
set logscale y
plot "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 2:3 w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 2:3 w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 2:3 w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Time evolution of the particle trajectory
set xlabel "x/d"
set ylabel "z/d"
set zlabel "y/d"
set xyplane 0
set xrange [-5:5]
set yrange [-5:5]
set zrange [0:*]
unset logscale y
splot "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "blue"      t "lmax=12", \
      "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "red"       t "lmax=13", \
      "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 14:16:(-($15-256)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~

~~~gnuplot Phase diagram of the vertical and horizontal velocities
set xtics -2,0.05,-1
set xlabel "u_{p,y}/u_{ref}"
set ylabel "sqrt(u_{p,x}^2 + u_{p,z}^2)/u_{ref}"
set xrange [-1.7:-1.3]
set yrange [*:*]
plot "case-250/CFL-0.25/lmax-12/cmax-1e-3/level-12.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "blue"      t "lmax=12", \
     "case-250/CFL-0.25/lmax-13/cmax-1e-3/level-13.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "red"       t "lmax=13", \
     "case-250/CFL-0.25/lmax-14/cmax-1e-3/level-14.dat" u 18:(sqrt($17*$17 + $19*$19)) w l lw 2 lc rgb "sea-green" t "lmax=14"
~~~
*/

/**
## References

~~~bib
@article{Uhlmann2014,
  title={The motion of a single heavy sphere in ambient fluid: a benchmark for interface-resolved particulate flow simulations with significant relative velocities},
  author={Uhlmann, Markus and Du{\v{s}}ek, Jan},
  journal={International journal of multiphase flow},
  volume={59},
  pages={221--243},
  year={2014},
  publisher={Elsevier}
}
~~~
*/
