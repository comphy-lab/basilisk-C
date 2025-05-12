/**
# Settling sphere in a box for $Re = 1.5$

This test case is inspired from [Wachs, 2011](#wachs2011) and the
experimental and numerical work of [ten Cate et al.,
2002](#tenCate2002). We investigae the problem of a settling sphere of
diameter $D$ in a cube box, whith a density ratio $\rho_p / \rho_f
\approx 1$. Note that in the experimental work of [ten Cate et al.,
2002](#tenCate2002), the domain is rectangular.

This test case is governed by the Reynolds number $Re =
\frac{UD}{\nu}$, where $U$ is the "steady" settling velocity, and the
Stokes number $St= \frac{1}{9}Re \frac{\rho_p}{\rho_f}$. The density
ratio $\frac{\rho_p}{\rho_f}$ is given by the Stokes number $St$ and
the viscosity is determined by using the following two expresions for
the drag coefficient: (i) a corrolation from [Abraham,
1970](#abraham1970):
$$
C_d = \frac{24}{(9.06)^2}\left(\frac{9.06}{\sqrt{Re}} +
1\right)^2,
$$
and (ii) a force balance:
$$
C_d = \frac{\rho_f V_p (\frac{\rho_p}{\rho_f} - 1)
g}{\frac{1}{2}\rho_f S_p U^2} = \frac{4}{3}
\frac{D^3}{(\nu Re)^2} (\frac{\rho_p}{\rho_f} - 1)g.
$$

We solve here the Navier-Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-particle.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.015)
#define grav (9.81) // Gravity

#if RE == 1 // Re = 4.1
#define Re   (4.1) // Reynolds number
#define St   (0.53) // Stokes number
#define nu   (2.2130e-04) // Viscosity
#define uref (0.060489) // Reference velocity
#elif RE == 2 // Re = 11.6
#define Re   (11.6) // Reynolds number
#define St   (1.5) // Stokes number
#define nu   (1.1713e-04) // Viscosity
#define uref (0.090579) // Reference velocity
#elif RE == 3 // Re = 31.9
#define Re   (31.9) // Reynolds number
#define St   (4.13) // Stokes number
#define nu   (6.0121e-05) // Viscosity
#define uref (0.12786) // Reference velocity
#else // Re = 1.5
#define Re   (1.5) // Reynolds number
#define St   (0.19) // Stokes number
#define nu   (3.65e-4) // Viscosity
#define uref (0.036500) // Reference velocity
#endif // RE

#define tref ((d)/(uref)) // Reference time, tref

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
We finally define the particle parameters. */

const double p_r = (9.*(St)/(Re)); // Ratio of solid and fluid density
const double p_v = (p_volume_sphere ((d))); // Particle volume
const coord  p_i = {(p_moment_inertia_sphere ((d), 9.*(St)/(Re))),
		    (p_moment_inertia_sphere ((d), 9.*(St)/(Re))),
		    (p_moment_inertia_sphere ((d), 9.*(St)/(Re)))}; // Particle moment of interia
const coord  p_g = {0., -(grav), 0.};

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters and vary the maximum level of
refinement. */

#define lmin (5) // Min mesh refinement level (l=5 is 3pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field
int lmax = 7; // Max mesh refinement level (l=7 is 12pt/d)

int main ()
{  
  /**
  The domain is $0.16^3$. */

  L0 = 0.16;
  size (L0);
  origin (-L0/2., 0., -L0/2.);

  /**
  We set the maximum timestep. Compared to other test cases, we reduce
  here the timestep to avoid oscillations of the settling velocity for
  high levels of refinement. */

  DT = 1.e-3*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);

#if PNEUMANN // Non-zero Neumann bc for pressure
  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
#else // Zero Neumann bc for pressure
  for (lmax = 7; lmax <= 9; lmax++) {
    
    /**
    We initialize the grid. */
  
    N = 1 << (lmin);
    init_grid (N);
    
    run();
  }
#endif // PNEUMANN
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
  We initialize the embedded boundary. */

  /**
  We first define the particle's initial position. */

  p_p.y = 0.12 + (d)/2.;
  
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

event coeffs (i++; t < 5.)
{
  char name1[80];
  sprintf (name1, "level-%d.dat", lmax);
  static FILE * fp = fopen (name1, "w");
  fprintf (fp, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	   i, t, dt,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.y, p_u.y,
	   p_u.x, p_u.z,
	   fabs(p_u.y)*(d)/(nu), 1./9.*fabs(p_u.y)*(d)/(nu)*(p_r)
	   );
  fflush (fp);

  double cell_wall = fabs (p_p.y - (d)/2.)/((L0)/(1 << (lmax)));
  if (cell_wall <= 1.)
    return 1; // stop
}

event logfile (t = end)
{
  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	   i, t, dt,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.y, p_u.y,
	   p_u.x, p_u.z,
	   fabs(p_u.y)*(d)/(nu), 1./9.*fabs(p_u.y)*(d)/(nu)*(p_r)
	   );
  fflush (stderr);
}

/**
## Animations
*/

event snapshots (t = {0.25, 0.5, 0.75, 1.})
{
  scalar omega[];
  vorticity (u, omega); // Vorticity in xy plane
  boundary ({omega});
  
  char name2[80];

  clear();
  view (fov = 22.5796,
	quat = {-0.0470654,0.377745,0.0187008,0.924523},
	tx = 0.00315453, ty = -0.476117, bg = {1.,1.,1.},
	width = 800, height = 800);
  
  draw_vof ("cs", "fs");
  squares ("omega", n = {0, 0, 1}, alpha = 1.e-12, map = cool_warm);
  cells (n = {1, 0, 0}, alpha = 1.e-12);
  sprintf (name2, "vorticity-level-%d-t-%.2f.png", lmax, t);
  save (name2);
}

/**
## Results

#### Settling velocity

We plot the time evolution of the settling velocity. We compare our
results to the experimental results of fig. 5b from [ten Cate et al.,
2002](#tenCate2002).

~~~gnuplot Time evolution of the settling velocity $u_{p,y}$
set terminal svg font ",16"
set key font ",10" bottom right spacing 0.7
set xlabel "t"
set ylabel "u_{p,y}"
set xrange [0:5]
set yrange [-0.14:0.01]
plot "../data/tenCate2002/tenCate2002-fig5b-Re1p5.csv"  u 1:2 w p ps 0.7 pt 7  lc rgb "black" t "fig. 5b, ten Cate et al., 2002, Re=1.5", \
     "../data/tenCate2002/tenCate2002-fig5b-Re4p1.csv"  u 1:2 w p ps 0.7 pt 5  lc rgb "black" t "Re=4.1", \
     "../data/tenCate2002/tenCate2002-fig5b-Re11p6.csv" u 1:2 w p ps 0.7 pt 2  lc rgb "black" t "Re=11.6", \
     "../data/tenCate2002/tenCate2002-fig5b-Re31p9.csv" u 1:2 w p ps 0.7 pt 12 lc rgb "black" t "Re=31.9", \
     "level-7.dat" u 2:15 w l lw 2 lc rgb "blue"      t "Basilisk, l=7", \
     "level-8.dat" u 2:15 w l lw 2 lc rgb "red"       t "Basilisk, l=8", \
     "level-9.dat" u 2:15 w l lw 2 lc rgb "sea-green" t "Basilisk, l=9"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,x}$
set key font ",14" top right spacing 1.1
set ylabel "u_{p,x}"
set yrange [-0.005:0.005]
plot "level-7.dat" u 2:16 w l lw 2 lc rgb "blue"      t "Basilisk, l=7", \
     "level-8.dat" u 2:16 w l lw 2 lc rgb "red"       t "Basilisk, l=8", \
     "level-9.dat" u 2:16 w l lw 2 lc rgb "sea-green" t "Basilisk, l=9"
~~~

~~~gnuplot Time evolution of the velocity $u_{p,z}$
set ylabel "u_{p,z}"
set yrange [-0.005:0.005]
plot "level-7.dat" u 2:16 w l lw 2 lc rgb "blue"      t "Basilisk, l=7", \
     "level-8.dat" u 2:16 w l lw 2 lc rgb "red"       t "Basilisk, l=8", \
     "level-9.dat" u 2:16 w l lw 2 lc rgb "sea-green" t "Basilisk, l=9"
~~~

Next, we plot the time evolution of the sphere's trajectory. We
compare our results to the experimental results of fig. 5a from [ten
Cate et al., 2002](#tenCate2002).

~~~gnuplot Time evolution of the settling trajectory
set ylabel "y/D"
set yrange [0:8]
plot "../data/tenCate2002/tenCate2002-fig5a-Re1p5.csv"  u 1:2 w p ps 0.7 pt 7  lc rgb "black" t "fig. 5a, ten Cate et al., 2002, Re=1.5", \
     "../data/tenCate2002/tenCate2002-fig5a-Re4p1.csv"  u 1:2 w p ps 0.7 pt 5  lc rgb "black" t "Re=4.1", \
     "../data/tenCate2002/tenCate2002-fig5a-Re11p6.csv" u 1:2 w p ps 0.7 pt 2  lc rgb "black" t "Re=11.6", \
     "../data/tenCate2002/tenCate2002-fig5a-Re31p9.csv" u 1:2 w p ps 0.7 pt 12 lc rgb "black" t "Re=31.9", \
     "level-7.dat" u 2:(($14)/0.015 - 0.5) w l lw 2 lc rgb "blue"      t "Basilisk, l=7", \
     "level-8.dat" u 2:(($14)/0.015 - 0.5) w l lw 2 lc rgb "red"       t "Basilisk, l=8", \
     "level-9.dat" u 2:(($14)/0.015 - 0.5) w l lw 2 lc rgb "sea-green" t "Basilisk, l=9"
~~~

Finally, we plot the time evolution of the number of mesh cells
between the bottom of the sphere and the wall.

~~~gnuplot Time evolution of the number of cells from the wall
set grid
set ylabel "N"
set yrange [0.5:100]
set logscale y

D = 0.015;
Dx = 0.16/(2**7);

plot "level-7.dat" u ($2):(($14 - D/2)/Dx) w l lw 2 lc rgb "blue"      t "Basilisk, l=7", \
     "level-8.dat" u ($2):(($14 - D/2)/Dx) w l lw 2 lc rgb "red"       t "Basilisk, l=8", \
     "level-9.dat" u ($2):(($14 - D/2)/Dx) w l lw 2 lc rgb "sea-green" t "Basilisk, l=9"
~~~

## Results for different $Re$

~~~gnuplot Time evolution of the settling velocity
reset
set terminal svg font ",16"
set key bottom right spacing 1.1
set xlabel "t"
set ylabel "u_{p,y}"
set xrange [0:5]
set yrange [-0.14:0.01]
plot "../sphere-settling/case-0/level-7.dat" u 2:15 w l lw 1.25 lc rgb "blue"      t "Basilisk, l=7", \
     "../sphere-settling/case-0/level-8.dat" u 2:15 w l lw 1.25 lc rgb "red"       t "l=8", \
     "../sphere-settling/case-0/level-9.dat" u 2:15 w l lw 1.25 lc rgb "sea-green" t "l=9", \
     "../sphere-settling/case-1/level-7.dat" u 2:15 w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-1/level-8.dat" u 2:15 w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-1/level-9.dat" u 2:15 w l lw 1.25 lc rgb "sea-green" notitle, \
     "../sphere-settling/case-2/level-7.dat" u 2:15 w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-2/level-8.dat" u 2:15 w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-2/level-9.dat" u 2:15 w l lw 1.25 lc rgb "sea-green" notitle, \
     "../sphere-settling/case-3/level-7.dat" u 2:15 w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-3/level-8.dat" u 2:15 w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-3/level-9.dat" u 2:15 w l lw 1.25 lc rgb "sea-green" notitle, \
     "../data/tenCate2002/tenCate2002-fig5b-Re1p5.csv" u 1:2 w p ps 0.7 pt 7 lc rgb "black"   t "Re=1.5", \
     "../data/tenCate2002/tenCate2002-fig5b-Re4p1.csv" u 1:2 w p ps 0.7 pt 5 lc rgb "black"   t "Re=4.1", \
     "../data/tenCate2002/tenCate2002-fig5b-Re11p6.csv" u 1:2 w p ps 0.7 pt 2 lc rgb "black"  t "Re=11.6", \
     "../data/tenCate2002/tenCate2002-fig5b-Re31p9.csv" u 1:2 w p ps 0.7 pt 12 lc rgb "black" t "fig. 5b, ten Cate et al., 2002, Re=31.9"
~~~

~~~gnuplot Time evolution of the settling trajectory
set key top right
set ylabel "y/D"
set yrange [0:8]
plot "../data/tenCate2002/tenCate2002-fig5a-Re1p5.csv"  u 1:2 w p ps 0.7 pt 7  lc rgb "black" t "fig. 5a, ten Cate et al., 2002, Re=1.5", \
     "../data/tenCate2002/tenCate2002-fig5a-Re4p1.csv"  u 1:2 w p ps 0.7 pt 5  lc rgb "black" t "Re=4.1", \
     "../data/tenCate2002/tenCate2002-fig5a-Re11p6.csv" u 1:2 w p ps 0.7 pt 2  lc rgb "black" t "Re=11.6", \
     "../data/tenCate2002/tenCate2002-fig5a-Re31p9.csv" u 1:2 w p ps 0.7 pt 12 lc rgb "black" t "Re=31.9", \
     "../sphere-settling/case-0/level-7.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "blue"      t "Basilisk, l=7", \
     "../sphere-settling/case-0/level-8.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "red"       t "l=8", \
     "../sphere-settling/case-0/level-9.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "sea-green" t "l=9", \
     "../sphere-settling/case-1/level-7.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-1/level-8.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-1/level-9.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "sea-green" notitle, \
     "../sphere-settling/case-2/level-7.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-2/level-8.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-2/level-9.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "sea-green" notitle, \
     "../sphere-settling/case-3/level-7.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "blue"      notitle, \
     "../sphere-settling/case-3/level-8.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "red"       notitle, \
     "../sphere-settling/case-3/level-9.dat" u 2:(($14)/0.015 - 0.5) w l lw 1.25 lc rgb "sea-green" notitle
~~~

## Results for $Re=31.9$

~~~gnuplot Time evolution of the settling velocity
reset
set terminal svg font ",16"
set key font ",16" top right spacing 1.1
set xlabel "t"
set ylabel "u_{p,y}"
set xrange [0:1.4]
set yrange [-0.14:0]
plot "../data/tenCate2002/tenCate2002-fig5b-Re31p9.csv" u 1:2 w p ps 0.7 pt 12 lc rgb "black" t "Re=31.9", \
     "../sphere-settling/case-3/level-7.dat" u 2:15 w l lw 1.25 lc rgb "blue"      t "Basilisk, l=7", \
     "../sphere-settling/case-3/level-8.dat" u 2:15 w l lw 1.25 lc rgb "red"       t "l=8", \
     "../sphere-settling/case-3/level-9.dat" u 2:15 w l lw 1.25 lc rgb "sea-green" t "l=9"
~~~

## References

~~~bib
@article{abraham1970,
  title={Functional Dependence of Drag Coefficient of a Sphere on Reynolds Number},
  author={Abraham, F.F.},
  journal={Physics of Fluids},
  volume={13},
  pages={2194--2195},
  year={1970}
}

@article{tenCate2002,
  title={Particle imaging velocimetry experiments and lattice-Boltzmann simulations on a single sphere settling under gravity},
  author={ten Cate, A. and Nieuwstad, C.H. and Derksen, J.J. and Van den Akker, H.E.A.},
  journal={Physics of Fluids},
  volume={14},
  pages={4012--4025},
  year={2002}
}

@article{wachs2011,
  title={PeliGRIFF, a parallel DEM-DLM/FD direct numerical simulation tool for 3D particulate flows},
  author={Wachs, A.},
  journal={J. Eng. Math},
  volume={71},
  pages={131--155},
  year={2011}
}
~~~
*/

