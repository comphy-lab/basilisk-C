/**
# Flow past a fixed sphere at $Re=20$

We solve here the Navier-Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myperfs.h"

#if BVIEW
#include "view.h"
#endif // BVIEW

/**
## Reference solution */

#define d    (1.)         // Diameter of the sphere
#define Re   (20.)        // Particle Reynolds number Re = ud/nu
#define uref (1.)         // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

/**
We also define the shape of the domain. Note here that *p_p* is the
position of the center of the sphere. */

#define sphere(x,y,z) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((d)/2.))
coord p_p;

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar psi[];
  foreach_vertex()
    psi[] = (sphere ((x - p.x), (y - p.y), (z - p.z)));
  boundary ({psi});
  fractions (psi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
Finally, we define the mesh adaptation parameters. */

#define lmin (5) // Min mesh refinement level (l=5 is 2pt/d for L0=16)
#define lmax (9) // Max mesh refinement level (l=9 is 32pt/d for L0=16)
#define cmax (1.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{
  /**
  The domain is $16\times 16 \times 16$. */

  L0 = 16.*(d);
  size (L0);
  origin (-(L0)/2., -(L0)/2., -(L0)/2.);
  
  /**
  We set the maximum timestep. */

  DT = 1./2560.*(tref);
 
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
  should also define the gradient of *u* at the full cell center of
  cut-cells. */
#endif // TREE

  /**
  We first define the particle's position. */

  foreach_dimension()
    p_p.x = 0.;

  /**
  We then define the initial mesh and the initial velocity. */

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
    u.x[] = cs[]*(uref);
  boundary ((scalar *) {u});

  /**
  We define the volume fraction at the previous timestep *csm1=cs*. */
    
  csm1 = cs;
    
  /**
  We define the boundary conditions for the velocity. */
   
  u.n[embed] = dirichlet (0);
  u.t[embed] = dirichlet (0);
  u.r[embed] = dirichlet (0);
  p[embed]   = neumann (0);

  uf.n[embed] = dirichlet (0);
  uf.t[embed] = dirichlet (0);
  uf.r[embed] = dirichlet (0);
  pf[embed]   = neumann (0);
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
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. This would affect in particular the computation of
  the pressure contribution to the hydrodynamic forces. */

  p_shape (cs, fs, p_p);
}
#endif // TREE

/**
## Profiling */

#if TRACE > 1
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs */

double CDm1 = 0., Tzm1 = 0.;

event logfile (i++; t < 10.*(tref))
{
  coord Fp, Fmu, Tp, Tmu;
  embed_force  (p, u, mu, &Fp, &Fmu);
  embed_torque (p, u, mu, (p_p), &Tp, &Tmu);

  // Drag and lifts
  double CD  = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);
  double CLy = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);
  double CLz = (Fp.z + Fmu.z)/(0.5*sq ((uref))*(M_PI)*sq ((d))/4.);

  // Torques and pitching torque
  double Tx = (Tp.x + Tmu.x)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);
  double Ty = (Tp.y + Tmu.y)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);
  double Tz = (Tp.z + Tmu.z)/(0.5*sq ((uref))*(M_PI)*sq ((d))/8.);

  double dF = fabs (CDm1 - CD);
  double dT = fabs (Tzm1 - Tz);
  CDm1 = CD, Tzm1 = Tz;  

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %ld %g\n",
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CLy, CLz, dF,
	   Tx, Ty, Tz, dT,
	   grid->tn, perf.t);
  fflush (stderr);
}

/**
## Visualization */

#if BVIEW
event snapshots (t = {0., 1.*(tref), 2.*(tref), 5.*(tref), 10.*(tref)})
{
  // norm of velocity
  scalar n2u[];
  foreach() {
    if (cs[] <= 0.)
      n2u[] = nodata;
    else
      n2u[] = sqrt (sq (u.x[]) + sq (u.y[]) + sq (u.z[]));
  }
  boundary ({n2u});

  stats sn2u = statsf (n2u);
  stats sux  = statsf (u.x);
  stats sp   = statsf (p);

  char name1[80];

  /**
  We first plot the whole domain. */

  clear ();
  view (fov = 20, camera = "front",
	tx = 1.e-12, ty = 1.e-12,
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells   (      n = {0,0,1}, alpha = 1.e-12);
  sprintf (name1, "vof-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 1);
  squares  ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name1, "nu-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 1);
  squares  ("u.x", n = {0,0,1}, alpha = 1.e-12, min = sux.min, max = sux.max, map = cool_warm);
  sprintf (name1, "ux-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 1);
  squares  ("p", n = {0,0,1}, alpha = 1.e-12, min = sp.min, max = sp.max, map = cool_warm);
  sprintf (name1, "p-t-%.0f.png", t/(tref));
  save (name1);

  /**
  We plot the mesh and embedded boundaries. */

  clear ();
  view (fov = 4, camera = "front",
	tx = -(p_p.x)/(L0), ty = -(p_p.y)/(L0),
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells   (      n = {0,0,1}, alpha = 1.e-12);
  sprintf (name1, "vof-xy-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 1);
  squares  ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name1, "nu-xy-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 1);
  squares  ("u.x", n = {0,0,1}, alpha = 1.e-12, min = sux.min, max = sux.max, map = cool_warm);
  sprintf (name1, "ux-xy-t-%.0f.png", t/(tref));
  save (name1);

  draw_vof ("cs", "fs", lw = 5);
  squares  ("p", n = {0,0,1}, alpha = 1.e-12, min = sp.min, max = sp.max, map = cool_warm);
  sprintf (name1, "p-xy-t-%.0f.png", t/(tref));
  save (name1);
}
#endif // BVIEW

/**
## Results for a sphere at $Re=20$

~~~gnuplot Time evolution of the drag coefficient $C_D$
reset
set terminal svg font ",16"
set key top right font ",16" spacing 1.1
set xtics 0,1,10
set ytics 0,1,100
set xlabel 't/(d/u_{inf})'
set ylabel 'C_D'
set xrange [0:10]
set yrange [1:5]
set title "Re=20, L0/d=16, 32pt/d, dt/tref=1/2560"

# Clift et al., 2005
Re = 20.;
CD = 24./Re*(1. + 0.1315*Re**(0.82 - 0.05*log10(Re)));

plot CD           w l lw 2 lc rgb "black" t "Clift et al., 2005",		\
     "log" u 2:14 w l lw 2 lc rgb "blue"  t "basilisk"
~~~

~~~gnuplot Time evolution of the error for the drag coefficient $C_D$
set ytics 0,0.5,100
set ylabel 'err (%)'
set yrange [0:7]
set title "Re=20, L0/d=16, 32pt/d, dt/tref=1/2560"
plot "log" u 2:(abs($14 - CD)/CD*100.) w l lw 2 lc rgb "blue" t "basilisk"
~~~

~~~gnuplot Time evolution of the number of cells
set ytics 0,0.5,100
set ylabel '# of cells (%)'
set yrange [0:10]
set title "Re=20, L0/d=16, 32pt/d, dt/tref=1/2560"

# Reference data
lmax = 9;

plot "log" u 2:($22/((2**lmax)**3)*100.) w l lw 2 lc rgb "blue" t "basilisk"
~~~

## References

~~~bib
@article{Clift2005,
  title={Bubbles, drops, and particles},
  author={Clift, Roland and Grace, John R and Weber, Martin E},
  year={2005},
  publisher={Courier Corporation}
}
*/
