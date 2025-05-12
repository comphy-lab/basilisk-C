/**
# Stokes flow past a fixed cylinder at different Bingham numbers

We solve here the Stokes equations and add the cylinder using an
[embedded boundary](/src/embed.h). */

#if _MPI
#define JACOBI 1
#endif // _MPI

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered2.h"
#include "../myviscosity-viscoplastic.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define nu   (1.) // Viscosity
#define uref (1.) // Reference velocity, uref
#define tref (sq (d)/(nu)) // Reference time, tref=d/u

#if BI // 1, 100, 1000
#define Bi   ((double) (BI))
#else // Re = 100
#define Bi   (100.) // Bingham number
#endif // BI

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
coord p_p;

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
be taken into account. We also need a field for the density and the
yield stress. */

scalar rhov[];
face vector muv[], Tv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
We define the mesh adaptation parameters. */

#define lmin (7)  // Min mesh refinement level (l=7 is 2pt/d)
#define lmax (12) // Max mesh refinement level (l=10 is 16pt/d)
#define cmax (1.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  L0 = 64.*(d);
  size (L0);
  origin (-L0/2., -L0/2.);

  /**
  We set the maximum timestep. We tune the maximum time step to ensure
  convergence of the Poisson solver for large Bingham numbers. */

  DT = max (1.e-6, 1.e-4/(Bi))*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  stokes       = true;
  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  NITERMAX     = 500;

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
  foreach()
    rhov[] = cm[];
  boundary ({rhov});
  
  foreach_face() {
    muv.x[] = (nu)*fm.x[];
    Tv.x[]  = (Bi)*(uref)/(d)*fm.x[];
  }
  boundary ((scalar *) {muv, Tv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  rho = rhov;
  mu  = muv;

  /**
  We also set the yield stress and the regularization parameters. */

  T = new face vector;
  Tv = T;

  eps = min (5.e-2, max (1.e-4, 5.e-4*(Bi)));

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p})
    s.third = false;
#else
  for (scalar s in {u, p})
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
  We initialize the embedded boundary. */

  foreach_dimension()
    p_p.x = 0.;
  
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
  We also define the volume fraction at the previous timestep
  *csm1=cs*. */

  csm1 = cs;
    
  /**
  We define the boundary conditions for the velocity. */
   
  u.n[embed] = dirichlet (0);
  u.t[embed] = dirichlet (0);
  p[embed]   = neumann (0);

  uf.n[embed] = dirichlet (0);
  uf.t[embed] = dirichlet (0);
  pf[embed]   = neumann (0);
  
  /**
  We initialize the velocity. */

  foreach() 
    u.x[] = cs[]*(uref);
  boundary ((scalar *) {u});

  /**
  We finally initialize the reference velocity field. */
  
  foreach()
    un[] = u.x[];
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u,yuyz}, (double[]) {1.e-2,(cmax),(cmax),1.e-30},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. This would affect in particular the computation of
  the pressure contribution to the hydrodynamic forces. */

  p_shape (cs, fs, p_p);
}
#endif // TREE

/**
## Restarts and dumps

Every few characteristic time, we also dump the fluid data for
post-processing and restarting purposes. */

#if DUMP
event dump_data (i += 1000)
{
  // Dump fluid
  char name [80];
  sprintf (name, "dump-%d", i);
  dump (name);
}

event dump_end (t = end)
{
  // Dump fluid
  dump ("dump-final");  
}
#endif // DUMP

/**
## Profiling */

#if TRACE > 1
event profiling (i += 50)
{
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs */

/**
## Outputs */

double CDm1 = 0.;

event logfile (i++; i <= 10000)
{
  double du = change (u.x, un);

  coord Fp, Fmu;
  embed_force  (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/((nu)*(d));
  double CL = (Fp.y + Fmu.y)/((nu)*(d));

  double dF = fabs (CDm1 - CD);
  CDm1 = CD;  

  fprintf (stderr, "%g %g %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
	   Bi, eps,
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD/(Bi), CL,
	   du, dF);
  fflush (stderr);
}

/**
## Snapshots */

event snapshot (i += 1000)
{
  char name [80];

  /**
  We recompute the yielded region as the value of viscosity may have
  changed. */
  
  yielded_region ();    
  
  /**
  We first plot the entire domain. */
 
  view (fov = 20,
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name, "mesh-%d.png", i);
  save (name);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name, "ux-%d.png", i);
  save (name);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name, "uy-%d.png", i);
  save (name);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name, "p-%d.png", i);
  save (name);

  squares ("yuyz", min = 0, max = 1);
  sprintf (name, "yuyz-%d.png", i);
  save (name);
  
  /**
  We then zoom on the particle. */

  view (fov = 3,
	tx = -(p_p.x)/(L0), ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name, "mesh-zoom-%d.png", i);
  save (name);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name, "ux-zoom-%d.png", i);
  save (name);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name, "uy-zoom-%d.png", i);
  save (name);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name, "p-zoom-%d.png", i);
  save (name);

  squares ("yuyz", min = 0, max = 1);
  sprintf (name, "yuyz-zoom-%d.png", i);
  save (name);
}

event movies (i += 100)
{
  /**
  We recompute the yielded region as the value of viscosity may have
  changed. */
  
  yielded_region ();
  
  /**
  We first plot the entire domain. */
 
  view (fov = 20,
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  save ("mesh.mp4");
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  save ("ux.mp4");

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  save ("uy.mp4");

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  save ("p.mp4");

  squares ("yuyz", min = 0, max = 1);
  save ("yuyz.mp4");

  /**
  We then zoom on the particle. */

  view (fov = 3,
	tx = -(p_p.x)/(L0), ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  save ("mesh-zoom.mp4");
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  save ("ux-zoom.mp4");

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  save ("uy-zoom.mp4");

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  save ("p-zoom.mp4");

  squares ("yuyz", min = 0, max = 1);
  save ("yuyz-zoom.mp4");
}

/**
## Results

![Time evolution of the yielded region](cylinder-unbounded/yuyz.mp4)(loop)

![Time evolution of the yielded region (zoom)](cylinder-unbounded/yuyz-zoom.mp4)(loop)

#### Drag coefficient for $Bi=100$

~~~gnuplot Time evolution of the drag coefficient $C_D$
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics format "%.0e" 0,5000,10000 
set xlabel 'i'
set ylabel 'C_{D}/Bi'
set xrange [0:*]
set yrange [0:20]

# Topkavi correlation
f(x) = 11.98 + 20.43*x**(-0.68)
Bi = 100

plot f(Bi)         w l lw 2 lc rgb "black" t "Topkavi et al., 2008", \
     'log' u 3:16  w l lw 2 lc rgb "blue" t "Basilisk"
~~~

~~~gnuplot Time evolution of the relative error for the drag coefficient $C_D$
set ylabel '|C_{D}/Bi - C_{D}^{*}|/C_{D}^{*}'
set yrange [*:100]
set logscale y
plot 'log' u 3:(100.*abs ($16 - f(Bi))/f(Bi)) w l lw 2 lc rgb "blue" t "Basilisk"
~~~

## Results on a supercomputer (cedar)

#### Drag coefficient for different Bingham number

~~~gnuplot Drag coefficient $C_D/Bi$ for different Bingham numbers $Bi$
set xtics format "%.0f" 0,10,20000
set xlabel 'Bi'
set ylabel 'C_{D}/Bi'
set xrange [*:20000]
set yrange [10:35]
set logscale x
unset logscale y
plot f(x)           w l lw 2         lc rgb "black" t "Topkavi et al., 2008", \
     'Bi.dat' u 1:5 w p pt 5 ps 0.75 lc rgb "blue"  t "Basilisk"
~~~

~~~gnuplot Time evolution of the drag coefficient $C_D/Bi$ for different Bingham number $Bi$
set xtics format "%.0e" 0,5000,10000
set xlabel 'i'
set xrange [0:*]
unset logscale
plot 'Bi-1/log'    u 3:16 w l lw 2 lc rgb "black"     t "Bi=1",		\
     'Bi-10/log'   u 3:16 w l lw 2 lc rgb "blue"      t "Bi=10",	\
     'Bi-50/log'   u 3:16 w l lw 2 lc rgb "red"       t "Bi=50",	\
     'Bi-100/log'  u 3:16 w l lw 2 lc rgb "sea-green" t "Bi=100",	\
     'Bi-250/log'  u 3:16 w l lw 2 lc rgb "coral"     t "Bi=250",	\
     'Bi-500/log'  u 3:16 w l lw 2 lc rgb "magenta"   t "Bi=500",	\
     'Bi-750/log'  u 3:16 w l lw 2 lc rgb "brown"     t "Bi=750",	\
     'Bi-1000/log' u 3:16 w l lw 2 lc rgb "gray"      t "Bi=1000",	\
     'Bi-2500/log' u 3:16 w l lw 2                    t "Bi=2500",	\
     'Bi-5000/log' u 3:16 w l lw 2                    t "Bi=5000",	\
     'Bi-7500/log' u 3:16 w l lw 2                    t "Bi=7500",	\
     'Bi-10000/log' u 3:16 w l lw 2                   t "Bi=10000"
~~~

~~~gnuplot Time evolution of the relative error for drag coefficient $C_D$ for different Bingham number $Bi$
set ylabel 'err (%)'
set yrange [0:10]
plot 'Bi-1/log'    u 3:(100.*abs (f(1)  - $16)/(f(1)))      w l lw 2 lc rgb "black"     t "Bi=1", \
     'Bi-10/log'   u 3:(100.*abs (f(10)  - $16)/(f(10)))    w l lw 2 lc rgb "blue"      t "Bi=10",	\
     'Bi-50/log'   u 3:(100.*abs (f(50)  - $16)/(f(50)))    w l lw 2 lc rgb "red"       t "Bi=50", \
     'Bi-100/log'  u 3:(100.*abs (f(100) - $16)/(f(100)))   w l lw 2 lc rgb "sea-green" t "Bi=100", \
     'Bi-250/log'  u 3:(100.*abs (f(250) - $16)/(f(250)))   w l lw 2 lc rgb "coral"     t "Bi=250", \
     'Bi-500/log'  u 3:(100.*abs (f(500) - $16)/(f(500)))   w l lw 2 lc rgb "magenta"   t "Bi=500", \
     'Bi-750/log'  u 3:(100.*abs (f(750) - $16)/(f(750)))   w l lw 2 lc rgb "brown"     t "Bi=750", \
     'Bi-1000/log' u 3:(100.*abs (f(1000) - $16)/(f(1000))) w l lw 2 lc rgb "gray"      t "Bi=1000", \
     'Bi-2500/log' u 3:(100.*abs (f(2500) - $16)/(f(2500))) w l lw 2                    t "Bi=2500", \
     'Bi-5000/log' u 3:(100.*abs (f(5000) - $16)/(f(5000))) w l lw 2                    t "Bi=5000", \
     'Bi-7500/log' u 3:(100.*abs (f(7500) - $16)/(f(7500))) w l lw 2                    t "Bi=7500", \
     'Bi-10000/log' u 3:(100.*abs (f(10000) - $16)/(f(10000))) w l lw 2                 t "Bi=10000"
~~~
*/
