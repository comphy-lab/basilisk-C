/**
# Sphere moving infinitely slowly towards to a plane wall in a viscoplastic Stokes flow

We consider here a semi-infinite domain and compute the hydrodynamic
normal force $F_n$ acting on a sphere held fixed but with a velicity
$U$ at its boundary. The value of the force $F_n$ depends on the "gap"
distance $\delta$ between the sphere and the wall [Brenner,
1961](#Brenner1961):
$$
F_n/F_s = \frac{a}{\delta} -\frac{1}{5}\log (\frac{\delta}{a}) + 0.97128,
$$
where $a$ is the radius of the sphere and $F_s$ is the Stokes drag
defined as [Stokes, 1851](#Stokes1851):
$$
F_s = 6 \pi \mu U a.
$$

We solve here the 3D Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h) on an adaptive grid. */

/**
## Notes

A priliminary study showed that 6 points per gap length was enough to
obtain a relative error of the order of 1%. Additionnally, for large
gap sizes, the size of the domain plays an important role in the
accuray of the solution as the long range effect become more
important. */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../../yield/src/yield-regularization.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)  // Diameter, 2*a
#define nu   (1.)  // Viscosity
#if GAP
#define gap  (GAP/1000.) // Dimensionless gap distance delta/a (0.4, 0.3, 0.2, 0.1, 0.05, 0.025 and 0.01)
#else // gap = 0.4
#define gap  (0.4)
#endif // GAP
#define uref (1.)
#define tref (sq (d)/(nu))
#define Ftol (1.e-4) // Tolerance on the normal force for the steady state

#define Fstokes  (6.*pi*(nu)*(uref)*(d)/2.) // Stokes drag
#define FnsFs(g) (1./(g) - 1./5.*log((g)) + 0.97128) // Normal force divided by Stokes drag

/**
We also define the shape of the domain. 

To avoid a high concentration of cells near the boundary and to
increase the accuracy of the boundary treatment, we define the wall
towards which the particle is "moving" using embedded boundaries
instead of using a boundary of the domain. */

#if DLENGTH
#define h (-((double) (DLENGTH))/2. + 1.1*(d) + 1.e-8) // Position of the wall
#else // L0=256
#define h (-(L0)/2. + 1.1*(d) + 1.e-8) // Position of the wall
#endif // DLENGTH
#define sphere(x,y,z,diam) (sq ((x)) + sq ((y)) + sq ((z)) - sq ((diam)/2.))
#define wall(x,w) ((x) - (w)) // + over, - under
coord p_p;

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = intersection (
    			  (sphere ((x - p.x), (y - p.y), (z - p.z), (d))),
    			  (wall ((x), (h)))
    			  );
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
We then need a field for the yield stress T. */

#define T1 (100./sqrt (2.)) // Yield stress
face vector Tv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
Finally, we define the mesh adaptation parameters.

We choose the adaptation criteria *cmax* as 0.5% of the settling
velocity *uref*. 

We first set the minimum level of refinement *lmin* */

#if DLENGTH == 8
#define lmin (4) // Min mesh refinement level (l=4 is 2pt/d for L0=8)
#elif DLENGTH == 16
#define lmin (5) // Min mesh refinement level (l=5 is 2pt/d for L0=16)
#elif DLENGTH == 32
#define lmin (6) // Min mesh refinement level (l=6 is 2pt/d for L0=32)
#elif DLENGTH == 64
#define lmin (7) // Min mesh refinement level (l=7 is 2pt/d for L0=64)
#elif DLENGTH == 128
#define lmin (8) // Min mesh refinement level (l=8 is 2pt/d for L0=128)
#elif DLENGTH == 256
#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d for L0=256)
#else // L0 = 256
#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d for L0=256)
#endif // DLENGTH

/**
We then set the maximum level of refinement *lmax. */

#if LMAX
#define lmax ((int) LMAX) // Max mesh refinement level (l=13 is 32pt/d for L0=256)
                    // For gap=0.4,   l=13 is 6pt/delta)
                    // For gap=0.3,   l=14 is 9pt/delta)
                    // For gap=0.2,   l=14 is 6pt/delta)
                    // For gap=0.1,   l=15 is 6pt/delta)
                    // For gap=0.05,  l=16 is 6pt/delta)
                    // For gap=0.025, l=17 is 6pt/delta)
                    // For gap=0.01,  l=18 is 5pt/delta)
#else
#if DLENGTH == 8
#define lmax (8)  // Max mesh refinement level (l=8  is 32pt/d for L0=8)
#elif DLENGTH == 16
#define lmax (9)  // Max mesh refinement level (l=9  is 32pt/d for L0=16)
#elif DLENGTH == 32
#define lmax (10) // Max mesh refinement level (l=10 is 32pt/d for L0=32)
#elif DLENGTH == 64
#define lmax (11) // Max mesh refinement level (l=11 is 32pt/d for L0=64)
#elif DLENGTH == 128
#define lmax (12) // Max mesh refinement level (l=12 is 32pt/d for L0=128)
#elif DLENGTH == 256
#define lmax (13) // Max mesh refinement level (l=13 is 32pt/d for L0=256)
#else // L0 = 256
#define lmax (13) // Max mesh refinement level (l=13 is 32pt/d for L0=256)
#endif // DLENGTH
#endif // LMAX
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field
#define lmem (6) // Initial mesh refinement level to avoid memory overload

/**
We finally define a useful function that allows us to define block
mesh refinement (BMR) around the sphere. */

#define cube(x,hx,y,hy,z,hz) (union (					\
				     union (				\
					    union ((x) - (hx)/2., -(x) - (hx)/2.), \
					    union ((y) - (hy)/2., -(y) - (hy)/2.)), \
				     union ((z) - (hz)/2., -(z) -(hz)/2.)))

int main ()
{  
  /**
  The domain is $256\times 256 \times 256$ by default. */

#if DLENGTH // 32, 64, 128, 256
  L0 = ((double) (DLENGTH))*(d);
#else // L0 = 256
  L0 = 256.*(d);
#endif // DLENGTH
  size (L0);
  origin (-(L0)/2., -(L0)/2., -(L0)/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-3*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  stokes       = true;
  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  NITERMAX     = 200; // Convergence is more difficult at smaller gap sizes
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmem);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions */

u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann (0);

u.n[right] = dirichlet (0);
u.t[right] = dirichlet (0);
u.r[right] = dirichlet (0);
p[right]   = neumann (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
u.r[bottom] = dirichlet (0);
p[bottom]   = neumann (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);
u.r[top] = dirichlet (0);
p[top]   = neumann (0);

u.n[back] = dirichlet (0);
u.t[back] = dirichlet (0);
u.r[back] = dirichlet (0);
p[back]   = neumann (0);

u.n[front] = dirichlet (0);
u.t[front] = dirichlet (0);
u.r[front] = dirichlet (0);
p[front]   = neumann (0);

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
  foreach_face() {
    muv.x[] = (nu)*fm.x[];
    Tv.x[]  = (T1)*fm.x[];
  }
  boundary ((scalar *) {muv, Tv});
}

/**
## Initial conditions 

We plot the initial mesh. */

event init (i = 0)
{
  char name2[80];

  /**
  We first plot the whole domain. */

  view (fov = 20,
	tx = 1.e-12, ty = 0.,
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "mesh-t-%.0f.png", t);
  save (name2);

  view (fov = 10,
	tx = ((L0)/4.)/(L0), ty = 0.,
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "mesh-zoom-t-%.0f.png", t);
  save (name2);
  
  /**
  We plot the mesh and embedded boundaries. */

  view (fov = 0.2*(256./(L0)),
	tx = -(p_p.x)/(L0), ty = -(p_p.y)/(L0),
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "vof-t-%.0f.png", t);
  save (name2);

  /**
  We finally plot the yielded region, defined by *yuyz=0*. */
  
  squares ("yuyz", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "yuyz-t-%.0f.png", t);
  save (name2);
}

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;

  /**
  We also set the yield stress and the regularization parameters. */

  T = new face vector;
  Tv = T;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p})
    s.third = false;
#else
  for (scalar s in {u, p})
    s.third = true;
#endif // ORDER2
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the full cell center of
  cut-cells. */
#endif // TREE

  /**
  We initialize the embedded boundary. */

  /**
  We first define the sphere's position. */
  
  p_p.x = (h) + ((d)/2. + (gap)*(d)/2.);
  p_p.y = 0.;
  p_p.z = 0.;

  /**
  If the simulation is not restarted, we define the initial mesh and
  the initial velocity. */

  if (!restore (file = "restart")) {

#if TREE
#if BMR
    /**
    When using BMR, we refine the mesh in user-defined rectangles
    around the embedded boundary. */

    // We induce here a spherical symmetry in the mesh that may help the Poisson solver converge.
    // Initial coarse grid refinement, depending on lmem
    if ((lmem) < (lmin)) {
      refine ((sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       2.*((L0)/(1 << (lmem))))) < 0. &&
	      level < (lmin));
      p_shape (cs, fs, p_p);
    }
    // Coarse grid far from the sphere
    unrefine ((sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       (3.*(d)))) > 0. &&
	      level > (1));
    p_shape (cs, fs, p_p);
    // Fine grid (lmax) around the sphere and coarse grid far from the sphere
    int lvl = (lmin);
    while (lvl < (lmax)) {
      refine ((sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       ((d) + 2.*(L0)/(1 << (lvl))))) < 0. &&
	      (sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       ((d) - 2.*(L0)/(1 << (lvl))))) > 0. &&
	      level < (lvl + 1));
      p_shape (cs, fs, p_p);
      lvl ++;
    }
    // Fine grid in the gap and around the sphere
    lvl = (lmin);
    while (lvl < (lmax)) {
      refine ((wall ((x),   ((h) - (d)/2.*(gap)))) >= 0. &&
	      (sphere ((x - ((h) - (d)/2.*(gap))),
		       (y - p_p.y),
		       (z - p_p.z),
		       ((d)*(4. + (1. - 4.)/((lmax) - (lmin))*((lvl + 1) - (lmin)))))) <= 0. &&
	      (sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       (d))) > 0. &&
	      level < (lvl + 1));
      p_shape (cs, fs, p_p);
      lvl ++;
    }
    // Unrefine in the sphere and in the wall
    unrefine ((
    	       (sphere ((x - p_p.x),
    			(y - p_p.y),
    			(z - p_p.z),
    			((d) - 2.*(L0)/(1 << (lmax))))) < 0.  ||
    	       (wall ((x), ((h) - 2.*(L0)/(1 << (lmax))))) < 0.) &&
    	      level > (1));
#else // AMR
    /**
    Before defining the embedded boundaries, we refine the mesh up to
    *lmin* in the vicinity of the particle. */

    if ((lmem) < (lmin)) {
      assert ((d)/2. + (gap)*(d)/2. < (L0)/(1 << (lmem)));
      refine ((sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       2.*((L0)/(1 << (lmem))))) < 0. &&
	      level < (lmin));
    }
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
#endif // BMR
#endif // TREE
  }

  /**
  After initialization or restart, we still need to define the face
  fraction *fs* as it is not dumped. */
  
  p_shape (cs, fs, p_p);
  
  /**
  Whether restarting or not, we define the volume fraction at the
  previous timestep *csm1=cs*. */
    
  csm1 = cs;
    
  /**
  We define the boundary condition $-U\mathbf{e_x}$ for the
  velocity. */
      
  u.n[embed] = dirichlet (x > (h) + (gap)*(d)/4. ? -(uref) : 0.);
  u.t[embed] = dirichlet (0);
  u.r[embed] = dirichlet (0);
  p[embed]   = neumann (0);

  uf.n[embed] = dirichlet (x > (h) + (gap)*(d)/4. ? -(uref) : 0.);
  uf.t[embed] = dirichlet (0);
  uf.r[embed] = dirichlet (0);

  /**
  We finally initialize, even when restarting, the reference velocity
  field. */
  
  foreach()
    un[] = u.y[];
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE && !BMR
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-1,(cmax),(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. */

  p_shape (cs, fs, p_p);
}
#endif // TREE && !BMR

/**
## Restarts and dumps

Every few characteristic time, we also dump the fluid data for
post-processing and restarting purposes. */

#if DUMP
event dump_data (t += 2.*(tref))
{
  // Dump fluid
  char name [80];
  sprintf (name, "dump-%g", t/(tref));
  dump (name);
}
#endif // DUMP

event dump_end (t = end)
{
  // Dump fluid
  dump ("dump-final");  
}

/**
## Profiling */

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs 

We define a color volume fraction in order to compute the force acting
only on the sphere. */

scalar col[];

int nc = 0;
double Fm1 = 0.;

event init (i = 0)
{
  trash ({col});
  nc = 0;
  Fm1 = 0.;
}

event logfile (i++; t <= 10.*(tref))
{
  /**
  We update the color volume fraction. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sphere ((x - p_p.x), (y - p_p.y), (z - p_p.z), (d)));
  boundary ({phi});
  fractions (phi, col);
  boundary ((scalar *) {col});

  /**
  We compute the hydrodynamic forces acting on the sphere. */

  coord Fp, Fmu;
  embed_color_force (p, u, mu, col, &Fp, &Fmu);
  double F  = fabs (Fp.x + Fmu.x)/(Fstokes);
  double dF = fabs (F - Fm1);
  Fm1 = F;
  
  fprintf (stderr, "%g %g %d %d %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
	   (gap), (L0)/(d), (lmin), (lmax),
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   (FnsFs ((gap))), F, fabs (F - (FnsFs ((gap)))), dF
	   );
  fflush (stderr);

  /**
  We stop the simulation if the relative error *dF* is smaller than
  *Ftol*, and if this occurs for 5 iterations in a row. */

  if (dF <= (Ftol))
    nc++;
  else
    nc = 0;

  if (t > 2.*(tref) && dF <= (Ftol) && nc > 5) {

    // Dump fluid
    char name [80];
    sprintf (name, "dump-final");
    dump (name);
    
    return 1; /* stop */
  }
}

/**
## Snapshots */

event snapshot (t += 1.*(tref))
{
  char name2[80];

  /**
  We first plot the whole domain. */

  clear ();
  view (fov = 20,
	tx = 0., ty = 0.,
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "mesh-t-%.0f.png", t/(tref));
  save (name2);
  
  /**
  We plot the mesh and embedded boundaries. */

  clear ();
  view (fov = 0.2*(256./(L0)),
	tx = -(p_p.x)/(L0), ty = -(p_p.y)/(L0),
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "vof-t-%.0f.png", t/(tref));
  save (name2);

  /**
  We finally plot the yielded region, defined by *yuyz=0*. */

  squares ("yuyz", n = {0,0,1}, alpha = 0, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 0);
  sprintf (name2, "yuyz-t-%.0f.png", t);
  save (name2);  
}

/**
## Results

#### Results for $\frac{\delta}{a} = 0.4$

![Initial mesh around the embedded boundaries](sphere-near-wall-viscoplastic/vof-t-0.png)

~~~gnuplot Time evolution of the normal force *Fn*
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,1,10
set ytics 0,2,10
set xlabel 't/(d^2/{\nu})'
set ylabel 'F_n/F_s'
set xrange [0:5]
set yrange [0:10]
plot 1           w l lw 2 lc rgb 'black' t 'Stokes', \
     'log' u 6:18 w l lw 2 lc rgb 'blue'  t 'Analytic, Brenner, 1961', \
     ''    u 6:19 w l lw 2 lc rgb 'red'   t 'Basilisk, {\delta}/a=0.4'
~~~

## References

~~~bib
@article{Stokes1851,
  title={On the effect of the internal friction of fluids on the motion of pendulums},
  author={Stokes, G.G.},
  year={1851},
  publisher={Cambridge Pitt Press}
}
@article{Brenner1961,
  title={The slow motion of a sphere through a viscous fluid towards a plane surface},
  author={Brenner, H.},
  journal={Chemical Engineering Science},
  volume={16},
  pages={242--251},
  year={1961}
}
~~~
*/
