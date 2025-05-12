/**
# Sphere moving infinitely slowly towards to a plane surface in a Stokes flow

We compute here the hydrodynamic normal force $F_{\Gamma,x}$ acting on
a sphere of radius $r$ (diameter $d$) approaching the plane surface
$x=-\frac{L_0}{2} + 1.1 d$ with a velocity $\mathbf{u}_{\Gamma} =
-u_{\mathrm{ref}} \mathbf{e_x}$ through a viscous fluid characterized
by a density $\rho = 1$ and a viscosity $\nu = 1$. We denote $\delta$
the distance between the plane and the bottom of the sphere.

This problem has been solved analytically in [Brenner,
1961](#Brenner1961) by considering the characteristic time scale
$t_{\mathrm{ref}} = \frac{\delta}{u_{\mathrm{ref}}}$ and assuming that
the Reynolds number verifies $Re = \frac{u_{\mathrm{ref}} r}{\nu} \ll
1$ and $Re \ll \frac{\delta}{r}$ to respectively neglect the inertia
and unsteady terms in the Navier-Stokes equations. See [Coolley et
al., 1969](#Cooley1969) for details. In this case, the normal force
acting on the sphere for any value of $\delta$ writes:

$$
F_{\Gamma,x} = \kappa F_{\mathrm{Stokes}},
$$
where $F_{\mathrm{Stokes}}$ is the well-known Stokes force defined in
\citep{Stokes1851} as:
$$
F_{\mathrm{Stokes}} = 6 \pi \mu u_{\mathrm{ref}} r,
$$
and $\kappa$ is a correction factor defined in \citep{Brenner1961} as:
$$
\kappa_{\mathrm{Brenner}} = \frac{4}{3} \sinh \alpha
\sum_{n=1}^{+\infty} \frac{n\left(n+1\right)}{\left(2n -
1\right)\left(2n + 3\right)} \left( \frac{2\sinh\left(\left(2n +
1\right)\alpha\right) + \left(2n + 1\right)\sinh
\left(2\alpha\right)}{4\sinh^2\left(\left(n +
\frac{1}{2}\right)\alpha\right) - \left(2n + 1\right)^2\sinh^2 \alpha}
- 1\right),
$$
where $\alpha = \cosh^{-1} \left( 1 + \frac{\delta}{r}\right)$. In
[Coolley et al., 1969](#Cooley1969), the authors have simplified the
previous expression by assuming that $\delta \ll r$ and have obtained
following expression for the correction factor $\kappa$:
$$
\kappa_{\mathrm{Cooley}} = \frac{r}{\delta} -\frac{1}{5}\log
(\frac{\delta}{r}) + 0.97128.
$$
Both expressions are equivalent in the limit of $\delta \ll r$ and the
authors report in [Coolley et al., 1969](#Cooley1969) a $3\%$
difference between the two expressions when $\delta /r \sim 0.5$.

We solve here the 3D Stokes equations and add the sphere using an
[embedded boundary](/src/embed.h) on an adaptive grid. */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myperfs.h"
#include "view.h"

/**
## Exact solution 

This table summerizes the results of [Brenner, 1961](#Brenner1961),
where the first column is the non-dimensional gap size $\delta/r$ and
the second column is the corresponding value of the drag correction
factor $\kappa$ (computed using the first 400 terms of the infinite
series). */

static double brenner[7][2] = {
  {0.4,   3.73562},
  {0.3,   4.61072},
  {0.2,   6.34089},
  {0.1,   11.4592},
  {0.05,  21.5858},
  {0.025, 41.7176},
  {0.01,  101.896}
};

double drag_correction_Brenner (double g)
{
  double kappa = 1.;
  for (int i = 0; i < 7; i++)
    if (g == brenner[i][0])
      kappa = brenner[i][1];
  return kappa;
}

/**
The following function defines the drag correction factor from [Cooley
et al., 1969](#Cooley1969). */

double drag_correction_Cooley (double g)
{
  return (1./(g) - 1./5.*log((g)) + 0.97128);
}

/**
## Reference solution */

#define d    (1.)  // Diameter, 2*r
#define nu   (1.)  // Viscosity

#if GAP
#define gap  (GAP/1000.) // Dimensionless gap distance delta/r (0.4, 0.3, 0.2, 0.1, 0.05, 0.025 and 0.01)
#else // gap = 0.4
#define gap  (0.4)
#endif // GAP

#define uref (1.)
#define tref (min (sq ((d)/(nu)), (gap)*(d)/2./(uref))) // tref = delta/u

/**
We also define the shape of the domain. 

To avoid a high concentration of cells near a domain boundary
(pathological case for Basilisk) and to increase the accuracy of the
boundary treatment, we define the wall towards which the particle is
"moving" using embedded boundaries instead of using a boundary of the
domain. */

#if DLENGTH
#define h (-((double) (DLENGTH))*(d)/2. + 1.1*(d) + 1.e-8) // Position of the wall
#else // L0 = 256
#define h (-(256.)*(d)/2. + 1.1*(d) + 1.e-8) // Position of the wall
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
We also define a reference velocity field. */

scalar un[];

/**
Finally, we define the mesh adaptation parameters. */

#define lmem (6) // Initial mesh refinement level to avoid memory overload

#if DLENGTH == 16
#define lmin (5) // Min mesh refinement level (l=7 is 2pt/d for L0=16)
#elif DLENGTH == 32
#define lmin (6) // Min mesh refinement level (l=8 is 2pt/d for L0=32)
#elif DLENGTH == 64
#define lmin (7) // Min mesh refinement level (l=7 is 2pt/d for L0=64)
#elif DLENGTH == 128
#define lmin (8) // Min mesh refinement level (l=8 is 2pt/d for L0=128)
#elif DLENGTH == 256
#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d for L0=256)
#else // L0 = 256
#define lmin (9) // Min mesh refinement level (l=9 is 2pt/d for L0=256)
#endif // DLENGTH

#if LMAX
#define lmax ((int) (LMAX)) // Max mesh refinement level (l=13 is 32pt/d for L0=256)
                    // For gap=0.4,   l=13 is 6pt/delta)
                    // For gap=0.3,   l=14 is 9pt/delta)
                    // For gap=0.2,   l=14 is 6pt/delta)
                    // For gap=0.1,   l=15 is 6pt/delta)
                    // For gap=0.05,  l=16 is 6pt/delta)
                    // For gap=0.025, l=17 is 6pt/delta)
                    // For gap=0.01,  l=18 is 5pt/delta)
#else
#define lmax (13) // Max mesh refinement level (l=13 is 32pt/d for L0=256)
#endif // LMAX

#if CMAX
#define cmax (((double) (CMAX))*1.e-3*(uref))
#else // cmax = 1.e-2
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field
#endif // CMAX

int main ()
{  
  /**
  The domain is $256\times 256 \times 256$ by default. */

#if DLENGTH // 16, 32, 64, 128, 256
  L0 = ((double) (DLENGTH))*(d);
#else // L0 = 256
  L0 = 256.*(d);
#endif // DLENGTH
  size (L0);
  origin (-(L0)/2., -(L0)/2., -(L0)/2.);

  /**
  We set the maximum timestep. */

#if DTMAX
  DT = ((double) (DTMAX))*1.e-3*(tref);
#else // DT = 1.e-2  
  DT = 1.e-2*(tref);
#endif // DTMAX
 
  /**
  We set the tolerance of the Poisson solver. */

  stokes       = true;
  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  NITERMAX     = 200; // Convergence is more difficult for smaller gap sizes
  
  /**
  We initialize the grid. */
  
  N = 1 << (lmem);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We impose no-slip boundary conditions on all boundaries as the fluid
is supposed to be at rest at an infinite distance from the plane. */

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
  for (scalar s in {u, p})
    s.third = false;
#else
  for (scalar s in {u, p})
    s.third = true;
#endif // ORDER2
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
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
  the initial velocity. 

  Note here that we do not dump the pressure so the first drag
  computation after a restart will be incorrect. */

  if (!restore (file = "restart")) {

#if TREE
    /**
    Before defining the embedded boundaries, we refine the mesh up to
    *lmin* in the vicinity of the particle. */

    if ((lmem) < (lmin)) {
      refine ((sphere ((x - p_p.x),
		       (y - p_p.y),
		       (z - p_p.z),
		       2.*((L0)/(1 << (lmem))))) < 0. &&
	      level < (lmin));
      p_shape (cs, fs, p_p);
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
  We define the boundary condition $-u_{\mathrm{ref}}\mathbf{e_x}$ for
  the velocity. */
      
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
    un[] = u.x[];
}

/**
## Embedded boundaries */

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-1,(cmax),(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. */

  p_shape (cs, fs, p_p);
}
#endif // TREE

/**
## Restarts and dumps

Every 100 time steps, we dump the fluid data for restarting
purposes. We use a relatively small time step interval as the
simulations with large values of *lmax* are relatively slow. */

#if DUMP
event dump_data (i += 100)
{
  // Dump fluid
  dump ("dump");
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
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs 

We define a color volume fraction in order to compute the force acting
only on the sphere. */

scalar col[];

double Fm1 = 0.;

event init (i = 0)
{
  trash ({col});
  Fm1 = 0.;
}

event logfile (i++; t <= 30.*(tref))
{
  double du = change (u.x, un);
  
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
  double F  = fabs (Fp.x + Fmu.x)/(6.*pi*(nu)*(uref)*(d)/2.);
  double dF = fabs (F - Fm1);
  Fm1 = F;

  /**
  We use the default *log* file to output the results. However, this
  file will be written over if the simulation is restarted. */
  
  fprintf (stderr, "%g %g %d %d %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   (gap), (L0)/(d), (lmin), (lmax),
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   Fp.x, Fp.y, Fp.z,
	   Fmu.x, Fmu.y, Fmu.z,
	   F,
	   (drag_correction_Brenner ((gap))), fabs (F - (drag_correction_Brenner ((gap)))),
	   (drag_correction_Cooley ((gap))),  fabs (F - (drag_correction_Cooley ((gap)))),
	   dF, du
	   );
  fflush (stderr);
}

/**
## Snapshots */

event snapshot (t += 5.*(tref))
{
  scalar n2u[];
  foreach() {
    if (cs[] <= 0.)
      n2u[] = nodata;
    else
      n2u[] = sqrt (sq (u.x[]) + sq (u.y[]) + sq (u.z[]));
  }
  boundary ({n2u});
  stats sn2u = statsf (n2u);
  
  char name2[80];

  /**
  We first plot the whole domain. */

  clear ();
  view (fov = 20, camera = "front",
	tx = 1.e-12, ty = 1.e-12,
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 1.e-12);
  sprintf (name2, "vof-t-%.0f.png", t/(tref));
  save (name2);

  squares ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name2, "nu-t-%.0f.png", t/(tref));
  save (name2);
  
  /**
  We plot the mesh and embedded boundaries. */

  //xy

  clear ();
  view (fov = 0.2*(256.*(d)/(L0)), camera = "front",
	tx = -(p_p.x)/(L0), ty = -(p_p.y)/(L0),
	bg = {1,1,1}, width = 800, height = 800);
  
  squares ("cs", n = {0,0,1}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells (n = {0,0,1}, alpha = 1.e-12);
  sprintf (name2, "vof-xy-t-%.0f.png", t/(tref));
  save (name2);
  
  squares ("n2u", n = {0,0,1}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name2, "nu-xy-t-%.0f.png", t/(tref));
  save (name2);

  // xz

  clear ();
  view (fov = 0.2*(256.*(d)/(L0)), camera = "top",
	tx = -(p_p.x)/(L0), ty = -(p_p.z)/(L0),
	bg = {1,1,1}, width = 800, height = 800);

  squares ("cs", n = {0,1,0}, alpha = 1.e-12, min = 0, max = 1, map = cool_warm);
  cells (n = {0,1,0}, alpha = 1.e-12);
  sprintf (name2, "vof-xz-t-%.0f.png", t/(tref));
  save (name2);

  squares ("n2u", n = {0,1,0}, alpha = 1.e-12, min = 0, max = sn2u.max, map = jet);
  sprintf (name2, "nu-xz-t-%.0f.png", t/(tref));
  save (name2);
}

/**
## Results

#### Results for $\frac{\delta}{r} = 0.3$ using *L0/d=128* and *lmax=12*

![Initial mesh and embedded boundaries](sphere-towards-wall/vof-t-0.png)

![Initial velocity](sphere-towards-wall/nu-xy-t-0.png)

~~~gnuplot Time evolution of the normal force $F_x$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,5,100
set ytics 0,2,10
set xlabel 't/t_{ref}'
set ylabel 'kappa'
set xrange [0:30]
set yrange [3:6]
set title "delta/r=0.3, dt=1e-2, delta/Delta_{max}=6"
plot 'log' u 6:25 w l lw 2 lc rgb 'blue'      t 'Brenner, 1961', \
     'log' u 6:27 w l lw 2 lc rgb 'red'       t 'Cooley et al., 1969', \
     ''    u 6:24 w l lw 2 lc rgb 'sea-green' t 'Basilisk'
~~~

~~~gnuplot Time evolution of the relative error for the normal force $F_x$
set ytics 1.e-8,1.e-1,100
set ylabel 'err (%)'
set yrange [1.e-3:1e2]
set logscale y
plot 'log' u 6:(100.*abs($27 - $25)/$25) w l lw 2 lc rgb 'black' t 'Error between Brenner and Cooley', \
     ''    u 6:(100.*$26/$25)            w l lw 2 lc rgb 'blue'  t 'Error with Brenner, 1961', \
     ''    u 6:(100.*$28/$27)            w l lw 2 lc rgb 'red'   t 'Error with Cooley et al., 1969'
~~~
*/

/**
## Plots specific to paper

#### Domain size study for $\delta/r=0.4$ and *DT=1e-2* and *cmax=1e-2*

~~~gnuplot Time evolution of the normal force $F_x$
reset
set terminal svg font ",16"
set key top right font ",14" spacing 1
set grid ytics
set xtics 0,5,100
set ytics 0,0.1,10
set xlabel 't/t_{ref}'
set ylabel 'kappa'
set xrange [0:30]
set yrange [3.6:4]
set title 'delta/r=0.4, dt=1e-2, delta/Delta_{max}=12, cmax=1e-2'

# Average correction factor kappa

ftitle(a) = sprintf(", kappa_{inf}=%2.4f", a);

f1(x) = a1 + b1*exp(-abs(c1)*x); # L0 = 16
f2(x) = a2 + b2*exp(-abs(c2)*x); # L0 = 32
f3(x) = a3 + b3*exp(-abs(c3)*x); # L0 = 64
f4(x) = a4 + b4*exp(-abs(c4)*x); # L0 = 128
f5(x) = a5 + b5*exp(-abs(c5)*x); # L0 = 256

set fit errorvariables
fit [20:*][3.5:4] f1(x) 'gap-400/dtmax-10e-3/L0-16/lmax-10/cmax-10e-3/log'  u 6:24 via a1,b1,c1; # L0 = 16
fit [20:*][3.5:4] f2(x) 'gap-400/dtmax-10e-3/L0-32/lmax-11/cmax-10e-3/log'  u 6:24 via a2,b2,c2; # L0 = 32
fit [20:*][3.5:4] f3(x) 'gap-400/dtmax-10e-3/L0-64/lmax-12/cmax-10e-3/log'  u 6:24 via a3,b3,c3; # L0 = 64
fit [20:*][3.5:4] f4(x) 'gap-400/dtmax-10e-3/L0-128/lmax-13/cmax-10e-3/log' u 6:24 via a4,b4,c4; # L0 = 128
fit [20:*][3.5:4] f5(x) 'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 6:24 via a5,b5,c5; # L0 = 256

set print "gap-400-dtmax-10e-3-cmax-10e-3-domain.dat"
print 16  , a1 , a1_err , 100.*abs (a1 - a5)/a5 , 100.*abs (a1 - 3.73562)/3.73562
print 32  , a2 , a2_err , 100.*abs (a2 - a5)/a5 , 100.*abs (a2 - 3.73562)/3.73562
print 64  , a3 , a3_err , 100.*abs (a3 - a5)/a5 , 100.*abs (a3 - 3.73562)/3.73562
print 128 , a4 , a4_err , 100.*abs (a4 - a5)/a5 , 100.*abs (a4 - 3.73562)/3.73562
print 256 , a5 , a5_err , 0 , 100.*abs (a5 - 3.73562)/3.73562
set print

plot 'gap-400/dtmax-10e-3/L0-16/lmax-10/cmax-10e-3/log'  u 6:25 w l lw 2 lc rgb 'brown'     t 'Brenner, 1961', \
     'gap-400/dtmax-10e-3/L0-16/lmax-10/cmax-10e-3/log'  u 6:27 w l lw 2 lc rgb 'gray'      t 'Cooley et al., 1969', \
     'gap-400/dtmax-10e-3/L0-16/lmax-10/cmax-10e-3/log'  u 6:24 w l lw 1 lc rgb 'black'     t '{\delta}/r=0.4, L0/d=16'.ftitle(a1), \
     'gap-400/dtmax-10e-3/L0-32/lmax-11/cmax-10e-3/log'  u 6:24 w l lw 1 lc rgb 'blue'      t 'L0/d=32'.ftitle(a2), \
     'gap-400/dtmax-10e-3/L0-64/lmax-12/cmax-10e-3/log'  u 6:24 w l lw 1 lc rgb 'red'       t 'L0/d=64'.ftitle(a3), \
     'gap-400/dtmax-10e-3/L0-128/lmax-13/cmax-10e-3/log' u 6:24 w l lw 1 lc rgb 'sea-green' t 'L0/d=128'.ftitle(a4), \
     'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 6:24 w l lw 1 lc rgb 'orange'    t 'L0/d=256'.ftitle(a5)
~~~

~~~gnuplot Evolution of the relative error for the normal force $F_x$ with the size of the domain *L0*
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 8,2,256
set ytics 1.e-8,1.e-1,100
set xlabel 'L_0/d'
set ylabel 'err (%)'
set xrange [8:512]
set yrange [1.e-4:1e1]
set logscale
set title 'delta/r=0.4, dt=1e-2, delta/Delta_{max}=12, cmax=1e-2'
plot 'gap-400-dtmax-10e-3-cmax-10e-3-domain.dat' u 1:4 w lp pt 6 ps 1.25 lw 1 lc rgb "black" t 'Error with reference solution for L_0/d=256', \
     'gap-400-dtmax-10e-3-cmax-10e-3-domain.dat' u 1:5 w lp pt 3 ps 1.25 lw 1 lc rgb "blue"  t 'Error with Brenner, 1961'
~~~

~~~gnuplot Order of convergence of the normal force $F_x$ with the size of the domain *L0*
set ytics -4,1,4
set grid ytics
set xlabel 'L_0/d'
set ylabel 'Order'
set xrange [16:256]
set yrange [0:4]
unset logscale y
set title 'delta/r=0.4, dt=1e-2, delta/Delta_{max}=12, cmax=1e-2'
plot '< sort -k 1,1n gap-400-dtmax-10e-3-cmax-10e-3-domain.dat | awk -f ../data/Brenner1961/order-domain.awk'  u 1:2 w lp ps 1.25 pt 7 lc rgb "black" notitle
~~~
*/

/**
#### Mesh and time convergence study for $\delta/r = 0.4$ and *L0=256*

~~~gnuplot Time evolution of the normal force $F_x$
reset
set terminal svg font ",16"
set key top right font ",8" spacing 0.7
set grid ytics
set xtics 0,5,100
set ytics 0,0.1,10
set xlabel 't/t_{ref}'
set ylabel 'kappa'
set xrange [0:30]
set yrange [3.6:4]
set title 'delta/r=0.4, L0=256'

# Average correction factor kappa

ftitle(a) = sprintf(", kappa_{inf}=%2.4f", a);

f11(x) = a11 + b11*exp(-abs(c11)*x); # L0 = 256, lmax = 12
f12(x) = a12 + b12*exp(-abs(c12)*x); # L0 = 256, lmax = 13
f13(x) = a13 + b13*exp(-abs(c13)*x); # L0 = 256, lmax = 14
f14(x) = a14 + b14*exp(-abs(c14)*x); # L0 = 256, lmax = 15

f21(x) = a21 + b21*exp(-abs(c21)*x); # L0 = 256, lmax = 12
f22(x) = a22 + b22*exp(-abs(c22)*x); # L0 = 256, lmax = 13
f23(x) = a23 + b23*exp(-abs(c23)*x); # L0 = 256, lmax = 14
f24(x) = a24 + b24*exp(-abs(c24)*x); # L0 = 256, lmax = 15

f31(x) = a31 + b31*exp(-abs(c31)*x); # L0 = 256, lmax = 12
f32(x) = a32 + b32*exp(-abs(c32)*x); # L0 = 256, lmax = 13
f33(x) = a33 + b33*exp(-abs(c33)*x); # L0 = 256, lmax = 14
f34(x) = a34 + b34*exp(-abs(c34)*x); # L0 = 256, lmax = 15

f41(x) = a41 + b41*exp(-abs(c41)*x); # L0 = 256, lmax = 12
f42(x) = a42 + b42*exp(-abs(c42)*x); # L0 = 256, lmax = 13
f43(x) = a43 + b43*exp(-abs(c43)*x); # L0 = 256, lmax = 14
f44(x) = a44 + b44*exp(-abs(c44)*x); # L0 = 256, lmax = 15

set fit errorvariables

fit [20:*][3.5:4] f11(x) 'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log' u 6:24 via a11,b11,c11; # L0 = 256, lmax = 12
fit [20:*][3.5:4] f12(x) 'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-10e-3/log' u 6:24 via a12,b12,c12; # L0 = 256, lmax = 13
fit [20:*][3.5:4] f13(x) 'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 6:24 via a13,b13,c13; # L0 = 256, lmax = 14
fit [20:*][3.5:4] f14(x) 'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:24 via a14,b14,c14; # L0 = 256, lmax = 15

fit [20:*][3.5:4] f21(x) 'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-10e-3/log'  u 6:24 via a21,b21,c21; # L0 = 256, lmax = 12
fit [20:*][3.5:4] f22(x) 'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-10e-3/log'  u 6:24 via a22,b22,c22; # L0 = 256, lmax = 13
fit [20:*][3.5:4] f23(x) 'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-10e-3/log'  u 6:24 via a23,b23,c23; # L0 = 256, lmax = 14
fit [1.:*][3.5:*] f24(x) 'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:24 via a24,b24,c24; # L0 = 256, lmax = 15

fit [20:*][3.5:4] f31(x) 'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-1e-3/log' u 6:24 via a31,b31,c31; # L0 = 256, lmax = 12
fit [20:*][3.5:4] f32(x) 'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-1e-3/log' u 6:24 via a32,b32,c32; # L0 = 256, lmax = 13
fit [20:*][3.5:4] f33(x) 'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-1e-3/log' u 6:24 via a33,b33,c33; # L0 = 256, lmax = 14
fit [20:*][3.5:4] f34(x) 'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log' u 6:24 via a34,b34,c34; # L0 = 256, lmax = 15

fit [20:*][3.5:4] f41(x) 'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-1e-3/log'  u 6:24 via a41,b41,c41; # L0 = 256, lmax = 12
fit [20:*][3.5:4] f42(x) 'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-1e-3/log'  u 6:24 via a42,b42,c42; # L0 = 256, lmax = 13
fit [20:*][3.5:*] f43(x) 'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-1e-3/log'  u 6:24 via a43,b43,c43; # L0 = 256, lmax = 14
fit [5.:*][3.5:*] f44(x) 'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:24 via a44,b44,c44; # L0 = 256, lmax = 15

set print "gap-400-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat"
print 1.e-2 , 0.4*0.5/256*2**12 , a11 , a11_err , 100.*abs (a11 - a44)/a44 , 100.*abs (a11 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**13 , a12 , a12_err , 100.*abs (a12 - a44)/a44 , 100.*abs (a12 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**14 , a13 , a13_err , 100.*abs (a13 - a44)/a44 , 100.*abs (a13 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**15 , a14 , a14_err , 100.*abs (a14 - a44)/a44 , 100.*abs (a14 - 3.73562)/3.73562
set print

set print "gap-400-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat"
print 1.e-3 , 0.4*0.5/256*2**12 , a21 , a21_err , 100.*abs (a21 - a44)/a44 , 100.*abs (a21 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**13 , a22 , a22_err , 100.*abs (a22 - a44)/a44 , 100.*abs (a22 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**14 , a23 , a23_err , 100.*abs (a23 - a44)/a44 , 100.*abs (a23 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**15 , a24 , a24_err , 100.*abs (a24 - a44)/a44 , 100.*abs (a24 - 3.73562)/3.73562
set print

set print "gap-400-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat"
print 1.e-2 , 0.4*0.5/256*2**12 , a31 , a31_err , 100.*abs (a31 - a44)/a44 , 100.*abs (a31 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**13 , a32 , a32_err , 100.*abs (a32 - a44)/a44 , 100.*abs (a32 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**14 , a33 , a33_err , 100.*abs (a33 - a44)/a44 , 100.*abs (a33 - 3.73562)/3.73562
print 1.e-2 , 0.4*0.5/256*2**15 , a34 , a34_err , 100.*abs (a34 - a44)/a44 , 100.*abs (a34 - 3.73562)/3.73562
set print

set print "gap-400-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat"
print 1.e-3 , 0.4*0.5/256*2**12 , a41 , a41_err , 100.*abs (a41 - a44)/a44 , 100.*abs (a41 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**13 , a42 , a42_err , 100.*abs (a42 - a44)/a44 , 100.*abs (a42 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**14 , a43 , a43_err , 100.*abs (a43 - a44)/a44 , 100.*abs (a43 - 3.73562)/3.73562
print 1.e-3 , 0.4*0.5/256*2**15 , a44 , a44_err , 0 , 100.*abs (a44 - 3.73562)/3.73562
set print

plot 'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log'  u 6:25 w l lw 2 lc rgb 'brown'     t 'Brenner, 1961', \
     'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log'  u 6:27 w l lw 2 lc rgb 'gray'      t 'Cooley et al., 1969', \
     'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log' u 6:24 w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'black'     t '{\delta}/r=0.4, dt=1e-2, cmax=1e-2, lmax=12'.ftitle(a11), \
     'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-10e-3/log' u 6:24 w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=13'.ftitle(a12), \
     'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 6:24 w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=14'.ftitle(a13), \
     'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:24 w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=15'.ftitle(a14), \
     'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-10e-3/log'  u 6:24 w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t '{\delta}/r=0.4, dt=1e-3, cmax=1e-2, lmax=12'.ftitle(a21), \
     'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-10e-3/log'  u 6:24 w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=13'.ftitle(a22), \
     'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-10e-3/log'  u 6:24 w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=14'.ftitle(a23), \
     'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:24 w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=15'.ftitle(a24), \
     'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-1e-3/log'  u 6:24 w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'red'       t '{\delta}/r=0.4, dt=1e-2, cmax=1e-3, lmax=12'.ftitle(a31), \
     'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-1e-3/log'  u 6:24 w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=13'.ftitle(a32), \
     'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-1e-3/log'  u 6:24 w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=14'.ftitle(a33), \
     'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:24 w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=15'.ftitle(a34), \
     'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-1e-3/log'   u 6:24 w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t '{\delta}/r=0.4, dt=1e-3, cmax=1e-3, lmax=12'.ftitle(a41), \
     'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-1e-3/log'   u 6:24 w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=13'.ftitle(a42), \
     'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-1e-3/log'   u 6:24 w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=14'.ftitle(a43), \
     'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'   u 6:24 w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=15'.ftitle(a44)
~~~u

~~~gnuplot Time evolution of the relative error for the nornal force $F_x$
set ytics 1.e-8,1.e-1,100
set ylabel 'err (%)'
set yrange [1.e-4:1e2]
set logscale y
set title 'delta/r=0.4, L0=256'
plot 'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log' u 6:(100.*abs($27 - $25)/$25) w l lw 2 lc rgb 'brown' t 'Error between Brenner and Cooley', \
     'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'black'     t '{\delta}/r=0.4, dt=1e-2, cmax=1e-2, lmax=12', \
     'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=13', \
     'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=14', \
     'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=15', \
     'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t '{\delta}/r=0.4, dt=1e-3, cmax=1e-2, lmax=12', \
     'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=13', \
     'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=14', \
     'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=15', \
     'gap-400/dtmax-10e-3/L0-256/lmax-12/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'red'       t '{\delta}/r=0.4, dt=1e-2, cmax=1e-3, lmax=12', \
     'gap-400/dtmax-10e-3/L0-256/lmax-13/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=13', \
     'gap-400/dtmax-10e-3/L0-256/lmax-14/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=14', \
     'gap-400/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=15', \
     'gap-400/dtmax-1e-3/L0-256/lmax-12/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t '{\delta}/r=0.4, dt=1e-3, cmax=1e-3, lmax=12', \
     'gap-400/dtmax-1e-3/L0-256/lmax-13/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=13', \
     'gap-400/dtmax-1e-3/L0-256/lmax-14/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=14', \
     'gap-400/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=15'
~~~

~~~gnuplot Evolution of the relative error for the normal force $F_x$ with the reference solution
reset
set terminal svg font ",16"
set key top right font ",14" spacing 1
set grid ytics
set xtics 1.5,2,21
set ytics 1.e-8,1.e-1,100
set xlabel 'delta/Delta_{max}'
set ylabel 'err (%)'
set xrange [1.5:48]
set yrange [1.e-3:1e2]
set logscale
plot 'gap-400-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat' u 2:5 w lp pt 6 ps 1.25 lw 1 lc rgb "black"     t 'Error with reference solution, {\delta}/r=0.4, dt=1e-2, cmax=1e-2', \
     'gap-400-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat' u 2:5 w lp pt 3 ps 1.25 lw 1 lc rgb "blue"      t 'dt=1e-3, cmax=1e-2', \
     'gap-400-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat' u 2:5 w lp pt 8 ps 1.25 lw 1 lc rgb "red"       t 'dt=1e-2, cmax=1e-3', \
     'gap-400-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat' u 2:5 w lp pt 1 ps 1.25 lw 1 lc rgb "sea-green" t 'dt=1e-3, cmax=1e-3'
~~~

~~~gnuplot Evolution of the relative error for the normal force $F_x$ with analytic solution
plot 'gap-400-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat' u 2:6 w lp pt 6 ps 1.25 lw 1 lc rgb "black"     t 'Error with Brenner, 1961, {\delta}/r=0.4, dt=1e-2, cmax=1e-2', \
     'gap-400-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat' u 2:6 w lp pt 3 ps 1.25 lw 1 lc rgb "blue"      t 'dt=1e-3, cmax=1e-2', \
     'gap-400-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat' u 2:6 w lp pt 8 ps 1.25 lw 1 lc rgb "red"       t 'dt=1e-2, cmax=1e-3', \
     'gap-400-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat' u 2:6 w lp pt 1 ps 1.25 lw 1 lc rgb "sea-green" t 'dt=1e-3, cmax=1e-3'
~~~
*/

/**
#### Mesh and time convergence study for $\delta/r = 0.05$ and *L0=256*

~~~gnuplot Time evolution of the normal force $F_x$
reset
set terminal svg font ",16"
set key top right font ",8" spacing 0.7
set grid ytics
set xtics 0,5,100
set ytics 0,0.5,10
set xlabel 't/t_{ref}'
set ylabel 'kappa'
set xrange [0:30]
set yrange [21:23]

# Average correction factor kappa

ftitle(a) = sprintf(", kappa_{inf}=%2.4f", a);

f11(x) = a11 + b11*exp(-abs(c11)*x); # L0 = 256, lmax = 15
f12(x) = a12 + b12*exp(-abs(c12)*x); # L0 = 256, lmax = 16
f13(x) = a13 + b13*exp(-abs(c13)*x); # L0 = 256, lmax = 17
f14(x) = a14 + b14*exp(-abs(c14)*x); # L0 = 256, lmax = 18

f21(x) = a21 + b21*exp(-abs(c21)*x); # L0 = 256, lmax = 15
f22(x) = a22 + b22*exp(-abs(c22)*x); # L0 = 256, lmax = 16
f23(x) = a23 + b23*exp(-abs(c23)*x); # L0 = 256, lmax = 17
f24(x) = a24 + b24*exp(-abs(c24)*x); # L0 = 256, lmax = 18

f31(x) = a31 + b31*exp(-abs(c31)*x); # L0 = 256, lmax = 15
f32(x) = a32 + b32*exp(-abs(c32)*x); # L0 = 256, lmax = 16
f33(x) = a33 + b33*exp(-abs(c33)*x); # L0 = 256, lmax = 17
f34(x) = a34 + b34*exp(-abs(c34)*x); # L0 = 256, lmax = 18

f41(x) = a41 + b41*exp(-abs(c41)*x); # L0 = 256, lmax = 15
f42(x) = a42 + b42*exp(-abs(c42)*x); # L0 = 256, lmax = 16
f43(x) = a43 + b43*exp(-abs(c43)*x); # L0 = 256, lmax = 17
f44(x) = a44 + b44*exp(-abs(c44)*x); # L0 = 256, lmax = 18

set fit errorvariables

fit [20:*][21:*] f11(x) 'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:24 via a11,b11,c11; # L0 = 256, lmax = 15
fit [20:*][21:*] f12(x) 'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-10e-3/log' u 6:24 via a12,b12,c12; # L0 = 256, lmax = 16
fit [20:*][21:*] f13(x) 'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-10e-3/log' u 6:24 via a13,b13,c13; # L0 = 256, lmax = 17
fit [2.:*][21:*] f14(x) 'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-10e-3/log' u 6:24 via a14,b14,c14; # L0 = 256, lmax = 18

fit [20:*][21:*] f21(x) 'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:24 via a21,b21,c21; # L0 = 256, lmax = 15
fit [10:*][21:*] f22(x) 'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-10e-3/log'  u 6:24 via a22,b22,c22; # L0 = 256, lmax = 16
fit [1.:*][21:*] f23(x) 'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-10e-3/log'  u 6:24 via a23,b23,c23; # L0 = 256, lmax = 17
fit [0.1:*][3.5:*] f24(x) 'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-10e-3/log'  u 6:24 via a24,b24,c24; # L0 = 256, lmax = 18

fit [20:*][21:*] f31(x) 'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log' u 6:24 via a31,b31,c31; # L0 = 256, lmax = 15
fit [20:*][21:*] f32(x) 'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-1e-3/log' u 6:24 via a32,b32,c32; # L0 = 256, lmax = 16
fit [5.:*][21:*] f33(x) 'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-1e-3/log' u 6:24 via a33,b33,c33; # L0 = 256, lmax = 17
fit [0.1:*][21:*] f34(x) 'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-1e-3/log' u 6:24 via a34,b34,c34; # L0 = 256, lmax = 18

fit [10:*][21:*] f41(x) 'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:24 via a41,b41,c41; # L0 = 256, lmax = 15
fit [1.:*][21:*] f42(x) 'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-1e-3/log'  u 6:24 via a42,b42,c42; # L0 = 256, lmax = 16
fit [0.1:*][21:*] f43(x) 'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-1e-3/log'  u 6:24 via a43,b43,c43; # L0 = 256, lmax = 17
#fit [0.1:*][21:*] f44(x) 'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-1e-3/log'  u 6:24 via a44,b44,c44; # L0 = 256, lmax = 18

set print "gap-50-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat"
print 1.e-2 , 0.05*0.5/256*2**15 , a11 , a11_err , 100.*abs (a11 - a34)/a34 , 100.*abs (a11 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**16 , a12 , a12_err , 100.*abs (a12 - a34)/a34 , 100.*abs (a12 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**17 , a13 , a13_err , 100.*abs (a13 - a34)/a34 , 100.*abs (a13 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**18 , a14 , a14_err , 100.*abs (a14 - a34)/a34 , 100.*abs (a14 - 21.5858)/21.5858
set print

set print "gap-50-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat"
print 1.e-3 , 0.05*0.5/256*2**15 , a21 , a21_err , 100.*abs (a21 - a34)/a34 , 100.*abs (a21 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**16 , a22 , a22_err , 100.*abs (a22 - a34)/a34 , 100.*abs (a22 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**17 , a23 , a23_err , 100.*abs (a23 - a34)/a34 , 100.*abs (a23 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**18 , a24 , a24_err , 100.*abs (a24 - a34)/a34 , 100.*abs (a24 - 21.5858)/21.5858
set print

set print "gap-50-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat"
print 1.e-2 , 0.05*0.5/256*2**15 , a31 , a31_err , 100.*abs (a31 - a34)/a34 , 100.*abs (a31 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**16 , a32 , a32_err , 100.*abs (a32 - a34)/a34 , 100.*abs (a32 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**17 , a33 , a33_err , 100.*abs (a33 - a34)/a34 , 100.*abs (a33 - 21.5858)/21.5858
print 1.e-2 , 0.05*0.5/256*2**18 , a34 , a34_err , 100.*abs (a34 - a34)/a34 , 100.*abs (a34 - 21.5858)/21.5858
set print

set print "gap-50-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat"
print 1.e-3 , 0.05*0.5/256*2**15 , a41 , a41_err , 100.*abs (a41 - a34)/a34 , 100.*abs (a41 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**16 , a42 , a42_err , 100.*abs (a42 - a34)/a34 , 100.*abs (a42 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**17 , a43 , a43_err , 100.*abs (a43 - a34)/a34 , 100.*abs (a43 - 21.5858)/21.5858
print 1.e-3 , 0.05*0.5/256*2**18 , a44 , a44_err , 0 , 100.*abs (a44 - 21.5858)/21.5858
set print

plot 'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:25 w l lw 2 lc rgb 'brown'     t 'Brenner, 1961', \
     'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:27 w l lw 2 lc rgb 'gray'      t 'Cooley et al., 1969', \
     'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:24 w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'black'     t '{\delta}/r=0.05, dt=1e-2, cmax=1e-2, lmax=15'.ftitle(a11), \
     'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-10e-3/log' u 6:24 w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=16'.ftitle(a12), \
     'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-10e-3/log' u 6:24 w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=17'.ftitle(a13), \
     'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-10e-3/log' u 6:24 w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=18'.ftitle(a14), \
     'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:24 w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t '{\delta}/r=0.05, dt=1e-3, cmax=1e-2, lmax=15'.ftitle(a21), \
     'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-10e-3/log'  u 6:24 w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=16'.ftitle(a22), \
     'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-10e-3/log'  u 6:24 w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=17'.ftitle(a23), \
     'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-10e-3/log'  u 6:24 w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=18'.ftitle(a24), \
     'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:24 w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'red'       t '{\delta}/r=0.05, dt=1e-2, cmax=1e-3, lmax=15'.ftitle(a31), \
     'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-1e-3/log'  u 6:24 w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=16'.ftitle(a32), \
     'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-1e-3/log'  u 6:24 w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=17'.ftitle(a33), \
     'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-1e-3/log'  u 6:24 w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=18'.ftitle(a34), \
     'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'   u 6:24 w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t '{\delta}/r=0.05, dt=1e-3, cmax=1e-3, lmax=15'.ftitle(a41), \
     'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-1e-3/log'   u 6:24 w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=16'.ftitle(a42), \
     'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-1e-3/log'   u 6:24 w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=17'.ftitle(a43), \
     'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-1e-3/log'   u 6:24 w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=18'.ftitle(a44)
~~~

~~~gnuplot Time evolution of the relative error for the nornal force $F_x$
set ytics 1.e-8,1.e-1,100
set ylabel 'err (%)'
set yrange [1.e-4:1e2]
set logscale y
plot 'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:(100.*abs($27 - $25)/$25) w l lw 2 lc rgb 'brown' t 'Error between Brenner and Cooley', \
     'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'black'     t '{\delta}/r=0.05, dt=1e-2, cmax=1e-2, lmax=15', \
     'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=16', \
     'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=17', \
     'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-10e-3/log' u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'black'     t 'lmax=18', \
     'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t '{\delta}/r=0.05, dt=1e-3, cmax=1e-2, lmax=15', \
     'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=16', \
     'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=17', \
     'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-10e-3/log'  u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'blue'      t 'lmax=18', \
     'gap-50/dtmax-10e-3/L0-256/lmax-15/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 400  lw 1 lc rgb 'red'       t '{\delta}/r=0.05, dt=1e-2, cmax=1e-3, lmax=15', \
     'gap-50/dtmax-10e-3/L0-256/lmax-16/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=16', \
     'gap-50/dtmax-10e-3/L0-256/lmax-17/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=17', \
     'gap-50/dtmax-10e-3/L0-256/lmax-18/cmax-1e-3/log'  u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 400  lw 1 lc rgb 'red'       t 'lmax=18', \
     'gap-50/dtmax-1e-3/L0-256/lmax-15/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 6 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t '{\delta}/r=0.05, dt=1e-3, cmax=1e-3, lmax=15', \
     'gap-50/dtmax-1e-3/L0-256/lmax-16/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 3 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=16', \
     'gap-50/dtmax-1e-3/L0-256/lmax-17/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 8 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=17', \
     'gap-50/dtmax-1e-3/L0-256/lmax-18/cmax-1e-3/log'   u 6:(100.*$26/$25) w lp pt 1 ps 1.25 pi 4000 lw 1 lc rgb 'sea-green' t 'lmax=18'
~~~

~~~gnuplot Evolution of the relative error for the normal force $F_x$ with the reference solution
reset
set terminal svg font ",16"
set key top right font ",14" spacing 1
set grid ytics
set xtics 1.5,2,21
set ytics 1.e-8,1.e-1,100
set xlabel 'delta/Delta_{max}'
set ylabel 'err (%)'
set xrange [1.5:48]
set yrange [1.e-2:1e2]
set logscale
plot 'gap-50-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat' u 2:5 w lp pt 6 ps 1.25 lw 1 lc rgb "black"     t 'Error with reference solution, {\delta}/r=0.05, dt=1e-2, cmax=1e-2', \
     'gap-50-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat' u 2:5 w lp pt 3 ps 1.25 lw 1 lc rgb "blue"      t 'dt=1e-3, cmax=1e-2', \
     'gap-50-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat' u 2:5 w lp pt 8 ps 1.25 lw 1 lc rgb "red"       t 'dt=1e-2, cmax=1e-3', \
     'gap-50-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat' u 2:5 w lp pt 1 ps 1.25 lw 1 lc rgb "sea-green" t 'dt=1e-3, cmax=1e-3'
~~~

~~~gnuplot Evolution of the relative error for the normal force $F_x$ with analytic solution
plot 'gap-50-dtmax-1e-2-L0-256-cmax-1e-2-convergence.dat' u 2:6 w lp pt 6 ps 1.25 lw 1 lc rgb "black"     t 'Error with Brenner, 1961, {\delta}/r=0.05, dt=1e-2, cmax=1e-2', \
     'gap-50-dtmax-1e-3-L0-256-cmax-1e-2-convergence.dat' u 2:6 w lp pt 3 ps 1.25 lw 1 lc rgb "blue"      t 'dt=1e-3, cmax=1e-2', \
     'gap-50-dtmax-1e-2-L0-256-cmax-1e-3-convergence.dat' u 2:6 w lp pt 8 ps 1.25 lw 1 lc rgb "red"       t 'dt=1e-2, cmax=1e-3', \
     'gap-50-dtmax-1e-3-L0-256-cmax-1e-3-convergence.dat' u 2:6 w lp pt 1 ps 1.25 lw 1 lc rgb "sea-green" t 'dt=1e-3, cmax=1e-3'
~~~
*/

/**
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

@article{Cooley1969,
  title={On the slow motion generated in a viscous fluid by the approach of a sphere to a plane wall or stationary sphere},
  author={Cooley, MDA and O'neill, ME},
  journal={Mathematika},
  volume={16},
  number={1},
  pages={37--49},
  year={1969},
  publisher={London Mathematical Society}
}
~~~
*/
