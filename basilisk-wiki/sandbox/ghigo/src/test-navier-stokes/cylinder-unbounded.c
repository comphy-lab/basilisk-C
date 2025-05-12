/**
# Flow past a fixed cylinder at different Reynolds number

We reproduce here the test case propsed in [Wu and Shu](#wu2009) and
solve the Navier-Stokes equations and add the cylinder using an
[embedded boundary](/src/embed.h). */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define uref (0.1) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u

#if RE == 1 // Re = 40
#define Re (40.)
#elif RE == 2 // Re = 100
#define Re (100.)
#else // Re = 20
#define Re (20.) // Particle Reynolds number Re = ud/nu
#endif // RE

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define EPS (1.e-12)
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
be taken into account. */

face vector muv[];

/**
We also define a reference velocity field. */

scalar un[];

/**
We define the mesh adaptation parameters. */

#define lmin (7) // Min mesh refinement level (l=7 is 3pt/d)
#define lmax (11) // Max mesh refinement level (l=11 is 50pt/d)
#define cmax (1.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
#if RE == 2 // Re = 100
  /**
  The domain is $50d\times 50d$. */

  L0 = 50.*(d);
#else // Re = 20 and 40
  /**
  The domain is $40d\times 40d$. */

  L0 = 40.*(d);
#endif // RE
  size (L0);
  origin (0., -L0/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);
  NITERMAX     = 150;

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

  /**
  We first define the cylinder's position, following [Wu and
  Shu](#wu2009). */

#if RE == 2 // Re = 100
  p_p.x = 20.*(d) + (EPS);
#else // Re = 20 and 40
  p_p.x = 16.*(d) + (EPS);
#endif // RE
  p_p.y = (EPS);
  
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
  We define the no-slip boundary conditions for the velocity. */
   
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
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We also reset the embedded fractions to avoid interpolation errors
  on the geometry. This would affect in particular the computation of
  the pressure contribution to the hydrodynamic forces. */

  p_shape (cs, fs, p_p);
}
#endif // TREE

/**
## Outputs */

double CDm1 = 0.;
double sdt = 0., CDavg = 0., CLavg = 0.;
double CDmin = HUGE, CDmax = -HUGE, CLmin = HUGE, CLmax = -HUGE;
double La = 0., Lamin = HUGE, Lamax = -HUGE;

event logfile (i++; t <= 200.*(tref))
{
  double du = change (u.x, un);
  
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(d));

  double dF = fabs (CDm1 - CD);
  CDm1 = CD;  

  fprintf (stderr, "%g %d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
	   (Re),
	   i, t/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CL,
	   du, dF);
  fflush (stderr);

  if ((du < 1.e-5 && dF < 1.e-4 && t > 100.*(tref)) || t > 170.*(tref)) {
    
    /**
    We compute the average, minimun and maximum drag and lift
    coefficients after an initial transitory state. */

    CDavg += CD*dt;
    CLavg += CL*dt;
    sdt   += dt;
    
    if (CD > CDmax)
      CDmax = CD;
    if (CD < CDmin)
      CDmin = CD;
    if (CL > CLmax)
      CLmax = CL;
    if (CL < CLmin)
      CLmin = CL;

    /**
    #### Recirculation length.

    We evaluate the lenght of the recirculation area downstream of the
    cylinder. */

    double step = ((L0)/(1 << (lmax)));
    
    La = -((p_p.x) + (d)/2);
    for (double x = ((p_p.x) + (d)/2) + step/2.; x <= (L0) - step/2.; x += step) {
      double o = interpolate (u.x, x, (p_p.y));
      if (o > 0.) {
	La += x;
	break;
      }
    }

    /**
    We scale the recirculation area by *d/2*. */

    La /= (d)/2.;

    if (La > Lamax)
      Lamax = La;
    if (La < Lamin)
      Lamin = La;
  }
}

event coeffs (t = end)
{
  char name1[80];
  sprintf (name1, "avg-min-max.dat");
  static FILE * fp = fopen (name1, "w");
  fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g\n",
	   (Re), (1 << (lmax)),
	   CDavg, CDmin, CDmax,
	   CLavg, CLmin, CLmax,
	   La, Lamin, Lamax);
  fflush (fp);
  fclose (fp);
}

/**
## Snapshots */

event snapshot (t = end)
{
  scalar omega[];
  vorticity (u, omega);
  
  char name2[80];

  /**
  We first plot the entire domain. */
 
  view (fov = 20, camera = "front",
	tx = -(L0)/2./(L0), ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh.png");
  save (name2);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux.png");
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy.png");
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name2, "p.png");
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega.png");
  save (name2);

  /**
  We then zoom on the cylinder. */

  view (fov = 6, camera = "front",
	tx = -(p_p.x + 3.)/(L0), ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh-zoom.png");
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux-zoom.png");
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy-zoom.png");
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("p", map = cool_warm);
  sprintf (name2, "p-zoom.png");
  save (name2);

  draw_vof ("cs", "fs", lw = 5);
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega-zoom.png");
  save (name2);
}

/**
## Results

#### Snapshots for level 8

<p>

<img src="cylinder-unbounded/mesh.png" alt="drawing" width="25%"/>
<img src="cylinder-unbounded/mesh-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-unbounded/ux.png" alt="drawing" width="25%"/>
<img src="cylinder-unbounded/ux-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-unbounded/uy.png" alt="drawing" width="25%"/>
<img src="cylinder-unbounded/uy-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-unbounded/p.png" alt="drawing" width="25%"/>
<img src="cylinder-unbounded/p-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-unbounded/omega.png" alt="drawing" width="25%"/>
<img src="cylinder-unbounded/omega-zoom.png" alt="drawing" width="25%"/> 

</p>

#### Drag and lift coefficients

We compare the steady values of the drag coefficient $C_D$ with those
obtained in the benchmark study of [Wu et al., 2009](#wu2009).

~~~gnuplot Time evolution of the drag coefficient $C_D$
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,20,200
set xlabel 't/(d/u)'
set ylabel 'C_{D}'
set xrange [0:200]
set yrange [1:3]
plot 2.091 w p pt 1 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=20", \
     1.565 w p pt 2 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=40", \
     1.364 w p pt 3 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=100",	\
     'log' u 3:15 w l lc rgb "blue" t "Basilisk, l=9"
~~~

~~~gnuplot Time evolution of the lift coefficient $C_L$
set ylabel 'C_{L}'
set yrange [-0.4:0.8]
plot 0     w p pt 1 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=20/40", \
     0.344 w p pt 2 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=100", \
     'log' u 3:16 w l lc rgb "blue" t "Basilisk, l=9"
~~~

~~~gnuplot Minimum and maximum drag coeffficient $C_D$
set xtics 128,2,8192
set xlabel 'N'
set ylabel 'C_{D}'
set xrange [128:8192]
set yrange [1:3]
set logscale x
plot 2.091 w lp pt 1 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=20", \
     1.565 w lp pt 2 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=40", \
     1.364 w lp pt 3 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=100",	\
     'avg-min-max.dat' u 2:3 w p pt 7  lc rgb "black" t "Basilisk, C_{D,avg}", \
     ''                u 2:4 w p pt 5  lc rgb "blue"  t "Basilisk, C_{D,min}", \
     ''                u 2:5 w p pt 2  lc rgb "red"   t "Basilisk, C_{D,max}"
~~~

~~~gnuplot Minimum and maximum lift coeffficient $C_L$
set ylabel 'C_{L}'
set yrange [-0.4:0.8]
plot 0     w lp pt 1 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=20/40", \
     0.344 w lp pt 2 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=100", \
     'avg-min-max.dat' u 2:6 w p pt 7  lc rgb "black" t "Basilisk, C_{L,avg}", \
     ''                u 2:7 w p pt 5  lc rgb "blue"  t "Basilisk, C_{L,min}", \
     ''                u 2:8 w p pt 2  lc rgb "red"   t "Basilisk, C_{L,max}"
~~~

#### Recirculation area

We do the same for the length of the recirculation area $L_a$.

~~~gnuplot Recirculation area
set key top right
set ylabel 'La'
set yrange [1:6]
plot 1.86 w lp pt 1 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=20", \
     4.62 w lp pt 2 ps 0.6 lc rgb "black" t "Wu et al., 2009, Re=40", \
     'avg-min-max.dat' u 2:9  w p pt 7  lc rgb "black" t "Basilisk, La_{avg}", \
     ''                u 2:10 w p pt 5  lc rgb "blue"  t "Basilisk, La_{min}", \
     ''                u 2:11 w p pt 2  lc rgb "red"   t "Basilisk, La_{max}"
~~~

## References

~~~bib
@article{williamson1988,
  title={Defining a universal and continuous Strouhal–Reynolds number relationship for the laminar vortex shedding of a circular cylinder},
  author={Williamson, C.H.K.},
  journal={The Physics of Fluid},
  volume={31},
  pages={2742--2744},
  year={1988},
  publisher={American Institute of Physics}
}

@article{wu2009,
  title={Implicit velocity correction-based immersed boundary-lattice Boltzmann method and its applications},
  author={Wu, J. and Shu, C.},
  journal={Journal of Computational Physics},
  volume={228},
  pages={1963--1979},
  year={2009},
  publisher={Elsevier}
}
~~~
*/
