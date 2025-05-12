/**
# Flow past a fixed confined cylinder at $Re = 20$

We reproduce here a benchmark solution proposed in [Schafer and
Turek](#schafer1996). In this benchmark case, we study the flow around
a cylinder positioned slightly off mid-stream. Instead of using a
square domain, we use embedded boundaries to define a rectangular
domain.

We solve here the Navier-Stokes equations and add the cylinder using
an [embedded boundary](/src/embed.h). */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "view.h"

/**
## Reference solution */

#define d    (0.1) // Cylinder characteristic length
#if RE // Re = 100
#define uref (2./3.*1.5) // Reference velocity, uref
#define Re   (100.)
#else // Re = 20
#define uref (2./3.*0.3) // Reference velocity, uref
#define Re   (20.)
#endif // RE
#define tref ((d)/(uref)) // Reference time, tref=d/u
#define um   (1.5*(uref)) // Maximum velocity

/**
We also define the shape of the domain. */

#define ht  (2.1*(d)) // Top-half width of the channel
#define hb  (2.*(d)) // Bottom-half width of the channel
#define EPS (1.e-12)

#define cylinder(x,y) (sq ((x)) + sq ((y)) - sq ((d)/2.))
#define wall(y,w)     ((y) - (w)) // + over, - under

coord p_p;

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = intersection (
			  (cylinder ((x - p.x), (y - p.y))),
			  intersection (
					-(wall (y,  0.5*((hb) + (ht)))),
					 (wall (y, -0.5*((hb) + (ht))))
					)
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
Since both the walls of the channel and the cylinder are described
using the same embedded volume fraction *cs*, we declare a color field
*p_col* that we use to color the cylinder only. */

scalar p_col[];

void p_shape_col (scalar c, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (cylinder ((x - p.x), (y - p.y)));
  boundary ({phi});
  fractions (phi, c);
}

/**
We define the mesh adaptation parameters. */

#define lmin (6)  // Min mesh refinement level (l=6 is 3pt/d)
#define lmax (11) // Max mesh refinement level (l=11 is 92pt/d)
#define cmax (5.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $22 d\times 22 d$. */

  L0 = 22.*(d);
  size (L0);
  origin (0., -L0/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);

  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions with a parabolic velocity profile. */

#define profile(y) ((y) <= 0.5*((hb) + (ht)) && (y) >= -0.5*((hb) + (ht)) ? \
		    4.*((y) + 0.5*((hb) + (ht)))/((hb) + (ht))*(1. - ((y) + 0.5*((hb) + (ht)))/((hb) + (ht))) : \
		    0.)

u.n[left] = dirichlet ((um)*(profile (y)));
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

uf.n[left]   = (um)*(profile (y));
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
  We first define the cylinder's position. */

  p_p.x = 2.*(d) + (EPS);
  p_p.y = (hb) - 0.5*((hb) + (ht)) + (EPS);
  
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
    u.x[] = cs[]*(um)*(profile (y));
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

double sdt = 0., CDm1 = 0.;
double CDavg = 0., CDmin = HUGE, CDmax = -HUGE;
double CLavg = 0., CLmin = HUGE, CLmax = -HUGE;
double La = 0., Lamin = HUGE, Lamax = -HUGE;
double DP = 0., DPmin = HUGE, DPmax = -HUGE;

event logfile (i++; t <= 200.*(tref))
{
  double du = change (u.x, un);
  
  coord Fp, Fmu;
  p_shape_col (p_col, p_p);
  boundary ({p_col});
  embed_color_force (p, u, mu, p_col, &Fp, &Fmu);

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
    We compute here the pressure drop between the front and back of
    the cylinder end evaluate the lenght of the recirculation area. */

    double step = ((L0)/(1 << (lmax)));

    /**
    #### Pressure drop. */

    double pa = interpolate (p, ((p_p.x) - (d)/2) - step/2., (p_p.y));
    double pe = interpolate (p, ((p_p.x) + (d)/2) + step/2., (p_p.y));
    DP = pa - pe;

    if (DP > DPmax)
      DPmax = DP;
    if (DP < DPmin)
      DPmin = DP;
    
    /**
    #### Recirculation length. */

    La = -((p_p.x) + (d)/2);
    for (double x = ((p_p.x) + (d)/2) + step/2.; x <= (L0) - step/2.; x += step) {
      double o = interpolate (u.x, x, (p_p.y));
      if (o > 0.) {
	La += x;
	break;
      }
    }

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
  fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   (Re), (1 << (lmax)),
	   CDavg, CDmin, CDmax,
	   CLavg, CLmin, CLmax,
	   La, Lamin, Lamax,
	   DP, DPmin, DPmax);
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

  view (fov = 3, camera = "front",
	tx = -(p_p.x + 4.*(d))/(L0), ty = 0.,
	bg = {1,1,1},
	width = 800, height = 250);

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

#### Snapshots for level 11

<p>

<img src="cylinder-confined/mesh.png" alt="drawing" width="25%"/>
<img src="cylinder-confined/mesh-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-confined/ux.png" alt="drawing" width="25%"/>
<img src="cylinder-confined/ux-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-confined/uy.png" alt="drawing" width="25%"/>
<img src="cylinder-confined/uy-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-confined/p.png" alt="drawing" width="25%"/>
<img src="cylinder-confined/p-zoom.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-confined/omega.png" alt="drawing" width="25%"/>
<img src="cylinder-confined/omega-zoom.png" alt="drawing" width="25%"/> 

</p>

#### Drag and lift coefficients

We compare the steady values of the drag and lift coefficients $C_D$
and $C_L$ with those obtained in the benchmark study of [Schafer and
Turek](#schafer1996).

~~~gnuplot Time evolution of the drag coefficient $C_D$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,20,200
set xlabel 't/(d/u)'
set ylabel 'C_{D}'
set xrange [0:200]
set yrange [0:11]
plot 5.57 w p pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{D,min}", \
     5.59 w p pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{D,max}", \
     3.22 w p pt 3 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{D,min}", \
     3.24 w p pt 4 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{D,max}", \
     'log' u 3:15 w l lc rgb "blue" t "Basilisk, l=11"
~~~

~~~gnuplot Time evolution of the lift coefficient $C_L$
set ylabel 'C_{L}'
set yrange [-2:4]
plot 0.0104 w p pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{L,min}", \
     0.0110 w p pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{L,max}", \
     0.99   w p pt 3 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{L,min}", \
     1.01   w p pt 4 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{L,max}", \
     'log' u 3:16 w l lc rgb "blue" t "Basilisk, l=11"
~~~

~~~gnuplot Minimum and maximum drag coeffficient $C_D$
set xtics 64,2,4096
set xlabel 'N'
set ylabel 'C_{D}'
set xrange [64:4096]
set yrange [1:9]
set logscale x
plot 5.57 w lp pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{D,min}", \
     5.59 w lp pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{D,max}", \
     3.22 w lp pt 3 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{D,min}", \
     3.24 w lp pt 4 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{D,max}", \
     'avg-min-max.dat' u 2:3 w p pt 7  lc rgb "black" t "Basilisk, C_{D,avg}", \
     ''                u 2:4 w p pt 5  lc rgb "blue"  t "Basilisk, C_{D,min}", \
     ''                u 2:5 w p pt 2  lc rgb "red"   t "Basilisk, C_{D,max}"
~~~

~~~gnuplot Minimum and maximum lift coeffficient $C_L$
set ylabel 'C_{L}'
set yrange [-2:4]
plot 0.0104 w lp pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{L,min}", \
     0.0110 w lp pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, C_{L,max}", \
     0.99   w lp pt 3 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{L,min}", \
     1.01   w lp pt 4 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, C_{L,max}", \
     'avg-min-max.dat' u 2:6 w p pt 7  lc rgb "black" t "Basilisk, C_{L,avg}", \
     ''                u 2:7 w p pt 5  lc rgb "blue"  t "Basilisk, C_{L,min}", \
     ''                u 2:8 w p pt 2  lc rgb "red"   t "Basilisk, C_{L,max}"
~~~

#### Pressure drop and recirculation area

We do the same for the pressure drop $\Delta p$ and the length of the
recirculation area $L_a$.

~~~gnuplot Recirculation area
set ylabel 'La'
set yrange [0:0.25]
plot 0.0842 w lp pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, La_{min}", \
     0.0852 w lp pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, La_{max}", \
     'avg-min-max.dat' u 2:9  w p pt 7  lc rgb "black" t "Basilisk, La_{avg}", \
     ''                u 2:10 w p pt 5  lc rgb "blue"  t "Basilisk, La_{min}", \
     ''                u 2:11 w p pt 2  lc rgb "red"   t "Basilisk, La_{max}"
~~~

~~~gnuplot Pressure drop
set ylabel '{D}p'
set yrange [-1:5]
plot 0.1172 w lp pt 1 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, {D}p_{min}", \
     0.1176 w lp pt 2 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=20, {D}p_{max}", \
     2.46   w lp pt 3 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, {D}p_{min}", \
     2.5    w lp pt 4 ps 0.6 lc rgb "black" t "Schafer et al., 1996, Re=100, {D}p_{max}", \
     'avg-min-max.dat' u 2:12 w p pt 7  lc rgb "black" t "Basilisk, DP_{avg}", \
     ''                u 2:13 w p pt 5  lc rgb "blue"  t "Basilisk, DP_{min}", \
     ''                u 2:14 w p pt 2  lc rgb "red"   t "Basilisk, DP_{max}"
~~~

## References

~~~bib
@article{schafer1996,
  title={Benchmark Computations of Laminar Flow Around a Cylinder},
  author={Schafer, M. and Turek., S.},
  journal={Flow Simulation with High-Performance Computers II},
  year={1996},
  publisher={}
}
~~~
*/
