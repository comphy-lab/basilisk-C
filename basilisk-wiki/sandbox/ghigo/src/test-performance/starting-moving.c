/**
# Starting flow around a moving cylinder */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "../myperfs.h"

/**
## Reference solution */

#define d    (1.)
#define uref (1.) // Reference velocity, uref
#define tref ((d)/(uref)) // Reference time, tref=d/u
#define Re   (1000.)

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq (x) + sq (y) - sq ((d)/2.))

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
Finally, we define the mesh adaptation parameters and vary the maximum
level of refinement.

High-resolution is needed to resolve the boundary layers
properly. [Mohaghegh et al., 2017](#mohaghegh2017) propose to use a
maximum resolution of order $D/10/\sqrt{Re}$, with $Re$ the Reynolds
number and $D$ the cylinder diameter. The cylinder diameter will be
set to unity, and the domain size to 18, so that the corresponding
levels of refinement are approximately 12 and 16 for $Re=1000$ and
$Re=9500$, respectively. */

#define lmin (6)  // Min mesh refinement level (l=6 is 3pt/d)
#define lmax (11) // Max mesh refinement level (l=12 is 227pt/d)
#define cmax (1.e-3*(uref)) // Absolute refinement criteria for the velocity field

/**
We finally define a useful function that allows us to define a
rectangular block mesh refinement (BMR) around the cylinder. */

#define rectangle(x,hx,y,hy) (union (					\
				     union ((x) - (hx)/2., -(x) - (hx)/2.), \
				     union ((y) - (hy)/2., -(y) - (hy)/2.)) \
			      )

int main ()
{  
  /**
  The domain is $18\times 18$ and only half the cylinder is
  represented. */

  L0 = 18.;
  size (L0);
  origin (-L0/2., 0.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4; // To reduce oscillations due to small timestep
  TOLERANCE_MU = 1.e-4*(uref);
  
  /**
  We initialize the grid. */
  
#if UNIFORM
  N = 1 << (lmax);
#else // ADAPT
  N = 1 << (lmin);
#endif // UNIFORM
  init_grid (N);
  
  run();
}

/**
## Boundary conditions */

u.n[left] = dirichlet (0);
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

uf.n[left]   = 0;
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
  We initialize the embedded boundary. */
  
#if TREE && !UNIFORM
#if BMR
  int lvl = (lmin);
  double wxmin = 8.*(d), wxmax = 5.*(d);
  double xmin = -(wxmin/2. - 2.*(d)), xmax = -(wxmax/2. - 1.*(d));
  double wymin = 8.*(d), wymax = 2.*(d);
  double ymin = 0., ymax = 0.;
  while (lvl < (lmax)) {
    refine (
	    (rectangle ((x - (p_p.x + (xmin + (xmax - xmin)/((lmax) - (lmin))*((lvl + 1) - (lmin))))),
			wxmin + (wxmax - wxmin)/((lmax) - (lmin))*((lvl + 1) - (lmin)),
			(y - (p_p.y + (ymin + (ymax - ymin)/((lmax) - (lmin))*((lvl + 1) - (lmin))))),
			wymin + (wymax - wymin)/((lmax) - (lmin))*((lvl + 1) - (lmin)))) <= 0. &&
	    level < (lvl + 1)
	    );
    p_shape (cs, fs, p_p);
    lvl ++;
  }
  
  /**
  In the case of the moving cylinder, we don't unrefine inside the
  cylinder. */
  
  unrefine ((rectangle ((x - (p_p.x + xmin)), wxmin,
			(y - (p_p.y + ymin)), wymin)) > 0. &&
	    level > 1);
#else // AMR
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
#endif // TREE && !UNIFORM
  
  p_shape (cs, fs, p_p);
  
  /**
  We initialize the particle's speed. */
  
  p_u.x  = -(uref);
}

/**
## Embedded boundaries 

The cylinder's position is advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  p_p.x -= (uref)*(dt);
}

/**
## Adaptive mesh refinement */

#if TREE && !BMR && !UNIFORM
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
#endif // TREE && !BMR && !UNIFORM

/**
## Profiling */

#if TRACE > 1
event profiling (i += 25) {
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

/**
## Outputs

#### Positions of the separation points

We would like to track the positions with time of the separation
points on the surface of the cylinder, as done by [K & L,
1995](#koumoutsakos1995).

We first define a function which computes the vorticity at the surface
of the cylinder and returns an array of $(\theta,\omega)$ pairs, with
$\theta$ the angular coordinate and $\omega$ the corresponding value
of vorticity. */

typedef struct {
  double theta, omega;
} ThetaOmega;

int compar_theta (const void * a, const void * b)
{
  const ThetaOmega * p1 = a, * p2 = b;
  return p1->theta > p2->theta ? 1 : -1;
}

ThetaOmega * theta_omega()
{
  Array * a = array_new();
  foreach(serial)
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      embed_geometry (point, &b, &n);
      double xe = x + b.x*Delta, ye = y + b.y*Delta;

      ThetaOmega t;
      t.omega = embed_vorticity (point, u, b, n);
      t.theta = atan2(ye - p_p.y, xe - p_p.x);
      array_append (a, &t, sizeof (ThetaOmega));
    }
  qsort (a->p, a->len/sizeof(ThetaOmega), sizeof(ThetaOmega), compar_theta);
  ThetaOmega t = {nodata, nodata};
  array_append (a, &t, sizeof (ThetaOmega));
  ThetaOmega * p = a->p;
  free (a);
  return p;
}

/**
The zeros of the function approximated by the $(\theta,\omega)$ array
are then recorded, together with the corresponding time, in the file
pointed to by *fp*. */

void omega_zero (FILE * fp)
{
  // fixme: this function will not work with MPI
  ThetaOmega * a = theta_omega();
  for (ThetaOmega * o = a; (o + 1)->theta != nodata; o++) {
    ThetaOmega * o1 = o + 1;
    if (o1->omega*o->omega < 0.)
      fprintf (fp, "%g %g\n", t,
	       o->theta + o->omega*(o1->theta - o->theta)/
	       (o->omega - o1->omega));
  }
  free (a);
  fflush (fp);
}

/**
#### Log files */

event logfile (i++)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
	   i, t, dt,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   Fp.x, Fp.y, Fmu.x, Fmu.y);
  fflush (stderr);

  static FILE * fo = fopen ("omega.dat", "w");
#if !_MPI
  omega_zero (fo);
#endif // !_MPI
}

/**
#### Surface vorticity profiles

We are also interested in the details on the surface of the cylinder,
in particular surface vorticity. */

void cpout (FILE * fp)
{
  foreach(serial)
    if (cs[] > 0. && cs[] < 1.) {

      coord b, n;
      double area = embed_geometry (point, &b, &n);
      double xe = x + b.x*Delta, ye = y + b.y*Delta;

      fprintf (fp, "%g %g %g %g %g %g\n",
	       xe - p_p.x, // 1
	       ye - p_p.y, // 2
	       atan2(ye - p_p.y, xe - p_p.x), // 3
	       embed_interpolate (point, p, b), // 4
	       area*Delta, // 5
	       embed_vorticity (point, u, b, n) // 6
	       );
      fflush (fp);
    }
}

event surface_profiles (t = {0.5,0.9,1.5,2.5})
{
  char name3[80];
  sprintf (name3, "cp-t-%g-pid-%d", t, pid());
  FILE * fp = fopen (name3, "w");
  cpout (fp);
  fclose (fp);
}

/**
## Results

~~~gnuplot Time history of drag coefficient. Re = 1000. See Fig. 1a of [Mohaghegh et al. 2017](#mohahegh2017).
reset
set terminal svg font ",16"
set key font ",14" top spacing 1
set grid
set xlabel 't/(d/u)'
set ylabel 'C_D'
set xrange [0:3]
set yrange [0:2.5]
plot '../data/starting/fig1a.SIM' u ($1/2.):2 w l lw 2     lc rgb 'black' t 'Mohaghegh et al., SIM, 2017, C_D', \
     '../data/starting/fig1a.KL'  u ($1/2.):2 pt 7 ps 0.7  lc rgb 'black' t 'K. and L., 1995, C_D', \
     '../data/starting/fig19.p'   u ($1/2.):2 pt 11 ps 0.7 lc rgb 'black' t 'K. and L., 1995, F_p', \
     '../data/starting/fig19.f'   u ($1/2.):2 pt 9 ps 0.7  lc rgb 'black' t 'K. and L., 1995, F_{mu}', \
     'log' u 2:(4.*($14+$16)) w l lw 1 lc rgb 'blue'      t 'Basilisk, l=12', \
     ''    u 2:(4.*$14)       w l lw 1 lc rgb 'red'       notitle, \
     ''    u 2:(4.*$16)       w l lw 1 lc rgb 'sea-green' notitle
~~~

Note that the points of Figure 4 of [K. & L. 1995](#koumoutsakos1995)
do not seem to match the data in Fig. 5a and 5b of the same paper,
which explains the disagreement in the figure below. This agreement
should be much better as can be seen on the more detailed surface
vorticity figures below.

~~~gnuplot Location of the points of zero surface vorticity. Re = 1000. See Fig. 4 of [K. & L. 1995](#koumoutsakos1995).
set key font ",16" bottom right spacing 1.1
set ylabel 'angle/pi'
set yrange [0:0.5]
plot '../data/starting/fig4.1000.KL' u ($1/2.):2 pt 7 ps 0.7 lc rgb "black" t 'K. and L., 1995', \
     'omega.dat' u 1:($2/pi) pt 5 ps 0.25 lc rgb "blue"      t 'Basilisk, l=12'
~~~

~~~gnuplot Surface vorticity at $t/(d/u)=0.5$. Re = 1000. See Fig. 5a of [Mohaghegh et al. 2017](#mohahegh2017).
set key top right
set xlabel 'theta/pi'
set ylabel 'omega_s/(d/u)'
set xrange [0:1]
set yrange [-100:40]
plot '../data/starting/fig5a.SIM' every 20 pt 7 ps 0.7 lc rgb "black" t 'Mohaghegh et al., SIM, 2017', \
     '../data/starting/fig5a.KL' w l lw 2 lc rgb "black" t 'K. and L., 1995', \
     '< cat cp-t-0.5-pid-* | sort -k3,4 | awk -f ../data/starting/surface.awk' w l lw 2 lc rgb "blue" \
     t 'Basilisk, l=12'
~~~

~~~gnuplot Surface vorticity at $t/(d/u)=1.5$. Re = 1000. See Fig. 5b of [Mohaghegh et al. 2017](#mohahegh2017).
set yrange [-100:100]
plot '../data/starting/fig5b.SIM' every 20 pt 7 ps 0.7 lc rgb "black" t 'Mohaghegh et al., SIM, 2017', \
     '../data/starting/fig5b.KL' w l lw 2 lc rgb "black" t 'K. and L., 1995', \
     '< cat cp-t-1.5-pid-* | sort -k3,4 | awk -f ../data/starting/surface.awk' w l lw 2 lc rgb "blue" \
     t 'Basilisk, l=12'
~~~

~~~gnuplot Time evolution of the timestep $\Delta t$
set xlabel "t/(d/u)"
set ylabel "dt"
set xrange [0:3]
set yrange [1.e-5:1.e-2]
set logscale y
plot 'log' u 2:3 w l lw 2 lc rgb 'blue' t 'Basilisk, l=12'
~~~

## References

~~~bib
@article{bouard1980,
  title={The early stage of development of the wake behind an 
         impulsively started cylinder for 40 < {Re} < 10^4^},
  author={Bouard, Roger and Coutanceau, Madeleine},
  journal={Journal of Fluid Mechanics},
  volume={101},
  number={3},
  pages={583--607},
  year={1980},
  publisher={Cambridge University Press}
}

@article{koumoutsakos1995,
  title={High-resolution simulations of the flow around an 
         impulsively started cylinder using vortex methods},
  author={Koumoutsakos, Petros and Leonard, A},
  journal={Journal of Fluid Mechanics},
  volume={296},
  pages={1--38},
  year={1995},
  publisher={Cambridge University Press}
}

@article{mohaghegh2017,
  title={Comparison of sharp and smoothed interface methods for simulation
         of particulate flows II: Inertial and added mass effects},
  author={Mohaghegh, Fazlolah and Udaykumar, HS},
  journal={Computers \& Fluids},
  volume={143},
  pages={103--119},
  year={2017},
  publisher={Elsevier}
}
~~~
*/
