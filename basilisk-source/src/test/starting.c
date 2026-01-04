/**
# Starting flow around a cylinder

This is a canonical case of complex boundary layer separation,
inspired by the experiments of [Bouard & Coutanceau,
1980](#bouard1980). Notable early numerical simulations include the
results of [Koumoutsakos and Leonard, 1995](#koumoutsakos1995),
hereafter K & L, which will be used in the comparisons below.

We will solve the Navier--Stokes equations and add the cylinder using
an [embedded boundary](/src/embed.h).

The "double projection" method is necessary to avoid noise in the
pressure field which would pollute the drag force results. */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "view.h"

/**
High-resolution is needed to resolve the boundary layers
properly. [Mohaghegh et al., 2017](#mohaghegh2017) propose to use a
maximum resolution of order $D/10/\sqrt{Re}$, with $Re$ the Reynolds
number and $D$ the cylinder diameter. The cylinder diameter will be
set to unity, and the domain size to 18, so that the corresponding
levels of refinement are approximately 12 and 16 for $Re=1000$ and
$Re=9500$, respectively. */

double D = 1.;
int maxlevel = 12;  // 15/16 for Re = 9500, 12 for Re = 1000
double Re = 1000;   // or 9500
double cmax = 3e-3; // 1e-3 for Re = 9500, 3e-3 for Re = 1000

/**
We need a field for viscosity, so that the embedded boundary metric
can be taken into account. */

face vector muv[];

/**
We set a unit velocity inflow on the left and an outflow on the
right. */

double u0 = 1.;

u.n[left] = dirichlet(u0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

/**
Command line arguments can be used to change the default
parameters. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    Re = atof (argv[2]);
  if (argc > 3)
    cmax = atof (argv[3]);
  
  /**
  The domain is $18\times 18$ and only half the cylinder is modelled. */
  
  size (18 [1]);
  origin (- L0/2.);

  /**
  We set the viscosity field and tune the Poisson solver. */
   
  mu = muv;
  TOLERANCE = 1e-4 [*];
  NITERMIN = 2;
  
  run();
}

/**
The viscosity field needs to take the embedded boundary metric into
account. Given that the diameter and velocity are set to one, the
viscosity is just $1/Re$. */

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*u0*D/Re;
}

/**
## Initial conditions
*/

event init (t = 0)
{

  /**
  We can restart from a dump file. */

  if (!restore ("restart")) {

    /**
    Otherwise, we first create a mesh initially refined only around
    the cylinder. */
    
    refine (level <= maxlevel*(1. - sqrt(fabs(sq(x) + sq(y) - sq(D/2.)))/2.));

    /**
    Then initialize the embedded boundary with cylinder. */
    
    solid (cs, fs, sq(x) + sq(y) - sq(D/2.));
    
    foreach()
      u.x[] = u0*cs[]; // fixme: with 1 this results in sub-optimal adaptation
  }
  else { // restart

    /**
    When we restart, we still need to restore the face fraction field
    *fs*, since it is not dumped. */
    
    solid (cs, fs, sq(x) + sq(y) - sq(D/2.));
    
  }

  /**
  The boundary condition is zero velocity on the embedded boundary. */
    
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);  
}

/**
## Positions of the separation points

We would like to track the positions with time of the separation points
on the surface of the cylinder, as done by [K & L,
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
  foreach (serial)
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      embed_geometry (point, &b, &n);
      x += b.x*Delta, y += b.y*Delta;

      ThetaOmega t;
      t.omega = embed_vorticity (point, u, b, n);
      t.theta = atan2(y, x);
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
pointed to by *fp*.  */

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
## Log files

The log file contains the time, timestep, pressure and viscous forces
components exerted by the fluid on the cylinder.

The "omega" file contains the locations where the surface vorticity
vanishes. */

event logfile (i += 10)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (stderr, "%d %g %g %.6f %.6f %.6f %.6f\n",
	   i, t, dt, Fp.x, Fp.y, Fmu.x, Fmu.y);

  static FILE * fp = fopen ("omega", "w");
  omega_zero (fp);
}

/**
## Images and animations

We display the vorticity field and the corresponding adaptive mesh. */

void display_omega (int width, int height)
{
  view (fov = 1.68, tx = -0.0251023,
	ty = 1e-12, // fixme: this is necessary to re-center the view
	width = width, height = height);
  draw_vof ("cs", filled = -1, fc = {1,1,1});
  draw_vof ("cs", lw = 6);
  squares ("omega", min = -12, max = 12, linear = true, map = blue_white_red);
  mirror ({0,1}) {
    draw_vof ("cs", lw = 6);
    squares (color = "level", min = 6, max = maxlevel);
  }
}

/**
Still images are created at times matching those in various papers. */

event snapshots (t = 1.; t <= 3.; t += 1.)
{
  scalar omega[];
  vorticity (u, omega);

  if (t == 1. || t == 3.) {
    // Fig. 4.
    view (fov = 0.470916, tx = -0.0197934, ty = -0.0258659,
	  width = 524, height = 480);
    draw_vof ("cs", filled = -1, fc = {1,1,1});
    draw_vof ("cs", filled = 0, lw = 1);
    double max = Re == 9500 ? 40 : 12;
    squares ("omega", min = -max, max = +max, linear = 1, map = blue_white_red);
    char name[80];
    sprintf (name, "zoom-%g.png", t);
    save (name);

    draw_vof ("cs", lw = 2);
    squares ("level");
    sprintf (name, "cells-%g.png", t);
    save (name);
  }
  
  if (t == 3.) {
    // Fig. 3.
    display_omega (640, 480);
    save ("omega-3.png");
  }
  
  p.nodump = false;
  char name [80];
  sprintf (name, "dump-%g", t);
  dump (name);
}

/**
An animation is created. The iteration interval is adjusted depending
on the maximum spatial resolution. */

event movie (i += 10*(1 << (maxlevel - 12)))
{
  scalar omega[];
  vorticity (u, omega);
  display_omega (1280, 960);
  save ("omega.mp4");
}

/**
## Surface vorticity profiles

We are also interested in the details on the surface of the cylinder,
in particular surface vorticity. */

void cpout (FILE * fp)
{
  foreach (serial)
    if (cs[] > 0. && cs[] < 1.) {
      coord b, n;
      double area = embed_geometry (point, &b, &n);
      x += b.x*Delta, y += b.y*Delta;

      fprintf (fp, "%g %g %g %g %g %g\n",
	       x, // 1
	       y, // 2
	       atan2(y, x), // 3
	       embed_interpolate (point, p, b), // 4
	       area*Delta, // 5
	       embed_vorticity (point, u, b, n)); // 6
    }
}

event surface_profiles (t = {0.5,0.9,1.5,2.5})
{
  char name[80];
  sprintf (name, "cp-%g", t);
  FILE * fp = fopen (name, "w");
  cpout (fp);
  fclose (fp);
}

#if 0
event snapshot (i += 10) {
  scalar omega[];
  vorticity (u, omega);
  p.nodump = false;
  dump();

  FILE * fp = fopen ("cp", "w");
  cpout (fp);
  fclose (fp);
}
#endif

/**
## Adaptive mesh refinement

The mesh is adapted according to the embedded boundary and velocity
field. */

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,cmax,cmax}, maxlevel, 5);
}

/**
## Results

### Re = 1000

![Animation of the vorticity field and adaptive mesh, Re = 1000.](starting/omega.mp4)(width="640" height="480")

The final state at $tU/D = 3$ can be compared with figure 3 (top row) of
[Mohahegh et al. 2017](#mohahegh2017).

~~~gnuplot Time history of drag coefficient. Re = 1000. See Fig. 1a of [Mohahegh et al. 2017](#mohahegh2017).
set xlabel 'tU/D'
set ylabel 'C_D'
set grid
set pointsize 0.5
plot [][0:2]\
     'fig1a.SIM' u ($1/2.):2 w l lw 2 t 'SIM (Mohahegh et al 2017)', \
     'fig1a.KL' u ($1/2.):2 pt 7 t 'K and L. 1995', \
     'fig19.f' u ($1/2.):2 pt 9 t 'friction, K and L. 1995', \
     'fig19.p' u ($1/2.):2 pt 11 t 'pressure, K and L. 1995', \
     'log' u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (12 levels)', \
     'level-13/log' u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (13 levels)', \
     '' u 2:(4.*$6) w l lw 2 t 'friction (13 levels)', \
     '' u 2:(4.*$4) w l lw 2 t 'pressure (13 levels)'     
~~~

Note that the points of Figure 4 of [K. & L. 1995](#koumoutsakos1995)
do not seem to match the data in Fig. 5a and 5b of the same paper,
which explains the disagreement in the figure below. This agreement
should be much better as can be seen on the more detailed surface
vorticity figures below.

~~~gnuplot Location of the points of zero surface vorticity. Re = 1000. See Fig. 4 of [K. & L. 1995](#koumoutsakos1995).
reset
set xlabel 'tU/D'
set ylabel 'Angle/pi'
set grid
set key bottom right
set pointsize 0.5
plot 'fig4.1000.KL' u ($1/2.):2 pt 7 t 'K and L. 1995', \
     'level-13/omega' u 1:($2/pi) pt 5 t 'Basilisk (13 levels)', \
     'omega' u 1:($2/pi) pt 5 t 'Basilisk (12 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=0.5$. Re = 1000. See Fig. 5a of [Mohahegh et al. 2017](#mohahegh2017).
reset
set xlabel 'theta/pi'
set ylabel 'omega_sD/U'
set grid
plot [0:1]'fig5a.SIM' every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     'fig5a.KL' w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 level-13/cp-0.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (13 levels)', \
     '< sort -k3,4 cp-0.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (12 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=1.5$. Re = 1000. See Fig. 5b of [Mohahegh et al. 2017](#mohahegh2017).
plot [0:1]'fig5b.SIM' every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     'fig5b.KL' w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 level-13/cp-1.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (13 levels)', \
     '< sort -k3,4 cp-1.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (12 levels)'
~~~

### Re = 9500

![Animation of the vorticity field and adaptive mesh, Re = 9500, 16 levels.](starting/level-16/omega.mp4)(width="640" height="480")

The final state at $tU/D = 3$ can be compared with figure 3 (bottom row) of
[Mohahegh et al. 2017](#mohahegh2017).

~~~gnuplot Time history of drag coefficient. Re = 9500. See Fig. 1b of [Mohahegh et al. 2017](#mohahegh2017).
set xlabel 'tU/D'
set ylabel 'C_D'
set grid
plot [][0:2.5] \
     'fig1b.SIM' u ($1/2.):2 w l lw 2 t 'SIM (Mohahegh et al 2017)', \
     'fig1b.KL' u ($1/2.):2 pt 7 t 'K and L. 1995', \
     'level-15/log' u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (15 levels)', \
     'level-16/log' u 2:(4.*($4+$6)) w l lw 2 t 'Basilisk (16 levels)'
~~~

~~~gnuplot Location of the points of zero surface vorticity. Re = 9500. See also Fig. 4 of [K. & L. 1995](#koumoutsakos1995).
reset
set term pngcairo enhanced font ",10"
set output 'loc9500.png'
set xlabel 'tU/D'
set ylabel 'Angle/pi'
set grid
set key bottom right
plot 'level-15/omega' u 1:($2/pi) pt 7 ps 0.25 t 'Basilisk (15 levels)', \
     'level-16/omega' u 1:($2/pi) pt 5 ps 0.25 t 'Basilisk (16 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=0.9$. Re = 9500. See Fig. 6a of [Mohahegh et al. 2017](#mohahegh2017).
reset
set term @SVG
set xlabel 'theta/pi'
set ylabel 'omega_sD/U'
set grid
plot [0:1]'fig6a.SIM' every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     'fig6a.KL' w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 level-15/cp-0.9 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (15 levels)', \
     '< sort -k3,4 level-16/cp-0.9 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (16 levels)'
~~~

~~~gnuplot Surface vorticity at $tU/D=2.5$. Re = 9500. See Fig. 6b of [Mohahegh et al. 2017](#mohahegh2017).
plot [0:1]'fig6b.SIM' every 20 pt 7 t 'SIM (Mohahegh et al 2017)', \
     'fig6b.KL' w l lw 2 t 'K and L 1995', \
     '< sort -k3,4 level-15/cp-2.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (15 levels)', \
     '< sort -k3,4 level-16/cp-2.5 | awk -f surface.awk' w l lw 2 \
     t 'Basilisk (16 levels)'
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
