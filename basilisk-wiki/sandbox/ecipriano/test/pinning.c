/**
# Pinning a Liquid Droplet

We want to test the module [pinning.h](../src/pinning.h), to suspend
a liquid droplet at a specific point of the domain. The aim is
to compute the equilibrium contact angle and to compare it with
the analytical solution:

$$
  \theta_{2D} = \arccos \left( \dfrac{\rho_l g \pi r_d^2}{2 \sigma}\right)
$$

$$
  \theta_{AXI} = \arccos \left( \dfrac{4\rho_l g r_d^3}{3 \sigma d_f}\right)
$$

given by the balance between gravity and surface tension force.
*/

#if AXI
# include "axi.h"
#endif
#include "navier-stokes/centered.h"
#include "pinning.h"
#define FILTERED 1
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"

/**
## Boundary Conditions

We set the velocity to zero on the `bottom` of the domain, in order to mimic
the effect of a solid fiber. */

u.n[bottom] = dirichlet (0.);
u.t[bottom] = dirichlet (0.);
p[bottom] = neumann (0.);
uf.n[bottom] = 0.;
uf.t[bottom] = 0.;

/**
## Simulation Setup

We initialize a volume fraction `fn` which is used to evaluate the convergence
of the numerical simulation, based on the small variation of volume fraction
between two consecutive time steps. */

#define CONVERGENCE_TOLERANCE 1.e-10
scalar fn[];

/**
We set the diameter of the droplet and the maximum level of refinement. */

int maxlevel = 7;
double D0 = 1.e-3 [*];

int main (void) {

  /**
  We set the density and viscosity to water-air properties, in order to ensure
  density and viscosity ratios which are representative of realistic
  situations. */

  rho1 = 1000., rho2 = 1.;
  mu1 = 1.e-3, mu2 = 1.e-5;

  /**
  We set the length of the domain, and we move the origin
  according to the fiber diameter along y (specific to AXI). */

  L0 = 4.*D0;
  X0 = -0.5*L0;
#if AXI
  Y0 = 0.08*D0;
#endif

  /**
  We set the gravity contribution along the x direction according to Basilisk
  convention for AXI simulations. */

  G.x = -9.81;

  /**
  Data for the `pinning` model: we set the pinning point to the initial
  radius of the droplet, considering the fiber diameter in AXI coordinates. */

  pinning.ap = sqrt (sq (0.5*D0) - sq (Y0));
  pinning.ac = 0.;

  /**
  We run the simulation at different values of surface tension. Shallower
  surface tension values would require the implementation of the normal
  components of the height function as explained in
  [contact.h](/src/contact.h). */

#if AXI
  for (f.sigma = 0.016; f.sigma <= 0.06; f.sigma += 0.004)
#else
  for (f.sigma = 0.006; f.sigma <= 0.06; f.sigma += 0.004)
#endif
  {
    init_grid (1 << maxlevel);
    run();
  }
}

/**
We initialize half liquid droplet on the bottom boundary. */

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

event init (t = 0) {
  fraction (f, circle (x, y, 0.5*D0));

  foreach()
    fn[] = f[];
}

/**
We make sure that the interface is never unrefined. */

event adapt (i++)
{
  scalar fa[];
  foreach() {
    if (interfacial (point, f))
      fa[] = rand();
    else
      fa[] = 0.;
  }
  adapt_wavelet ({fa,u.x,u.y}, (double[]){1.e-3,1.e-3,1.e-3}, maxlevel);
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Maximum Velocity

We output the maximum velocity to verify that the velocity relaxes toward a
null value. */

event logger (i += 100)
{
  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (stdout, "%g %g\n", t, normf(un).max);
}

/**
### Contact Angle and Droplet Shape

We write a function that calculates and write the numerical contact angle from
the height--functions, and we write the droplet shape using the facets. */

void write_theta (void) {
  double ca = 0.;
  foreach_point (pinning.ap, Y0)
    ca = atan (1./(h.x[0,-1] - h.x[]));

  fprintf (stdout, "\n\n");
  fflush (stdout);

  fprintf (stderr, "%g %g\n", f.sigma, fabs (ca*180./pi));
  fflush (stderr);
}

void write_facets (void) {
  char name[80];
  sprintf (name, "facets-%.3f", f.sigma);

  FILE * fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
}

/**
We stop the simulation if the variation of volume fraction between two
consecutive time steps is smaller than 1.e-10, similarly to what is done in
[spurious.c](/src/test/spurious.c). */

event stop (i++)
{
  double dc = change (f, fn);

  if (i > 1 && dc < CONVERGENCE_TOLERANCE) {
    write_theta();
    write_facets();
    return 1;
  }
}

/**
At equilibrium (t = 10 seems sufficient), we compute the
numerical contact angle. */

event end (t = 2.)
{
  write_theta();
  write_facets();
}

/**
## Results

We compare the steady state shape of the droplet at different surface tension
values.

~~~gnuplot Droplet shape at different surface tension values 2D case
reset
set term push
set term @SVG size 640,180
set size ratio -1
set xrange[-1.1e-3:1.1e-3]
set yrange[0:0.5e-3]
unset xtics
unset ytics
unset border
plot 0 lt -1 notitle, \
     "facets-0.006" w l lw 2 t "sigma = 0.006", \
     "facets-0.010" w l lw 2 t "sigma = 0.010", \
     "facets-0.014" w l lw 2 t "sigma = 0.014", \
     "facets-0.030" w l lw 2 t "sigma = 0.030", \
     "facets-0.058" w l lw 2 t "sigma = 0.058"
set term pop
~~~

~~~gnuplot Droplet shape at different surface tension values AXI case
reset
set term push
set term @SVG size 640,180
set size ratio -1
set xrange[-1.1e-3:1.1e-3]
set yrange[0:0.5e-3]
unset xtics
unset ytics
unset border
plot 0.08*1.e-3 lt -1 notitle, \
     "../pinning-axi/facets-0.016" w l lw 2 t "sigma = 0.016", \
     "../pinning-axi/facets-0.020" w l lw 2 t "sigma = 0.020", \
     "../pinning-axi/facets-0.024" w l lw 2 t "sigma = 0.024", \
     "../pinning-axi/facets-0.040" w l lw 2 t "sigma = 0.040", \
     "../pinning-axi/facets-0.056" w l lw 2 t "sigma = 0.056"
set term pop
~~~

We plot the comparison between the equilibrium contact angle
and the one obtained from the numerical simulation.

~~~gnuplot Steady contact angle 2D case
reset
set yrange[20:100]
set xrange[0:0.06]
set xlabel "Surface tension (sigma)"
set ylabel "Contact Angle [deg]"
set grid

f(x) = acos(pi*(0.5e-3)**2*1000.*9.81/2/x)*180/pi # 2D

plot f(x) w l lw 2 t "Theoretical", \
     "log" u 1:2 w p pt 8 ps 1.5 t "Numerical"
~~~

~~~gnuplot Steady contact angle AXI case
reset
set yrange[0:100]
set xrange[0:0.06]
set xlabel "Surface tension (sigma)"
set ylabel "Contact Angle [deg]"
set grid

df = 2.*0.08*1e-3
f(x) = acos(4.*1000*9.81*(0.5e-3**3)/(3*x*df))*180/pi # 3D

plot f(x) w l lw 2 t "Theoretical", \
     "../pinning-axi/log" u 1:2 w p pt 8 ps 1.5 t "Numerical"
~~~

We plot the trend of the maximum velocity which relaxes toward a null value
during the simulation time.

~~~gnuplot Relaxation of the maximum velocity 2D case
reset
set xrange[0:2]
set xlabel "time [s]"
set ylabel "Maximum velocity norm [m/s]"
set grid
set logscale y

plot "out" i 0  w l t "sigma = 0.006", \
     "out" i 4  w l t "sigma = 0.022", \
     "out" i 6  w l t "sigma = 0.03" , \
     "out" i 13 w l t "sigma = 0.058"
~~~

~~~gnuplot Relaxation of the maximum velocity AXI case
set xrange[0:2]

plot "../pinning-axi/out" i 0  w l t "sigma = 0.016", \
     "../pinning-axi/out" i 4  w l t "sigma = 0.032", \
     "../pinning-axi/out" i 7  w l t "sigma = 0.044" , \
     "../pinning-axi/out" i 10 w l t "sigma = 0.056"
~~~
*/

