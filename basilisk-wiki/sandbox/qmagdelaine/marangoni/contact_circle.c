/**
# Example for the integral formulation of the surface tension with contact angles

We show here an example using an integral formulation of the surface tension
force, as depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/),
adapted to VOF.

For the moment, a lot of limitations remain. Nevertheless, this example already
shows that this new implementation of the surface tension allow to describe
contact angles. 

![Pressure field](contact_circle/pressure.mp4)

![Horizontal velocity field](contact_circle/u_x.mp4)


*/

#define VIDEO 1

/**
## Includes */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

/**
In this example, we use [integral_vof_quasi1D_circle.h](integral_vof_quasi1D_circle.h),
which extrapolates the interfaces as circle arcs. With this implementation,
we recover pecfectly the theoretical shape of the interface and converge
nicely toward a still interface: see the results section at the bottom. It is
not so nice if the interfaces are extrepolated with straight lines,
see [contact_line.c](contact_line.c). However, this implementation using circles
is just a first try: a lot of problem remains, notably it is not possible yet to
have straight interface for the moment. */

#include "integral_vof_quasi1D_circle.h"

#include "view.h"
#define BG 0.8 // light gray for background
#define DG 0. // dark gray

/**
## Paremeters */

#define L 4. // size of the box
#define LEVEL 6
#define T_END 100. //4.
#define DELTA_T (T_END/100.)
#define F_ERR 1e-10

/**
We set the value of the contact angle of the liquid on the walls. */

#define liquid_theta (45.*pi/180.)

/**
We define the geometry of the liquid between the two walls. */

#define Rf (L/(2.*cos(liquid_theta)))
#define Yc (1. + L/pow(2.*cos(liquid_theta), 2)*(pi/2. - liquid_theta + cos(liquid_theta)*sin(liquid_theta)))
#define shape(x, y) (Yc - pow((sq(Rf) - sq(x - L/2.)), 0.5) - y)

/**
We define a function to set and update the values of the surface tension on the
interface. We also define a scalar field to store it. */

void surface_tension (scalar gamma) {
  foreach()
    gamma[] = (interfacial(point, f) ? 1. : nodata);
  boundary({gamma});
}
scalar liquid_sigma[];

/**
As with [contact.h](/src/contact.h), we define vector field to store the height
functions assciated to the VOF tracer f. */

vector h[];
FILE * fp = NULL;


int main() {
  size (L);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);

  NITERMIN = 2;

  rho1 = 1.;
  rho2 = rho1/1.;
  mu1 = 1.;
  mu2 = mu1/1.;

  /**
  We associate the height functions $h$ to the VOF tracer $f$. */

  f.height = h;
  run();
}

/**
## Boundary conditions

We set a free flow boundary conditions on top. */

u.n[top] = neumann(0.);
uf.n[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

/**
## Initialization

In the init event, we initialize the surface tension and set different
attributes to the VOF tracer $f$:

* its potentially variable surface tension,
* its contact angle on the boundaries,
* the boundaries where we want to have a contact angle condition
(for the moment, it is only avaible on the left and on the right). */

scalar fi[];

event init (i = 0) {
  fraction (f, shape(x, y));
  boundary({f});

  surface_tension (liquid_sigma);
  f.sigma = liquid_sigma;
  f.theta_left = liquid_theta;
  f.theta_right = liquid_theta;

  foreach()
    fi[] = f[];
  boundary ({fi});

  fp = fopen ("velocity", "w");
}

event vof (i++) {
  if (i == 0) {
    static FILE * fpf = fopen ("init_facets", "w");
    output_facets(f, fpf);
    fflush(fpf);
    static FILE * fps = fopen ("init_shape", "w");
    coord D = {0., 1.};
    scalar position_y[];
    position (f, position_y, D, add=false);
    foreach()
      if (interfacial(point, f)) {
        fprintf(fps, "%g %g\n", x, position_y[]);
        fflush(fps);
      }
  }
}

/**
## Surface tension update

Since the surface tension is defined only on the interface (in view of having 
a non-uniform surface tension), we have to update it at each time step. */

event acceleration (i++) {
  surface_tension (liquid_sigma);
}

/**
## Movie and outputs */

event logfile (i++; t <= T_END)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  
  double df = change (f, fi);
  if (i > 1 && df < F_ERR)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach()
    un[] = norm(u);
  fprintf (fp, "%g %g %g\n", t, normf(un).max, df);
}

event error (t = end) {
  
  /**
  We recompute the reference solution. */
  
  scalar fref[];
  fraction (fref, shape(x, y));
  boundary({fref});
  
  /**
  And compute the maximum error on the curvature *ekmax*, the norm of
  the velocity *un* and the shape error *ec*. */
  
  double ekmax = 0.;
  scalar un[], ef[], kappa[];
  curvature (f, kappa);
  foreach() {
    un[] = norm(u);
    ef[] = f[] - fref[];
    if (kappa[] != nodata) {
      double ek = fabs (kappa[] - (/*AXI*/ + 1.)/Rf);
      if (ek > ekmax)
	ekmax = ek;
    }
  }
  
  /**
  We output these on standard error (i.e. the *log* file). */

  norm ne = normf (ef);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", 
	         LEVEL, liquid_theta, normf(un).max, 
	         ne.avg, ne.rms, ne.max, ekmax);
}

#if VIDEO
event movie (t = 0.; t += DELTA_T) {
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  int ratio = 2;
  
  /**
  We write a first video with the pressure field. */
  
  view (fov = 19.1/ratio, width = ratio*320, height = 320, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5/ratio + Y0/L0), relative = false);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("p", min = -1., max = 1., linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("pressure.mp4");
  
  /**
  We write a second video with the horizontal velocity field. */
  
  view (fov = 19.1/ratio, width = ratio*320, height = 320, samples = 1,
        bg = {BG, BG, BG}, tx = - (0.5 + X0/L0),
        ty = - (0.5/ratio + Y0/L0), relative = false);
  clear();
  draw_vof("f", edges = true, lw = 2., lc = {DG, DG, DG}, filled = 0);
  squares ("u.x", min = -1e-6, max = 1e-6, linear = false, map = cool_warm);
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("u_x.mp4");
}
#endif

/**
At the end of the video, we save the final shape of the interface, using the
two descriptions:

* the PLIC reconstruction,
* the height functions. */

event final_shape (t = end) {
  static FILE * fps = fopen ("final_shape", "w");
  static FILE * fpf = fopen ("final_facets", "w");
  coord D = {0., 1.};
  scalar position_y[];
  position (f, position_y, D, add=false);
  foreach()
    if (interfacial(point, f)) {
      fprintf(fps, "%g %g\n", x, position_y[]);
      fflush(fps);
    }
  output_facets(f, fpf);
  fflush(fpf);
}

/**
# Results

~~~gnuplot Final shapes using the height functions

darkgray="#666666"
purple="#9E5DB6"
blue="#5082DC"
turquoise="#008C7D"
forest="#149632"
orange="#FF780F"
raspberry="#FA0F50"
set style line 1 pt 7 ps 0.7

set terminal @PNG enhanced size 640,640 font ",8"
set output '_final_shape.png'
set key 
set border
set tics
set xlabel "x"
set ylabel "y"
set size ratio 1

h0 = 1.
L0 = 4.
theta = pi/4.
Yc = h0 + L0/(4.*cos(theta))*((pi/2. - theta)/cos(theta) + sin(theta))
r = L0/(2.*cos(theta))

plot \
  '../contact_tension/final_shape' u 1:2 w l lc rgb blue t 'tension.h' , \
  'final_shape' u 1:2 w l lc rgb raspberry t 'integral', \
  'final_shape' u 1:(Yc - (r**2 - ($1 - L0/2.)**2)**(0.5)) \
    w l lc rgb forest t 'model', \
  'init_shape' u 1:2 w l lc rgb orange t 'init'
~~~

The interface looks good for both implementations of the surface tension (and
contact angle).

~~~gnuplot Error on the final shapes with the height functions
set terminal @PNG enhanced size 640,640 font ",8"
set output '_delta_final_shape.png'
set xlabel "x"
set ylabel "y"
set size ratio 1

plot \
  '../contact_tension/final_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'tension.h - model' , \
  'final_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb raspberry t 'integral - model', \
  'init_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb orange t 'init - model'
~~~

Looking at the error (the difference between the obtained interface and the
theoretical position of the interface), we see that [tension.h](tension.h) and
[contact.h](contact.h) does not give exactly the right shape, whereas the
integral formulation seems perfect.

~~~gnuplot max velo
set terminal @PNG enhanced size 640,640 font ",8"
set output '_velocity.png'
set xlabel "t"
set ylabel "u"
set logscale y
set format y "%.2e"
set size ratio 1
plot '../contact_tension/velocity' u 1:2 w p ls 1 lc rgb blue t 'tension.h', \
     'velocity' u 1:2 w p ls 1 lc rgb raspberry t 'integral'
~~~

~~~gnuplot vof changes
set terminal @PNG enhanced size 640,640 font ",8"
set output '_vof_changes.png'
set xlabel "t"
set ylabel "df"
set size ratio 1
plot '../contact_tension/velocity' u 1:3 w p ls 1 lc rgb blue t 'tension.h', \
     'velocity' u 1:3 w p ls 1 lc rgb raspberry t 'integral'
~~~

Both implementations of the surface tension seem to have nice convergence toward
a perfect still meniscus. It is not the case if the interface is extrapolated
with straight lines, see [contact_line.h](contact_line.h). */
