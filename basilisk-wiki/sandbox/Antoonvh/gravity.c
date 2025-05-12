/**

![Venus and Earth are attracted to each other due to gravity. Image
 courtesy of [universe
 today](http://www.universetoday.com).](https://www.universetoday.com/wp-content/uploads/2009/12/Venus.png)

# Two Liquid Planets in 2D

 On this page we look at the implementation of a body force due to
gravity. The set-up *closely* mimics three aqua planets that float
around in a two dimensional space that interact with each other via
the virtue of gravity.

##Set-up

To model the flow inside the liquid planets and the effect of gravity
we use a continuum discription. Therefore we invoke the Navier-Stokes
solver. To distinguish the planets from their vacuum surroundings, we
use the "two-phase" scheme including surface tension.
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
/**
Because our planets are small in the vast emptiness of space, we need
a resolution corresponding to $8192$ grid cells in each direction. The
number of planets (*np*) is set to three. Furthermore, we declare
fields for a bodyforce (`grav`) that will be calculated from a gravity
potential (`G`). The later requies obvious physically-consistent
boundary conditions.
*/
int maxlevel = 13;
int np = 3;
face vector grav[];
scalar G[];
G[left]   = dirichlet (0);
G[right]  = dirichlet (0);
G[top]    = dirichlet (0);
G[bottom] = dirichlet (0);
/**
The grid is initialized with size $1024\times 1024$ a.u., we let the
reader decide if that is short for astronomical units or arbitrary
units. We also set the fluid properties of the two phases. To omit the
most prominent computational challenges, a small offer is made with
respect to the relevant dimensionless numbers that discribe the
aforementioned water-vacuum system.
*/
int main() {
  N = 256;
  L0 = 1024;
  X0 = Y0 = -L0/2;
  rho1 = 1.;
  rho2 = 0.001;
  mu1 = 0.01;
  mu2 = 0.0;
  f.sigma = 0.05;
  a = grav;
  run(); 
}
/**
### Initialization

To mimic the details of the big bang, we initialize our planets
randomly. Therefore, we set a random seed. Special care is required to
make sure all processors are on the same page when using
non-shared-memory parellization (e.g. MPI).
*/
event init (t = 0) {
  srand (time(NULL));
  double p[4]; // for xp, yp, R, U
  refine (level < 10 && fabs(x) < 20. && fabs(y) < 20.);
  for (int j = 1; j <= np; j++) {
# if _MPI
    if (pid()==0) {
# endif
      p[0] = 20.*noise();
      p[1] = 20.*noise();
      p[2] = 3.*fabs(noise()) + 1.5;
      p[3] = sign(noise())* (0.5 + 20./sq(p[2])*fabs(noise()));
# if _MPI
    }
    MPI_Bcast (&p, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
# endif
    double xp = p[0];
    double yp = p[1];
    double R = p[2];
    double U = p[3];
    refine (level < maxlevel && sq(x - xp)+sq(y - yp) < sq(R + 1) &&
	    sq(x - xp) + sq(y - yp)>sq(R - 1.));
    scalar f1[];
    fraction (f1, sq(R) - sq(x - xp) - sq(y - yp));
    foreach() {
      f[] += f1[];
      u.x[] += U*f1[] *y/(sqrt(sq(x) + sq(y)));
      u.y[] += U*f1[]*-x/(sqrt(sq(x) + sq(y)));;
    }
  }
  foreach()
    f[] = min(f[], 1.);
  boundary (all);
}
/**
### The Gravity Body Force

As mentioned before, a gravitational potential can be evaluated from
the (mass) density field. The relation between the two is discribed by
a Poisson equation that can conviniently be solved for with the
[dedicated solver](http://www.basilisk.fr/src/poisson.h). The body
force is then simply, $F_g \propto -\nabla G$. Remarkably, when using
a domain decomposition for parralel simulations, this method
faciliates the enheritly non-local gravity force to be evaluated
globally, well done Multigrid strategy!
*/
event acceleration (i++) {
  poisson (G, f); 
  boundary ({G});
  foreach_face() 
    grav.x[] -= f[]* (G[] - G[-1])/Delta;
  }
/**
### Sanity test

In order to check if the equations of motions are solved accurately
and if the set-up is physically consistent, we evaluate the centre of
mass of the system and write the result to a file. For a 'realistic'
system this postion should only move along a straight path, since no
external forces are applied. Again, i.e. for two planents in an
infinite and frictionless domain.
*/
event centre_of_mass (t += 1.) {
  double xm = 0;
  double ym = 0;
  double w = 0;
  foreach (reduction(+:xm) reduction(+:ym) reduction(+:w)) {
    xm += x*f[]*sq(Delta);
    ym += y*f[]*sq(Delta);
    w += f[]*sq(Delta);
  }
  static FILE *  fp = fopen ("pos.dat", "w");
  fprintf (fp, "%g\t%g\n", xm/w, ym/w);
}
/**
### Grid Adaptation

Rather unsophisticated, refinement is only based on the estimated
discretization error in the representation of the *planet-fraction*
field.
*/
event adapt (i++)
  adapt_wavelet({f}, (double []){0.01}, maxlevel, 6);
/**
### Movie

We generate a movie.
*/
double min;
event set_min (i = 0; i < 10; i+= 5.)
  min = statsf (G).min;

event bviewer (t += 0.1; t < 50.) {
  view (fov = 3.);
  draw_vof ("f", filled = 1, fc = {0.9, 0.1, 0.9});
  squares ("G", linear = true, min = min , max = min/2.);
  save ("mov.mp4");
}

/**
## Results

The movie below shows the gravity potential and the liquid planets in
magenta;

![Results](gravity/mov.mp4)

Great succes? We check if the centre of mass has moved in a straight
line or not.

~~~gnuplot
 set xlabel 'x[a.u]'
 set ylabel 'y [a.u]'
 set key left top box 1
 set size ratio -1
 plot 'pos.dat' using 1:2 t 'Centre of mass'
~~~

Appearently, the results are not as acurate as we would have liked them
to be.
*/
