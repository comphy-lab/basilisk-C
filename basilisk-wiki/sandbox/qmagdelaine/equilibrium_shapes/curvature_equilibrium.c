/**
# Finding equilibrium shapes - curvature driven flow

We want to find equilibrium shapes without solving the whole Navier-Stokes
equation. Nevertheless we want to keep the incompressible property of the fluid.
Therefore we indroduce an artificial pressure and use the projection method to
ensure that the velocity field is divergence free.

The method can be simplified. The velocity field is directly given by the
forces (to reduce locally the energy) and a pressure gradient :
$$
\alpha\, \mathbf{u} = \mathbf{f} - \nabla{p}
$$
To ensure the incompressibility, we need to have :
$$
\Delta p = \nabla \cdot \mathbf{f}
$$
We solve this equation using the [poisson solver of Basilisk](/src/poisson.h)
and then compute the velocity field with the previous equation. The VOF tracer
is advected thanks to [advection.h](/src/advection.h).

We take here the example of [meanflow.c](/sandbox/popinet/meanflow.c) : a
sinusoidally perturbated circle that want to go back a circle to minimize 
curvature. In [meanflow.c](/sandbox/popinet/meanflow.c), the volume increases. 
We apply the strategy presented above to keep the volume of the drop constant
during the evolution of its shape.

*/

#include "grid/multigrid.h"
#include "advection.h"
#include "vof.h"
#include "curvature.h"
#include "poisson.h"
#include "view.h"

scalar f[];
scalar * tracers = NULL, * interfaces = {f};
double tend = sq(0.35);

int main()
{
  origin (0., 0.);
  N = 128;
  init_grid(N);
  // this is the diffusive timestep limit (with diffusion coeff unity)
  DT = sq(L0/N)/2.;
  run();
}

double circle (double x, double y)
{
  double theta = atan2(y,x);
  double r = sqrt(sq(x) + sq(y));
  return 0.4*(1. + 0.2*cos(8.*theta)) - r;
}

scalar p[];
double initial_volume;

event init (i = 0)
{
  fraction (f, circle(x, y));
  initial_volume = statsf(f).sum;
  foreach()
    p[] = 0.;
}

// compute the local normal of norm 2

coord normal (Point point, scalar c) {
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

// compute the local normal of norm 2 all over the domain

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}

scalar kappa[], divf[], res[];
face vector force[];

event velocity (i++) {

  /**
  ## Curvature induced force computation */
  curvature (f, kappa);
  vector n[];
  compute_normal (f, n);
  face vector face_normal[];
  foreach_face() {
    double k = 0.;
    if (f[] > 0. && f[] < 1.) {
      k = kappa[];
      face_normal.x[] = n.x[];
      if (f[-1] > 0. && f[-1] < 1.) {
        face_normal.x[] = (face_normal.x[] + n.x[-1])/2.;
	      k = (k + kappa[-1])/2.;
      }
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      k = kappa[-1];
      face_normal.x[] = n.x[-1];
    }
    force.x[] = - k*face_normal.x[];
  }
  boundary((scalar *){force});

  /**
  ## Force divergence and pressure computation */

  foreach() {
    divf[] = 0.;
    foreach_dimension()
      divf[] += (force.x[1] - force.x[])/Delta;
  }
  boundary({p, divf});
  poisson (a = p, b = divf, tolerance = 1e-8, nrelax = 1, & res);
  boundary({p, res});
  
  /** 
  ## Pressure gradient and velocity field computation */
  
  face vector gradp[];
  foreach_face()
    gradp.x[] = (p[] - p[-1])/Delta;
  boundary((scalar *){gradp});
  foreach_face()
    u.x[] = force.x[] - gradp.x[]/* - gradp2.x[]*/;
  boundary((scalar *){u});
}

/**
We save the interface */

event output_interface (t = 0.; t += tend/10.; t < tend) {
  output_facets (f, stdout);
}

/**
We save some variables :

* the current time,
* the difference between the current volume of the drop and the initial one,
* the difference between the extrema of the poisson solver residual, 
* the difference between the extrema of the pressure,
* the mean of the norm 2 of the velocity field,

and generate a movie of the interface shape with the VOF tracer, the pressure, 
the residual of the poisson solver and the x component of the velocity field. */

event viewing (t = 0.; t += tend/100.; t < tend) {

  double volume = statsf(f).sum;
  double delta_error = statsf(res).max - statsf(res).min;
  double delta_p = statsf(p).max - statsf(p).min;
  scalar norm_u[];
  foreach()
    norm_u[] = sqrt(sq(uf.x[]) + sq(uf.y[])); 
  double velocity = normf(norm_u).avg;

  static FILE * fpw = fopen ("volume", "w");
  fprintf(fpw, "%g %.3e %.3e %g %g\n", t, volume-initial_volume, delta_error, delta_p,
          velocity);
  
  view (width = 600, height = 600, fov = 35, quat = {0, 0, -0.707, 0.707});

  clear();
  squares ("f", linear = false, map = jet);
  draw_vof("f", lw = 1.5);
  mirror (n = {1., 0., 0.}, alpha = 0.) {
    squares ("p", linear = false, map = jet);
    draw_vof("f", lw = 1.5);
  }
  mirror (n = {0., 1., 0.}, alpha = 0.) {
    squares ("res", linear = false, map = jet);
    draw_vof("f", lw = 1.5);
    mirror (n = {1., 0., 0.}, alpha = 0.) {
      squares ("uf.x", linear = false, map = jet);
      draw_vof("f", lw = 1.5);
    }
  }
  save ("movie.mp4");
}

/**
# Results

![Animation of the interface shape. Top-left : VOF tracer, top-right : pressure,
bottom-left : residual of the poisson solver, bottom-right : x component of the
velocity field.](curvature_equilibrium/movie.mp4)

~~~gnuplot Gain of volume
set format y "%.1e"
plot 'volume' u 1:2 w l t 'volume' lc rgb "#6495ED"
~~~

~~~gnuplot Evolution of the error of the poisson solver
plot 'volume' u 1:3 w l t 'error' lc rgb "#DC143C"
~~~

~~~gnuplot Evolution of pressure
set format y "%g"
plot 'volume' u 1:4 w l t 'pressure' lc rgb "coral"
~~~

~~~gnuplot Evolution of velocity norm average
plot 'volume' u 1:5 w l t 'velocity' lc rgb "coral"
~~~

~~~gnuplot Evolution of the interface
set size ratio -1
plot 'out' w l t '' lc rgb "#6495ED"
~~~
*/