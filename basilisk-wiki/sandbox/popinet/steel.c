/**
# How to impose boundary conditions on a masked box
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"

// NOTE: this uses global coordinates for u, i.e. for any condition below
// u.n -> u.x
// u.t -> u.y
// u.r -> u.z

// top is an outflow
u.t[top] = neumann(0);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
u.n[top] = neumann(0);
u.r[top] = neumann(0);

// bottom is a slip wall
u.t[bottom] = dirichlet(0);
u.n[bottom] = neumann(0);
u.r[bottom] = neumann(0);

// right and left are outflows
u.n[left] = neumann(0);
p[left]   = dirichlet(0.);
pf[left]  = dirichlet(0.);
u.t[left] = neumann(0.);
u.r[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
u.t[right] = neumann(0.);
u.r[right] = neumann(0.);

// front and back are slip walls
u.r[front] = dirichlet(0.);
u.t[front] = neumann(0.);
u.n[front] = neumann(0.);

u.r[back] = dirichlet(0.);
u.t[back] = neumann(0.);
u.n[back] = neumann(0.);

// the steel band u.y = 1, u.x = u.z = 0
bid steel;
u.n[steel] = dirichlet(0.);
u.r[steel] = dirichlet(0.);
u.t[steel] = dirichlet(1.);

int main()
{
  init_grid (16);

  origin (-L0/2.,-L0/2.,-L0/2.);
  
  DT = 0.01;
  mu = unityf;
  run ();
}

event init (i = 0) {
  mask (fabs(x) < 0.25 && fabs(z) < 0.25 ? steel : none);

  foreach()
    u.x[] = u.y[] = u.z[] = -1;

  for (scalar s in {u})
    foreach_dimension()
      s.v.x.i = -1; // u is now not a vector

  boundary ((scalar  *){u});

  foreach_boundary (steel)
    printf ("steel %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);

  foreach_boundary (left)
    printf ("left %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
  foreach_boundary (right)
    printf ("right %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
  foreach_boundary (top)
    printf ("top %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
  foreach_boundary (bottom)
    printf ("bottom %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
  foreach_boundary (front)
    printf ("front %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
  foreach_boundary (back)
    printf ("back %g %g %g %g %g %g\n", x, y, z,
	    (u.x[] + u.x[ghost])/2.,
	    (u.y[] + u.y[ghost])/2.,
	    (u.z[] + u.z[ghost])/2.);
}

/**
The values on the boundaries look OK.

~~~gnuplot
set xlabel 'x'
set xlabel 'y'
set xlabel 'z'
scale = 0.1
set view 120, 326, 1, 1
splot '< grep steel out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'steel', '< grep left out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'left', '< grep right out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'right', '< grep top out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'top', '< grep bottom out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'bottom', '< grep front out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'front', '< grep back out' u 2:3:4:($5*scale):($6*scale):($7*scale) w vector t 'back'
~~~
*/
