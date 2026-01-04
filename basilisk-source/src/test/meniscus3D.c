/**
# 3D meniscus

An initial perturbed free-surface is contained within a cube. The
contact angles on the faces of the cube are set by imposing a Neumann
boundary condition on $\eta$ as
$$
\frac{\partial \eta}{\partial n} = \frac{1}{\sqrt{R^2 - x^2 - y^2}},
$$
with $n$ the direction normal to the face and $R = \sqrt{3}$. These
boundary conditions correspond to the static equilibrium solution for
the intersection of a sphere of radius $R$ with the cube. The contact
angle varies between $\text{atan} \left( 1 / \sqrt{2} \right) \approx
35$ degrees on the axis of symmetry and 45 degrees ($\partial_n \eta =
1$) in the corners. A constant surface tension is imposed and the
viscous fluid relaxes toward its equilibrium position.

![Animation of the free-surface relaxation](meniscus3D/movie.mp4)

Time-implicit timestepping is used and greatly decreases the number of
timesteps required to reach convergence.

Both the maximum velocity and standard deviation on the curvature
decrease exponentially.

~~~gnuplot Convergence of the maximum velocity and standard deviation on curvature
set xlabel 'Time'
set logscale y
set format y "10^{%L}"
set xrange [0:5]
plot 'log' u 1:3 w l t 'Maximum velocity', \
     'log' u 1:5 w l t 'Standard deviation of curvature'
~~~

The final value of the mean curvature is close to the theoretical
value $1/R = 1/\sqrt{3}$.

~~~gnuplot Convergence of the mean curvature toward its theoretical value
set ylabel 'Mean curvature - 1/R'
plot 'log' u 1:($4/2. - 1./sqrt(3.)) w l t ''
~~~
*/

#include "grid/multigrid.h"
#include "layered/hydro-tension.h"
#include "layered/implicit.h"
#include "view.h"

int main()
{
  L0 = 2.;
  origin (-L0/2., -L0/2.);

  G = 0.;
  nu = 1.;

  /**
  We use a large timestep and decrease the tolerance to get clear
  convergence of the maximum velocity and standard deviation of
  curvature. */
  
  TOLERANCE = 1e-8;
  CFL_H = 40;
  
  N = 64;
  nl = 1;

  run();
}

/**
## Contact angle boundary conditions 

This is necessary to balance the surface tension acceleration on the
side walls. This needs to match the contact angle boundary conditions
below. */

event half_advection (i++)
{
  foreach_dimension() {
    ha.n[left] = 0.;
    ha.n[right] = 0.;
    hu.n[left] = 0.;
    hu.n[right] = 0.;
    hf.n[left] = 0.;
    hf.n[right] = 0.;
  }
}

event init (i = 0)
{
  foreach_dimension() {
    eta[left]  = neumann (1./sqrt(3. - x*x - y*y));
    eta[right] = neumann (1./sqrt(3. - x*x - y*y));
  }
    
  /**
  ## Initial conditions

  The initial free-surface shape. */

  foreach() {
    double H = 0.5 + 0.1*cos(2.*pi*x)*cos(2.*pi*y);
    zb[] = - 0.5;
    foreach_layer()
      h[] = H/nl;
  }
}

/**
## Outputs

We compute the interface curvature (times $\sigma$, which is
unity). This needs to be done in the `pressure()` event since the
$\sigma$ fields used by the `sigma_kappa()` macro are temporary. */

scalar kappa[];

event pressure (i++)
{
  foreach()
    kappa[] = sigma_kappa (eta, 0);
}
  
event logfile (i++)
{
  stats s = statsf (kappa);
  fprintf (stderr, "%g %g %g %g %g %d\n", t, dt, normf(u.x).max,
	   s.sum/s.volume, s.stddev, mgH.i);
}

event movie (i++; t <= 5)
{
  view (fov = 22, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = 0.02, ty = 0., width = 800, height = 600);
  char s[80];
  sprintf (s, "t = %.2f", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.3, max = 0.3);
  save ("movie.mp4");
}
