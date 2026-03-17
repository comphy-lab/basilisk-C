/**
# Turbulent flow in a channel with a "wall model"

This is a variant of this [example](loglayer.c) which solves the
entire turbulent flow profile in a channel with a viscous sublayer and
Prandtl's "mixing length" turbulence model.

The goal is to exclude the small scales and strong gradients
associated with the viscous sublayer and its transition toward the log
profile, to drastically reduce the necessary mesh refinement.

The configuration is the same but the left boundary of the channel is
reduced by a small amount $y_c$ on which a special boundary condition
is applied: a "wall model". */

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

/**
The parameters are the acceleration $G$, the fluid viscosity $\nu$ and
the width of the channel $W$ (including the viscous sublayer etc.). */

const double G = 1., nu = 1e-4, W = 1. [1];
const double kappa = 0.41 [0]; // von Karman constant

/**
The left boundary is set at $y_c$, a small value compared to $W$ but
sufficient to exclude the viscous sublayer, since $y_+ =
y_c\nu/u_\star = 10$ here. */

double yc = 1e-3;

/**
The numerical parameters are the minimum and maximum level of
refinement. */

int minlevel = 5, maxlevel = 11;

/**
## Wall model

To construct the left boundary condition at $y = y_c$, we start from
the full analytical solution (see [Lagrée, 2026](#lagree2026) p.34)
$$
u_+ = \frac{\frac{1}{y_+} - \frac{\sqrt{4y^2_+\kappa^2+1}}{y_+} + 
      2\kappa\sinh^{-1}(2y_+\kappa)}{2\kappa^2}
$$
*/

double u_plus (double y_plus)
{
  return (1./y_plus - sqrt(4.*sq(y_plus*kappa) + 1.)/y_plus +
          2.*kappa*asinh(2.*y_plus*kappa))/(2.*sq(kappa));
}

/**
The gradient of this function is
$$
\partial_{y_+} u_+ = \frac{\sqrt{4\kappa^2y_+^2 + 1} - 1}{2\kappa^2y_+^2}
$$

We want to construct a slip length $\lambda$ relating $U$ and its gradient
$\partial_yU$ at location $y_c$ i.e.
$$
U|_{y_c} = \lambda\partial_yU|_{y_c}
$$
and assuming that $U$ follows the analytical solution above. By
definition we have
$$
\partial_y U|_{y_c} = \frac{u_\star}{y_\star} \partial_{y_+} u_+|_{y_c}
$$
with $u_\star$ the "friction velocity", $y_\star = \nu/u_\star$ and
$y_+ = y/y_\star$. Using the derivative above we can write
$$
\partial_y U|_{y_c} = \frac{u^2_\star}{\nu}  \frac{\sqrt{4 \kappa^2 
  (y_c u_\star / \nu)^2 + 1} - 1}{2 \kappa^2  (y_c u_\star / \nu)^2}
$$
which gives after some manipulations
$$
  u_\star = \frac{\nu}{2 \kappa y_c}  \sqrt{\left( 1 + \frac{2 \kappa^2
  y_c^2}{\nu} \partial_y U|_{y_c} \right)^2 - 1}
$$
*/

double slip (double yc, double dudy, double * u_star)
{
  *u_star = nu/(2.*kappa*yc)*sqrt(sq(1. + 2.*sq(kappa*yc)/nu*fabs(dudy)) - 1.);

  /**
  Knowing $u_\star$ we can now compute $y_{c+} = y_c/y_\star$ and
  $$
  \lambda = \frac{u_\star u_+ (y_{c +})}{\partial_y U|_{y_c}}
  $$
  */
  
  double yc_plus = yc*(*u_star)/nu;
  return fabs(dudy) > 0. ? (*u_star)*u_plus (yc_plus)/fabs(dudy) : 0.;
}

/**
## Setup

The timestep cannot be too large. The tolerance controls the accuracy
of the implicit diffusion solution. The initial mesh has
$2^\text{minlevel}$ grid points. */

int main() {
  DT = 1 [0,1];
  N = 1 << minlevel;

  /**
  The domain size is $W - y_c$ and the left boundary is at $x =
  y_c$. */
  
  size (W - yc);
  origin (yc);
  run();

  /**
  We run a second computation with $y_c = 10^{-2}$ which corresponds to
  $y_{c+} = 100$. We decrease the maximum resolution from 11 to 9. */
  
  yc = 1e-2;
  size (W - yc);
  origin (yc);
  maxlevel = 9;
  N = 1 << minlevel;
  run();
}

/**
## Boundary conditions

The only unknown is the velocity $U$. On the right boundary the
default Neumann zero condition is used. On the left boundary we apply
a Navier/Robin boundary condition
$$
U = \lambda\partial_x U
$$
using the slip length $\lambda$ which will be computed with the
[slip()](#slip) function. */

scalar U[];
double lambda = 0, u_star;
U[left] = navier (0., lambda);

/**
## Time integration 

In this event we solve the diffusion equation
$$
\partial_t U = \partial_x(D\partial_x U) + G
$$
*/

mgstats mg;

event integration (i++)
{

  /**
  The timestep is computed taking account events if necessary and with
  a maximum value bounded by `DT`. */
  
  dt = dtnext (DT);

  /**
  Field `a` is just the acceleration. */
  
  scalar a[];
  foreach()
    a[] = G;

  /**
  We update the friction velocity $u_\star$ and the slip length
  $\lambda$ using only the value of the velocity gradient on the left
  boundary. */
  
  foreach_boundary (left)
    lambda = slip (yc, (U[] - U[-1])/Delta, &u_star);

  /**
  To initialize the diffusion coefficient $D$, we need a model for the
  turbulent viscosity $\nu_t$. We use Prandtl's "mixing length" model
  which assumes that an important characteristic length of the flow is
  the distance to the boundary. The mixing length $l$ is then defined
  as $\kappa x$ and the turbulent viscosity is
  $$
  \nu_t = l^2 \partial_x U
  $$
  The diffusion coefficient is the sum of the fluid and the turbulent
  viscosity. */
  
  face vector D[];
  foreach_face(x) {
    double l = kappa*x;
    double nu_t = sq(l)*fabs (U[] - U[-1])/Delta;
    D.x[] = nu + nu_t;
  }
  mg = diffusion (U, dt, D, r = a);
}

/**
## Mesh adaptation

We do not adapt at every timestep to minimize adaptation noise. */

event adapt (i += 10)
  adapt_wavelet ({U}, {1e-2}, maxlevel, minlevel);  

/**
## Convergence

400 timesteps are necessary/sufficient to reach a stationary
solution. */

event logfile (i++; i <= 400) {
  if (i == 0)
    printf ("\n# yc = %g\n", yc);
  printf ("%g %d %g %g %g\n", t, mg.i, perf.t, lambda, u_star);
}

/**
## Results

We output the final velocity profile. The (analytical) characteristic
"friction velocity" is $u_\star = \sqrt{GW}$. */

event profile (t = end) {
  double u_star = sqrt(G*W), y_star = nu/u_star;
  fprintf (stderr, "\n# yc = %g\n", yc);
  foreach (serial)
    fprintf (stderr, "%g %g %g %d\n", x/y_star, U[]/u_star, Delta/y_star, level);
}

/**
The solutions with $y_c = 10^{-3}$ and $y_c = 10^{-2}$, corresponding
to $y_{c+} = 10$ and $y_{c+} = 100$ agree closely with the analytical
solution.

~~~gnuplot Numerical and analytical velocity profiles
set xlabel 'y_+'
set ylabel 'u_+'
kappa = 0.41
set logscale x
u_plus(y) = (1./y - sqrt(4.*(y*kappa)**2 + 1.)/y + 2.*kappa*asinh(2.*y*kappa))/(2.*kappa**2)
u_log(y) = (log(y) + log(kappa) + log(4.) - 1.)/kappa
set key bottom right
set grid
plot [1:][0:21]\
     'log' index 'yc = 0.001' u 1:2 pt 5 t 'yc = 0.001', \
     'log' index 'yc = 0.01' u 1:2 pt 7 t 'yc = 0.01', \
     u_plus(x) lt 2 t 'Full analytical solution', \
     x lt 3 t 'Viscous sublayer', \
     u_log(x) lt 4 t 'Log layer'
~~~

To highlight the small difference in the center of the channel (which
is not described by the analytical solution) we remove the logscale
and add the numerical result obtained when [also resolving the viscous
sublayer](loglayer.c).

~~~gnuplot Velocity profiles
unset logscale x
plot '../loglayer/log' u 1:2 w l t 'Full numerical solution', \
     'log' index 'yc = 0.001' u 1:2 w l t 'yc = 0.001', \
     'log' index 'yc = 0.01' u 1:2 w l t 'yc = 0.01'
~~~

We can also check that the local friction velocity converges toward
the "global" analytical friction velocity $u_\star = \sqrt{GW} = 1$.

~~~gnuplot
set xlabel 'Time'
set ylabel 'u_*'
plot 'out' index 'yc = 0.001' u 1:5 w l t 'yc = 0.001', \
     'out' index 'yc = 0.01' u 1:5 w l t 'yc = 0.01'
~~~

The distribution of mesh size across the channel illustrates the large
gains obtained when avoiding the resolution of the viscous sublayer.

~~~gnuplot
set xlabel 'y_+'
set ylabel '\Delta_+'
set logscale
set key top left
plot '../loglayer/log' u 1:3 w l t 'Full numerical solution', \
     'log' index 'yc = 0.001' u 1:3 w l t 'yc = 0.001', \
     'log' index 'yc = 0.01' u 1:3 w l lw 2 t 'yc = 0.01'
~~~

The convergence of the multigrid diffusion solver is also much faster
when the large gradients in the viscous sublayer are avoided.

~~~gnuplot
reset
set ylabel '# of multigrid iterations'
set xlabel 'Time'
plot '../loglayer/out' u 1:2 w l t 'Full numerical solution', \
     'out' index 'yc = 0.001' u 1:2 w l t 'yc = 0.001', \
     'out' index 'yc = 0.01' u 1:2 w l lw 2 t 'yc = 0.01'
~~~

All this results in large gains in runtimes.

~~~gnuplot
set ylabel 'Runtime (seconds)'
set xlabel 'Time'
set logscale y
set yrange [1e-4:]
plot '../loglayer/out' u 1:3 w l t 'Full numerical solution', \
     'out' index 'yc = 0.001' u 1:3 w l t 'yc = 0.001', \
     'out' index 'yc = 0.01' u 1:3 w l lw 2 t 'yc = 0.01'
~~~

## See also

* [Viscous sublayer and log layer for turbulent flow in a channel](loglayer.c)

## References

~~~bib
@misc{lagree2026,
  author = {Pierre-Yves Lagrée},
  title = {{Équations de Saint Venant et application aux mouvements de fonds érodables. 
            ``MU4MEF04 - Ondes et Écoulements en milieu naturel", M1 SU}},
  year = 2026,
  url = {http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf}
}
~~~
*/
