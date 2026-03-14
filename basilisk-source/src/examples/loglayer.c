/**
# Viscous sublayer and log layer for turbulent flow in a channel

We consider a channel bounded by a no-slip wall on the left and a
"slip wall" (i.e. a symmetry plane) on the right.

The flow is assumed to be fully turbulent and the time-averaged
velocity is uniform along the flow direction and parallel to the
walls.

We follow Boussinesq's hypothesis that a "turbulent viscosity" exists
so that the equations of motion (i.e. "Reynolds-averaged
Navier-Stokes") reduce to
$$
\partial_x(D\partial_x U) + G = 0
$$
with $D$ the sum of the fluid's viscosity $\nu$ and a turbulent
viscosity $\nu_t$ and $G$ the driving acceleration.

This is a one-dimensional diffusion problem which can be solved on
adaptive 1D grid (a bi-tree). */

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

/**
The parameters are the acceleration $G$ and the fluid viscosity
$\nu$. Note that we need a small enough value of the viscosity to get
a separation of scales sufficient to capture both the viscous sublayer
and the log layer. $\kappa$ is the von Karman constant. */

const double G = 1., nu = 1e-4;
const double kappa = 0.41 [0];

/**
The numerical parameters are the minimum and maximum level of
refinement. We need a large range of scales i.e. $2^{16} \approx
O(10^5)$ to capture both the viscous sublayer and the log layer. */

const int minlevel = 5, maxlevel = 16;

/**
The timestep cannot be too large. The tolerance controls the accuracy
of the implicit diffusion solution. The initial mesh has
$2^\text{minlevel}$ grid points. */

int main() {
  DT = 1 [0,1];
  N = 1 << minlevel;
  TOLERANCE = 1e-4;
  run(); 
}

/**
The only unknown is the velocity $U$, with the a Dirichlet zero
boundary condition on the left and the default Neumann zero condition
on the right. */

scalar U[];
U[left] = dirichlet(0);

/**
`Un` is just used to track convergence. */

scalar Un[];

event init (t = 0) {
  foreach()
    Un[] = U[];
}

/**
## Time integration 

In this event we solve the diffusion equation
$$
\partial_x(D\partial_x U) + G = 0
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
    double nu_t = sq(l)*(U[] - U[-1])/Delta;
    D.x[] = nu + nu_t;
  }
  mg = diffusion (U, dt, D, r = a);
}

/**
## Mesh adaptation

We do not adapt at every timestep to minimize adaptation noise. */

event adapt (i += 10)
  adapt_wavelet ({U}, {1e-4}, maxlevel, minlevel);  

/**
## Convergence

If the change in $U$ between timesteps is small enough we stop the
simulation. */

event logfile (i++; i <= 400) {
  double du = change (U, Un);
  if (i > 0 && du < 1e-4)
    return 1; /* stop */
  fprintf (stderr, "%d %g %g %g %d\n", i, t, dt, du, mg.i);
}

/**
## Results

We output the final velocity profile. The characteristic "friction
velocity" is $u_\star = \sqrt{GL_0}$ where $L_0$ is the channel width. A
characteristic length scale of the viscous sublayer is $y_\star = \nu/u_\star$. */

event profile (t = end) {
  printf ("\n");
  const double u_star = sqrt(G*L0), y_star = nu/u_star;
  foreach (serial)
    printf ("%g %g %g %d\n", x/y_star, U[]/u_star, Delta/y_star, level);
}

/**
We can then compare the numerical solution with the analytical
solution (see p. 34 of [Lagrée, 2026](#lagree2026)).

~~~gnuplot Numerical and analytical velocity profiles
set xlabel 'y_+'
set ylabel 'u_+'
kappa = 0.41
set logscale x
u_plus(y) = (1./y - sqrt(4.*(y*kappa)**2 + 1.)/y + 2.*kappa*asinh(2.*y*kappa))/(2.*kappa**2)
u_log(y) = 1./kappa*log(y) + (log(kappa) - 1 + log(4.))/kappa
set key bottom right
plot [:][0:21]'out' u 1:2 pt 6 t 'Numerical', \
     u_plus(x) t 'Full analytical solution' lw 2, \
     x t 'Viscous sublayer', \
     u_log(x) t 'Log layer'
~~~

Although the mesh adaptation is not very visible on the graph above,
it is indeed important as illustrated below.

~~~gnuplot Mesh size
set ylabel '\Delta_+'
# unset logscale x
set logscale y
plot 'out' u 1:3 w l t ''
~~~

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
