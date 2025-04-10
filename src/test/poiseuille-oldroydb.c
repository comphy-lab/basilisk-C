/**
# Transient planar Poiseuille flow for a viscoelastic Oldroyd-B or FENE-P fluid

The onset of a Poiseuille flow from a fluid at rest after the sudden
application of a pressure gradient is characterized, for a Newtonian fluid, 
by a monotonic exponential increase of the axial velocity.

In the case of a viscoelastic fluid, an elastic oscillation is
superposed on the exponential increase. The problem has an analytical
solution due to [Waters \& King (1970)](#waters1970),

$$
U(Y,T) = 1.5 (1-Y^2) - 48 \sum_{k=1}^\infty \frac{\sin((1+Y)n/2)}{n^3} 
e^{\alpha_n T/2} G(T) 
$$
with $n=(2k-1)\pi$, $\alpha_n = 1 + \beta \,E \, n^2 /4$ and
$$
G(T) =\sinh(\beta_n T/2) + \frac{\gamma_n}{\beta_n} \cosh(\beta_n T/2)
$$
with 
$$
\beta_n = \sqrt{\alpha_n^2 - E \, n^2} \quad \mathrm{and} \quad 
\gamma_n = 1 - \frac{2-\beta}{4} \,E \, n^2
$$

$E$ is the elastic number given by $E = \lambda \mu_o/(\rho h^2)$. 
The time is made dimensionless with $\lambda$, $T = t/\lambda$, 
and the velocity with the average steady velocity,
$$
\bar{u}_\infty = \frac{-\Delta p}{\Delta x} \frac{h^2}{3 \mu_o} 
$$ 

The first modes of the analytical solution imply complex values of
$\beta_n$, which is why we include the complex library. */

#include <complex.h>

#define DT_MAX 0.001
// Total viscoelastic viscosity 
#define MU0 1.          
// Ratio of the solvent viscosity to the total viscosity
#define BETA (1/9.)
// Polymer viscosity
#define MUP ((1. - BETA)*MU0)
// Solvent viscosity 
#define MUS (BETA*MU0)
// Relaxation viscoelastic time
#define LAM 1.
// Average velocity steady flow
#define UAVG (1./(3.*MU0))

/**
Note that we are using as dimensioning scales the width of the gap,
the density of the fluid and the gradient of pressure (rather than the
average velocity for the steady flow). */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#if FENE_P
# include "fene-p.h"
#else
# include "log-conform.h"
# define L2 1.
#endif

int lev;

int main()
{
  size (1.[0]);
  periodic (right);
  p[left] = dirichlet(1.);
  p[right] = dirichlet(0.);

  DT = DT_MAX;
  lambda[] = LAM;
  mup[] = MUP;
  mu[] = {MUS,MUS};
#if FENE_P
  lev = 5;
  init_grid (1 << lev);
  L2 = 10.;   run();
  L2 = 50.;   run();
  L2 = 1000.; run();
#else // Oldroyd-B
  for (lev = 4; lev <= 6; lev++) {
    init_grid (1 << lev);
    run();
  }
#endif
}

u.t[top] = dirichlet(0);
u.t[bottom] = neumann(0);

event init (i = 0) 
{
  scalar s = tau_p.y.y;
  s[top] = dirichlet(0.);
  s = tau_p.x.y;
  s[bottom] = dirichlet(0.);
  foreach()
    u.x[] = 0.;
}

double analytical (double Y, double T, int KF)
{
  double E = LAM*MU0;
  double U = 0.;
  for (int k = 1; k <= KF; k++) {
    double n = (2*k - 1)*M_PI;
    double alpha = 1. + 0.25*BETA*E*sq(n);
    complex double beta = csqrt(sq(alpha) - E*sq(n));
    double gamma = 1. - 0.25*(2. - BETA)*E*sq(n);
    double G = creal(ccosh(0.5*beta*T) + gamma/beta*csinh(0.5*beta*T));
    U += 1./n/sq(n)*sin(0.5*(1. + Y)*n)*exp(-0.5*alpha*T)*G;
  }
  return 1.5*(1. - sq(Y)) - 48.*U;
}

event uaxis_evolution (t += 0.2; t <= 10.) {

  /**
  The  velocity is made dimensionless with the average velocity. 
  when the flow is fully developed, $u_{axis}/ \bar{u}_\infty$ = 3/2. */

  fprintf (stderr, "%g %g %g %d %g\n", t/LAM, 
	   interpolate(u.x, 0.5, 0.)/UAVG, analytical (0, t/LAM, 8), lev, L2);
}

#if FENE_P
event uprofile (t = end) {
  int np = 20;
  double hh = 0.999999999/np;
  for (int j = 0; j <= np; j++)
    printf ("%g %g %d %g\n", hh*j, interpolate(u.x, 0.5, hh*j)/UAVG, lev, L2);
}
#endif // FENE_P

/**
For an Oldroyd-B fluid the solution converges to the theoretical
solution with grid resolution.

~~~gnuplot Temporal evolution of the axial velocity
set xlabel 't'
set ylabel '{u_x(0,t)}'
set key left
plot "< grep '4 1$' log" u 1:2 t 'Basilisk (16x16)', '' u 1:2 t 'Theory' w lp
~~~

~~~gnuplot Temporal evolution of the error for different grids
set xlabel 't'
set ylabel '{/Symbol e}'
set key top right
p "< grep '4 1$' log" u 1:($2 - $3) w l t '16x16', \
  "< grep '5 1$' log" u 1:($2 - $3) w l t '32x32' , \
  "< grep '6 1$' log" u 1:($2 - $3) w l t '64x64'
~~~

For a FENE-P fluid, the solution converges toward the Oldroyd-B
solution when $L^2 \rightarrow \infty$.

~~~gnuplot Temporal evolution of the axial velocity
set ylabel '{u_x(0,t)}'
p "< grep '1000$' ../poiseuille-fenep/log" w l lc 7 t 'Analytical', \
  "< grep '10$' ../poiseuille-fenep/log" w lp ps 0.5 pt 5 t'Fene-P L^2 = 10', \
  "< grep '1000$' ../poiseuille-fenep/log"  w p ps 0.5 pt 5 t 'Fene-P L^2 = 1000'
~~~

~~~gnuplot Velocity profiles at instant t=10
set xlabel '{u_x(y,t=10)}'
set ylabel 'y'
set parametric
set trange [0:1]
set key top right
y(t)=1.5*(1-t*t)
x(t)=t
p   y(t),x(t) lc 7 t 'Newtonian', \
  "< grep '10$' ../poiseuille-fenep/out" u 2:1 w lp pt 5 ps 0.5 t 'L^2 = 10', \
  "< grep '50$' ../poiseuille-fenep/out" u 2:1 w lp pt 5 ps 0.5 t 'L^2 = 50', \
  "< grep '1000$' ../poiseuille-fenep/out" u 2:1 w p pt 5 ps 0.5 t 'L^2 = 1000'
~~~

## References
 
~~~bib
@article{waters1970,
  title={Unsteady flow of an elastico-viscous liquid},
  author={Waters, ND and King, MJ},
  journal={Rheologica Acta},
  volume={9},
  number={3},
  pages={345--355},
  year={1970},
  publisher={Springer}
}

@MastersThesis{gros2013,
  author = {E. Gros},
  title = {Comparing In-House Numerical Simulations to the Analytical
  Solution of the Poiseuille Flow for the {Oldroyd-B} Model},
  school = {Aachen University},
  address = {Germany},
  year = {2013},
  url = {http://www.cats.rwth-aachen.de:8080/theses/ba-gros-2013.pdf}
}
~~~
*/
