/**
# Marangoni-induced translation due to a temperature gradient

This reproduces the test case in section 3.4 of [Al Saud et al.,
2018](#alsaud2018) which should be consulted for more details. 

![Final velocity field, interface, surface tension gradient and adaptive
 mesh](marangoni/fields.png){ width="100%" }

We reproduce below Figure 13.a of [Al Saud et al., 2018](#alsaud2018),
but for a different capillary number (also note the different vertical
scale).

~~~gnuplot Ratio of the translation velocity and theoretical velocity for different resolutions.
set xlabel 't*'
set ylabel 'u_{drop}/U_{drop}'
set yrange [0.8:]
set key bottom right
plot 'out' i 2 u 5:6 w l t '{/Symbol D} = 1/32 R', \
     'out' i 1 u 5:6 w l t '{/Symbol D} = 1/16 R', \
     'out' i 0 u 5:6 w l t '{/Symbol D} = 1/8 R'
~~~

The rate of convergence and the level of error is also good compared
to Table 4 of [Al Saud et al., 2018](#alsaud2018).

~~~gnuplot Rate of convergence toward the theoretical terminal velocity
reset
set xlabel 'Number of grid points per radius'
set ylabel 'Relative error'
set logscale
set xtics 8,2,32
set xrange [6:40]
set grid
plot 'log' u 1:(($3 - $2)/$3) t '', 2.3*x**-2 t 'x^{-2}'
~~~
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "integral.h"
#include "view.h"

int LEVEL = 8;

/**
See section 3.4 of [Al Saud et al., 2018](#alsaoud2018). Note that we
use a capillary number *Ca* 10 times larger than in Al Saud et al. to
make the computation faster, but the results are good also for
$\text{Ca} = 0.066$. */

const double R = 1. [1], NablaT = 1., Mu = 1., Rho = 1. [0];
const double Re = 0.066, Ca = 0.66;
const double Gamma_T = Re*sq(Mu)/(Rho*sq(R)*NablaT);
const double Gamma_0 = (Gamma_T*R*NablaT)/Ca;
const double t0 = Mu/(Gamma_T*NablaT);
const double Cdrop = 1., Cbulk = 1.;
double U_drop;

/**
We need a field to store the variable surface tension coefficient. */

scalar sigmaf[];

int main()
{
  size (16*R);
  origin (- L0/2.);
  rho1 = rho2 = Rho;
  mu1 = mu2 = Mu;
  d.sigmaf = sigmaf;
  
  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4 [*];
  
  U_drop = - 2./((2. + 3.*mu2/mu1)*(2. + Cdrop/Cbulk))*Gamma_T*R*NablaT/mu1;

  for (LEVEL = 7; LEVEL <= 9; LEVEL++) {
    N = 1 << LEVEL;
    run();
  }
}

/**
We initialize the signed distance *d* and the surface tension gradient. */

event init (t = 0)
{
  foreach() {
    d[] = sqrt (sq(x) + sq(y)) - R;
    sigmaf[] = Gamma_0 + Gamma_T*NablaT*x;
  }
}

/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stdout, "%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	     mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	     mg.nrelax);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

double u_drop = 0.;

event logfile (i += 10)
{
  double xb = 0., vb = 0., sb = 0.;
  static double xb0 = 0., previous = 0.;
  if (t == 0.)
    previous = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vb += u.x[]*dv;
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    sb0 = sb;
    fprintf (stdout, "\nt dsb xb vb/U_drop ta u_drop/U_drop dt perf.t perf.speed\n");
  }
  u_drop = t > previous ? (xb/sb - xb0)/(t - previous) : 0.;
  fprintf (stdout, "%g %g %g %g %g %g %g %g %g ",
	   t/t0, (sb - sb0)/sb0, xb/sb, vb/sb/U_drop,
	   (t + previous)/2./t0, u_drop/U_drop,
	   dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fputc ('\n', stdout);
  xb0 = xb/sb, previous = t;
  fflush (stdout);
}

event graphics (t = 3.*t0)
{
  double U = - Gamma_T*R*NablaT/Mu;
  fprintf (stderr, "%d %g %g\n", N/16, u_drop/U, U_drop/U);
  if (LEVEL == 8) {
    view (fov = 30, near = 0.01, far = 1000,
	  tx = 0.009, ty = -0.076, tz = -0.291,
	  width = 1239, height = 575);
    draw_vof (c = "f", filled = - 1, fc = {1,1,1});
    draw_vof (c = "f", lw = 2);
    squares (color = "sigmaf", spread = 0.8, linear = true);
    vectors (u = "u", scale = 1);
    cells ();
    save ("fields.png");
  }
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-2, 1e-5, 1e-5}, LEVEL);
}
#endif

/**
## References

~~~bib
@hal{alsaud2018, hal-01706565}
~~~
*/
