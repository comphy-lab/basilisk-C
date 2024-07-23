/**
# Stommel gyre

The theory of [Stommel, 1948](#stommel1948) is one of the first model
of large-scale oceanic circulation. Surface wind stress is balanced by
(linear) bottom friction and gives raise to an asymmetric circulation
(i.e a "gyre"), due to earth's rotation. A narrow, fast "western
boundary current" is balanced by a broad, slow eastern return current
in [geostrophic
balance](https://en.wikipedia.org/wiki/Geostrophic_current), described
by [Sverdrup's theory](#sverdrup1947).

The numerical and analytical solutions for the streamfunction are
displayed below. The upper half is driven by a westerly wind/stress
(i.e. mid-latitude winds) and the bottom half by an easterly
wind/stress (i.e. "trade winds").

~~~gnuplot Numerical (black) and exact (blue) streamfunctions for $N = 64$.
unset key
set size ratio -1
unset surface
set pm3d map
set contour base
set cntrparam levels 10
set cntrlabel onecolor
unset xtics
unset ytics
set zrange [0:2.1e-5]
splot 'out' u 1:2:3 w l lc rgbcolor "black", 'out' u 1:2:4 w l lc rgbcolor "blue"
~~~

The error is mainly localised in the western boundary current and is
controlled by the boundary condition used for the free-surface
gradient. The slope of the free-surface perpendicular to the western
boundary is not zero since the corresponding pressure gradient term
must balance the Coriolis acceleration. Since the discretisation of
this gradient is first-order near the boundary, we expect first-order
convergence, which is indeed what we obtain.

This could be improved by using more sophisticated schemes (e.g.
asymmetric stencils) for computing the gradient near the boundary (see
'2nd-order gradient' on the graph below), however we do not think it
is worth the trouble since this test case over-emphasizes the
importance of boundaries. Real oceans are not bounded by deep vertical
walls, but by shallow seas, usually with very mild bottom slopes, for
which the Coriolis force is balanced by friction terms rather than
pressure against a vertical wall.

This test also checks that a basin limited by a narrow border ('border
= 1') of dry terrain behaves similarly to that limited only by the
numerical boundaries ('border = 0').

~~~gnuplot Error convergence on *u.y*.
reset
ftitle(a,b) = sprintf("%.2g/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) 'log' index 'border = 0' u (log($1)):(log($4)) via a,b
f1(x) = a1 + b1*x
fit f1(x) 'log' index 'border = 1' u (log($1)):(log($4)) via a1,b1
f2(x) = a2 + b2*x
fit f2(x) 'log' index 'border = 0' u (log($1)):(log($3)) via a2,b2
set xlabel 'Resolution'
set logscale
set xtics 32,2,512
set ytics format "% .0e"
set grid ytics
set cbrange [1:2]
set xrange [43:384]
set ylabel 'Error'
set yrange [*:*]
set key above
plot 'log' index 'border = 0' u 1:4 pt 6 t 'max (border = 0)', \
     exp(f(log(x))) t ftitle(a,b),		     \
     'log' index 'border = 0' u 1:3 t 'rms (border = 0)',      \
     exp(f2(log(x))) t ftitle(a2,b2) ,               \
     'log' index 'border = 1' u 1:4 pt 8 t 'max (border = 1)', \
     exp(f1(log(x))) t ftitle(a1,b1),		     \
     'log' index 'border = 1' u 1:3 t 'rms (border = 1)', \
     'stommel-ml.2nd' index 'border = 1' u 1:3 t 'rms (2nd-order gradient)'
~~~

## References

~~~bib
@article{sverdrup1947,
  title={Wind-driven currents in a baroclinic ocean; with application
  to the equatorial currents of the eastern {P}acific},
  author={Sverdrup, Harald Ulrich},
  journal={Proceedings of the National Academy of Sciences of the
  United States of America},
  volume={33},
  number={11},
  pages={318},
  year={1947},
  publisher={National Academy of Sciences}
}

@article{stommel1948,
  title={The westward intensification of wind-driven ocean currents},
  author={Stommel, Henry},
  journal={Eos, Transactions American Geophysical Union},
  volume={29},
  number={2},
  pages={202--206},
  year={1948},
  publisher={Wiley Online Library}
}

@book{pedlosky2013,
  title={Ocean circulation theory},
  author={Pedlosky, Joseph},
  year={2013},
  publisher={Springer Science \& Business Media},
  note={p. 41}
}
~~~

## Numerical setup
*/

#include "grid/multigrid.h"

/**
The only control parameter is the relative width of the western
boundary layer (see [Pedlosky, 2013](#pedlosky2013), page 41)
$$
\delta_S = \frac{r}{\beta}
$$
where $r$ is the linear friction coefficient and $\beta$ the Coriolis
$\beta$-plane coefficient. */

#define DELTA 0.05
#define MAXLEVEL 8
#define MAXEU 1e-7

double tau0 = 1e-5, Beta = 1.;

/**
We use the time-implicit version of the multilayer solver. */

#include "layered/hydro.h"
#include "layered/implicit.h"

/**
We include Coriolis acceleration with a $\beta$-plane variation and a
linear friction coefficient $K_0$. Better results are obtained with a
more implicit discretisation than the default ($\alpha_H = 1/2$). */

#define F0() (Beta*(y - 0.5))
#define K0() (DELTA*Beta)
#define alpha_H 1.
#include "layered/coriolis.h"

double border = 0.;

int main()
{

  /**
  We switch off advection of momentum. */

  linearised = true;
  DT = 0.5;
  TOLERANCE = 1e-9;
  for (border = 0; border <= 1; border++) {
    fprintf (stderr, "\n\n# border = %g\n", border);
    for (N = 64; N <= 256; N *= 2) {
      size (1. + 2.*border/N);
      origin (- border/N, - border/N);
      run();
    }
  }
}

/**
## Wind stress 

We add the zonal wind stress, which only depends on "latitude"
$y$. */

event acceleration (i++)
{
  foreach_face(x)
    if (hf.x[] > dry)
      ha.x[] -= tau0*cos(pi*y);
}

/**
## Theoretical solutions

This is Sverdrup's solution for the streamfunction. */

void sverdrup (scalar psi)
{
  foreach()
    psi[] = tau0*pi*(1. - x)*sin(pi*y);
}

/**
And Stommel's. */

double stommel (double x, double y)
{
  double alpha = sqrt(1./sq(2.*K0()) + sq(pi));
  return tau0/(K0()*pi)*sin(pi*y)*(1. - (exp((1. - x)/(2.*K0()))*sinh(alpha*x) +
				       exp(-x/(2.*K0()))*sinh(alpha*(1. - x)))/
				 sinh(alpha));
}

double stommel_u (double x, double y) {
  return (stommel(x, y - 1e-6) - stommel(x, y + 1e-6))/2e-6;
}

double stommel_v (double x, double y) {
  return (stommel(x + 1e-6, y) - stommel(x - 1e-6, y))/2e-6;
}

/**
## Initial and boundary conditions */

event init (i = 0)
{

  /**
  The "outside" is a dry, high terrain. Note that this is useful only
  when a "border" is not included. */

  foreach_dimension() {
    h[left] = 0.;
    eta[left] = dirichlet(1.);
    h[right] = 0.;
    eta[right] = dirichlet(1.);
  }

  /**
  The "inside" is a basin of constant depth unity. The initial velocity
  field is Stommel's analytical solution. */

  foreach() {
    zb[] = x > 0. && x < 1. && y > 0. && y < 1. ? -1. : 1.;
    h[] = max(0., - zb[]);
    if (h[] > dry) {
      u.x[] = stommel_u (x, y);
      u.y[] = stommel_v (x, y);
    }
  }
}

#if 0 // TREE
event adapt (i++)
{
  adapt_wavelet ((scalar *){u}, (double[]){MAXEU,MAXEU}, MAXLEVEL, 5);
}
#endif
 
/**
## Outputs

We compute the numerical and analytical streamfunctions. */

scalar un[];

void streamfunctions (scalar psi, scalar psim)
{
  foreach_dimension() {
    psi[left] = dirichlet(0);
    psi[right] = dirichlet(0);
  }
  scalar omega[];
  vorticity (u, omega);
  foreach()
    psi[] = 0.;
  poisson (psi, omega);

  foreach()
    psim[] = stommel(x, y)*(x > 0 && x < 1 && y > 0 && y < 1);
}

/**
We stop at convergence on the maximum meridional velocity. */

event logfile (i += 10; t <= 200)
{
  double du = change (u.y, un);

#if 0
  scalar ev[], hw[];
  foreach() {
    ev[] = (u.y[] - stommel_v(x, y))*(x > 0 && x < 1 && y > 0 && y < 1);
    hw[] = h[] > dry ? fabs (h[] - 1.) : nodata;
  }
  norm n = normf (ev);  
  //  fprintf (stderr, "%g %g %g %g %g\n", t, dt, du, n.rms, n.max);

  scalar psi[], psim[];
  streamfunctions (psi, psim);
  dump();
#endif
  
  if (i > 0 && du < 1e-10)
    return 1; /* stop */
}

/**
We compute the error between the analytical and numerical meridional
velocity component. */

event snapshot (t = end)
{
  scalar eu[], ev[];
  foreach() {
    eu[] = (u.x[] + (stommel(x, y + 1e-6) - stommel(x, y - 1e-6))/2e-6)*
      (x > 0 && x < 1 && y > 0 && y < 1);
    ev[] = (u.y[] - (stommel(x + 1e-6, y) - stommel(x - 1e-6, y))/2e-6)*
      (x > 0 && x < 1 && y > 0 && y < 1);
  }
  norm nv = normf (ev), nu = normf (eu);
  fprintf (stderr, "%d %g %.4g %.4g %.4g %.4g\n",
	   N, t, nv.rms, nv.max, nu.rms, nu.max);

  if (N == 64 && border == 0) {
    scalar psi[], psim[];
    streamfunctions (psi, psim);
    dump();
    output_field ({psi, psim});
  }
}
