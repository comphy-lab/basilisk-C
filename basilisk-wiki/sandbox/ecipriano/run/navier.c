/**
# Navier condition on an embedded boundary

We solve the Laplace equation outside a circle of radius $R$, imposing
the Navier (Robin) condition
$$
s_b + \lambda \partial_n s|_b = 0
$$
on the embedded boundary ($\mathbf{n}$ pointing out of the fluid,
i.e. towards the center of the circle) using the
[navier()](/sandbox/ecipriano/src/gradients.h)
function. The exact solution is
$$
s = A + \ln r, \quad A = \frac{\lambda}{R} - \ln R
$$
which is imposed as Dirichlet condition on the domain boundaries,
except on the left one, where the Navier condition is used too. */

#include "grid/multigrid.h"
#include "embed.h"
#include "gradients.h"
#include "poisson.h"

#define R 0.25
#define LAMBDA 0.1
#define A (LAMBDA/R - log(R))
#define exact(x,y) (A + 0.5*log(sq(x) + sq(y)))

int main() {
  for (int N = 32; N <= 256; N *= 2) {
    origin (-L0/2., -L0/2.);
    init_grid (N);

    /**
    The embedded boundary is a circle of radius $R$. */

    vertex scalar phi[];
    foreach_vertex()
      phi[] = sqrt (sq(x) + sq(y)) - R;
    fractions (phi, cs, fs);
    restriction ({cs, fs});
    cm = cs, fm = fs;

    /**
    We set the boundary conditions. On the left boundary $\mathbf{n} =
    -\mathbf{e}_x$, so $\partial_n s = -x/r^2$. */

    scalar s[], rhs[];
    s[embed] = navier (0., LAMBDA);
    s[left] = navier (exact(x,y) - LAMBDA*x/(sq(x) + sq(y)), LAMBDA);
    s[right] = dirichlet (exact(x,y));
    s[top] = dirichlet (exact(x,y));
    s[bottom] = dirichlet (exact(x,y));
    s.third = true;

    foreach()
      s[] = cs[] > 0. ? 0. : nodata, rhs[] = 0.;
    mgstats mg = poisson (s, rhs, alpha = fs, tolerance = 1e-9);

    /**
    We compute the maximum error in the whole domain and in the cut
    cells. */

    double emax = 0., ecut = 0.;
    foreach (reduction(max:emax) reduction(max:ecut))
      if (cs[] > 0.) {
        double e = fabs (s[] - exact(x,y));
        emax = max (emax, e);
        if (cs[] < 1.)
          ecut = max (ecut, e);
      }
    fprintf (stderr, "%d %d %.3e %.3e\n", N, mg.i, emax, ecut);
  }
}

/**
## Results

~~~gnuplot Convergence of the maximum error
set xlabel "Resolution"
set ylabel "Error"

set xr[2**4:2**9]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "log" u 1:3 via a,b

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", a, b)

plot "log" u 1:3 title "max", \
     "log" u 1:4 title "max (cut cells)", \
     f(x) w l title ftitle(a, b)
~~~
*/
