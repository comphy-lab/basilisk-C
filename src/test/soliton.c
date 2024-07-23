/**
# Green-Naghdi soliton

The [Green-Naghdi](/src/green-naghdi.h) system of equations admits
solitary wave solutions of the form
$$
\eta(\psi) = ah_0\mathrm{sech}^2\left(\psi\frac{\sqrt{3ah_0}}
 {2h_0\sqrt{h_0(1+a)}}\right)
$$
$$
u(\psi) = \frac{c\eta}{h_0+\eta}
$$
with $\psi = x - ct$ and the soliton velocity $c^2=g(1+a)h_0$. */

#include "grid/multigrid1D.h"
#include "green-naghdi.h"

/**
The domain is 700 metres long, the acceleration of gravity is 10
m/s^2^. We need to set the dispersion parameter $\alpha_d$ to one.
We compute the solution in one dimension for a number of grid
points varying between 128 and 1024. */

int main()
{
  X0 = -200.;
  L0 = 700.;
  G = 10.;
  alpha_d = 1.;
  for (N = 128; N <= 1024; N *= 2)
    run();
}

/**
We follow [Le Métayer et al, 2010](/src/references.bib#lemetayer2010)
(section 6.1.2) for the values of $h_0$ and $a$. */

double h0 = 10, a = 0.21;

double sech2 (double x) {
  double a = 2./(exp(x) + exp(-x));
  return a*a;
}

double soliton (double x, double t)
{
  double c = sqrt(G*(1. + a)*h0), psi = x - c*t;
  double k = sqrt(3.*a*h0)/(2.*h0*sqrt(h0*(1. + a)));
  return a*h0*sech2 (k*psi);
}

event init (i = 0)
{
  double c = sqrt(G*(1. + a)*h0);
  foreach() {
    double eta = soliton (x, t);
    h[] = h0 + eta;
    u.x[] = c*eta/(h0 + eta);
  }
}

/**
We output the profiles and reference solution at regular intervals. */

event output (t = {0,7.3,14.6,21.9,29.2}) {
  if (N == 256) {
    foreach()
      fprintf (stdout, "%g %g %g %g\n", x, h[] - h0, u.x[], soliton (x, t));
    fprintf (stdout, "\n");
  }
}

/**
We compute the error between the theoretical and numerical solutions
at $t = 29.2$. */

event error (t = end) {
  scalar e[];
  foreach()
    e[] = h[] - h0 - soliton (x, t);
  fprintf (stderr, "%d %g\n", N, normf(e).max/(a*h0));
}

/**
~~~gnuplot Depth profiles for N = 256.
set grid
set xlabel 'x'
set ylabel 'z'
plot 'out' u 1:4 w l t 'exact', 'out' w l t 'numerical'
~~~

The method has a second-order rate of convergence as expected.

~~~gnuplot Relative error as a function of resolution.
set logscale
set xlabel 'N'
set ylabel 'max|e|/a'
set xtics 128,2,1024
set cbrange [1:2]
set grid
fit [5:] a*x+b 'log' u (log($1)):(log($2)) via a,b
plot [100:1250]'log' u 1:2 pt 7 t '', \
     exp(b)*x**a t sprintf("%.0f/N^{%4.2f}", exp(b), -a)
~~~
*/
