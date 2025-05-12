/**
# Analytic solution proposed by O. Delestre

We solve the inviscid 1D blood flow equations in a straight
artery. The computed blood flow is compared to an
analytic solution. */

/**
## Analytic solution

We search for a simple solution of the form:
$$
\left\{
\begin{aligned}
& A\left(x,\, t \right) = a\left(t\right)x + b\left(t\right) \\
& U\left(x,\, t \right) = c\left(t\right)x + d\left(t\right).
\end{aligned}
\right.
$$
Injecting these expressions in the inviscid 1D blood flow equations
written in non-conservative form, we obtain the following ordinary
differential equations for the coefficients $a$, $b$, $c$ and $d$:
$$
\left\{
\begin{aligned}
& a = 0 \\
& d^{'} + bc = 0 \\
& c^{'} + c^2 = 0 \\
& d^{'} + cd = 0.
\end{aligned}
\right.
$$
We finally obtain the following analytic solution for the
cross-sectional area $A$ and the flow rate $Q$:
$$
\left\{
\begin{aligned}
& A\left(x,\, t\right) = \frac{C_2}{C_1 + t}\\
& U\left(x,\, t\right) = \frac{C_3 + x}{C_1 + t} \\
& Q\left(x,\, t\right) = AU,
\end{aligned}
\right.
$$
where $C_1$, $C_2$ and $C_3$ are constants chosen through the inlet
and outlet boundary conditions. */

#include "grid/multigrid1D.h"
#include "../bloodflow-hr.h"

#define BGHOSTS 2

/**
We define the artery's geometrical and mechanical properties. */

#define R0 (1.)
#define K0 (1.e4)

int main()
{
  /**
  The domain is 10.*/

  L0 = 10.;
  size (L0);
  origin (0.);

  DT = 1.e-5;

  /**
  We run the computation for different grid sizes. */
  
  for (N = 32; N <= 256; N *= 2) {
    init_grid (N);
    run();
  }
}

/**
## Boundary conditions

We impose the flow rate at the the cross-sectional area. Otherwise, we
use homogeneous Neumann boundary conditions. */

#define C1 (-1.)
#define C2 (-pi)
#define C3 (5.)

#define analytic_a(t) ((C2)/((C1) + (t)))
#define analytic_u(t,x) (((C3) + (x))/((C1) + (t)))
#define analytic_q(t,x) ((analytic_a ((t)))*(analytic_u ((t),(x))))

eta[left] = (K0)*sqrt (analytic_a (t)) - (K0)*sqrt (pi)*(R0);
q[left] = analytic_q (t, x);

eta[right] = (K0)*sqrt (analytic_a (t)) - (K0)*sqrt (pi)*(R0);
q[right] = analytic_q (t, x);

/**
## Defaults conditions
*/

event defaults (i = 0)
{
  gradient = minmod;
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We initialize the variables *k*, *zb*, *a* and *q*. */
  
  foreach() {
    k[] = (K0);
    zb[] = k[]*sqrt (pi)*(R0);
    a[] = (analytic_a (0.));
    q[] = (analytic_q (0., x));
  }
}

/**
## Post-processing

We output the computed fields. */

event field (t = {0., 0.1, 0.2, 0.3, 0.4})
{
  if (N == 128) {
    char name[80];
    sprintf (name, "fields-%.1f-pid-%d.dat", t, pid());
    FILE * ff = fopen (name, "w");
    
    foreach()
      fprintf (ff, "%g %g %g %g %g %g %g\n",
	       x, k[], sq (zb[]/k[]),
	       (analytic_a (t))/(pi*sq ((R0))), a[]/(pi*sq ((R0))),
	       (analytic_q (t,x)), q[]
	       );
  }
}

/**
Next, we compute the spatial error for the cross-sectional area $a$
and flow rate $q$. */

event error (t = 0.4)
{
  scalar err_a[], err_q[];
  foreach() {
    err_a[] = fabs (a[] - (analytic_a (t)));
    err_q[] = fabs (q[] - (analytic_q (t, x)));
  }
  boundary ((scalar *) {err_a, err_q});

  norm na = normf (err_a);
  norm nq = normf (err_q);
  
  fprintf (ferr, "%d %g %g %g %g %g %g\n",
	   N,
	   na.avg, na.rms, na.max,
	   nq.avg, nq.rms, nq.max);
}

/**
## End of simulation */

event end (t = 0.5)
{
  return 0;
}

/**
# Results for second order

#### Cross-sectional area and flow rate

We first plot the spatial evolution of the cross-sectional area $a$ at
$t={0, 0.1, 0.2, 0.3, 0.4}$ for $N=128$.

~~~gnuplot $a/a_0$ for $N=128$. 
reset
set xlabel 'x'
set ylabel 'a/a_0'
plot '< cat fields-0.0-pid-*' u 1:4 w l lw 3 lc rgb "black" t 'analytic', \
     '< cat fields-0.1-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.2-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.3-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.4-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:5 w l lw 2 lc rgb "blue" t 't=0', \
     '< cat fields-0.1-pid-*' u 1:5 w l lw 2 lc rgb "red" t 't=0.1', \
     '< cat fields-0.2-pid-*' u 1:5 w l lw 2 lc rgb "sea-green" t 't=0.2', \
     '< cat fields-0.3-pid-*' u 1:5 w l lw 2 lc rgb "coral" t 't=0.3', \
     '< cat fields-0.4-pid-*' u 1:5 w l lw 2 lc rgb "dark-violet" t 't=0.4'
~~~

~~~gnuplot $q$ for $N=128$. 
reset
set key bottom right
set xlabel 'x'
set ylabel 'q'
plot '< cat fields-0.0-pid-*' u 1:6 w l lw 3 lc rgb "black" t 'analytic', \
     '< cat fields-0.1-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.2-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.3-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.4-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:7 w l lw 2 lc rgb "blue" t 't=0', \
     '< cat fields-0.1-pid-*' u 1:7 w l lw 2 lc rgb "red" t 't=0.1', \
     '< cat fields-0.2-pid-*' u 1:7 w l lw 2 lc rgb "sea-green" t 't=0.2', \
     '< cat fields-0.3-pid-*' u 1:7 w l lw 2 lc rgb "coral" t 't=0.3', \
     '< cat fields-0.4-pid-*' u 1:7 w l lw 2 lc rgb "dark-violet" t 't=0.4'
~~~

#### Convergence

Finally, we plot the evolution of the error for the cross-sectional
area $a$ and the flow rate $q$ with the number of cells $N$.  We
compare the results with those obtained with a [first-order
scheme](delestre.c).

~~~gnuplot Spatial convergence for $a$
reset
set xlabel 'N'
set ylabel 'L_1(a),L_2(a),L_{max}(a)'
set format y '%.1e'
set logscale

ftitle(a,b) = sprintf('order %4.2f', -b)

f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x
fit f1(x) '../delestre/log' u (log($1)):(log($2)) via a1, b1
fit f2(x) '../delestre/log' u (log($1)):(log($3)) via a2, b2
fit f3(x) '../delestre/log' u (log($1)):(log($4)) via a3, b3

f11(x) = a11 + b11*x
f22(x) = a22 + b22*x
f33(x) = a33 + b33*x
fit f11(x) 'log' u (log($1)):(log($2)) via a11, b11
fit f22(x) 'log' u (log($1)):(log($3)) via a22, b22
fit f33(x) 'log' u (log($1)):(log($4)) via a33, b33

plot '../delestre/log' u 1:2 w p pt 6 ps 1.5 lc rgb "blue" t '|a|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:2 w p pt 6 ps 1.5 lc rgb "sea-green" t '|a|_1, '.ftitle(a11, b11), \
     exp (f11(log(x))) ls 1 lc rgb "coral" notitle, \
     '../delestre/log' u 1:3 w p pt 7 ps 1.5 lc rgb "navy" t '|a|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:3 w p pt 7 ps 1.5 lc rgb "dark-green" t '|a|_2, '.ftitle(a22, b22), \
     exp (f22(log(x))) ls 1 lc rgb "coral" notitle, \
     '../delestre/log' u 1:4 w p pt 5 ps 1.5 lc rgb "skyblue" t '|a|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:4 w p pt 5 ps 1.5 lc rgb "forest-green" t '|a|_{max}, '.ftitle(a33, b33), \
     exp (f33(log(x))) ls 1 lc rgb "coral" notitle
~~~

~~~gnuplot Spatial convergence for $q$
reset
set xlabel 'N'
set ylabel 'L_1(q),L_2(q),L_{max}(q)'
set format y '%.1e'
set logscale

ftitle(a,b) = sprintf('order %4.2f', -b)

f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x
fit f1(x) '../delestre/log' u (log($1)):(log($5)) via a1, b1
fit f2(x) '../delestre/log' u (log($1)):(log($6)) via a2, b2
fit f3(x) '../delestre/log' u (log($1)):(log($7)) via a3, b3

f11(x) = a11 + b11*x
f22(x) = a22 + b22*x
f33(x) = a33 + b33*x
fit f11(x) 'log' u (log($1)):(log($5)) via a11, b11
fit f22(x) 'log' u (log($1)):(log($6)) via a22, b22
fit f33(x) 'log' u (log($1)):(log($7)) via a33, b33

plot '../delestre/log' u 1:5 w p pt 6 ps 1.5 lc rgb "blue" t '|q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:5 w p pt 6 ps 1.5 lc rgb "sea-green" t '|q|_1, '.ftitle(a11, b11), \
     exp (f11(log(x))) ls 1 lc rgb "coral" notitle, \
     '../delestre/log' u 1:6 w p pt 7 ps 1.5 lc rgb "navy" t '|q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:6 w p pt 7 ps 1.5 lc rgb "dark-green" t '|q|_2, '.ftitle(a22, b22), \
     exp (f22(log(x))) ls 1 lc rgb "coral" notitle, \
     '../delestre/log' u 1:7 w p pt 5 ps 1.5 lc rgb "skyblue" t '|q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:7 w p pt 5 ps 1.5 lc rgb "forest-green" t '|q|_{max}, '.ftitle(a33, b33), \
     exp (f33(log(x))) ls 1 lc rgb "coral" notitle
~~~
*/
