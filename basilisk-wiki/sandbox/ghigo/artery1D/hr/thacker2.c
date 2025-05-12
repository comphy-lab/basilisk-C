/**
# Thacker solution

We here solve, using the 1D blood flow equations, the inviscid
pressure oscillations in a parabolic aneurysm. This test case is
greatly inspired by the solutions of [Thacker and Sampson](#ref),
well-known in the shallow water community.

## Analytic solution

We choose the following spatial variations of the neutral
cross-sectional area $A_0$ and the arterial wall rigidity $K^\prime$
in the domain $x\in \left[-a,\, a \right]$:

$$
\left\{
\begin{aligned}
& R_0 = \bar{R_0}\left[1 + \Delta
\mathcal{G} \left[1 - \frac{x^2}{a^2} \right] \right] && \text{ with }
\bar{R_0} >0 ,\, a >0 ,\, \Delta \mathcal{G} > -1 \\
& K^\prime =
\mathrm{cst} && \text{ with } K^\prime >0,
\end{aligned}
\right.
$$
where $\Delta \mathcal{G} > 0$. Simple analysis allows us to obtain
the following analytic solution:

$$
\left\{
\begin{aligned}
& U_0 = U \sin \left( \frac{t}{\tau}
\right) && \text{ with } U = \mathrm{cst} \\
& p = - \frac{1}{4} \rho
U^2 \cos \left( 2 \frac{t}{\tau} \right) - \rho \frac{x}{\tau} U \cos
\left( \frac{t}{\tau} \right),
\end{aligned}
\right.
$$
where: 
$$
\tau = \lvert \Delta \mathcal{G} \omega^2 \rvert^{-\frac{1}{2}} .
$$
Finally, we choose $U$ with respect to the following nonlinear
stability arguments, namely that the radius $R$ remains positive and
that the flow remains subcritical. */

#include "grid/multigrid1D.h"
#include "../bloodflow-hr.h"

#define BGHOSTS 2

/**
We define the artery's geometrical and mechanical properties. */

#define R0 (1.)
#define K0 (1.e4)

/**
Next, we define the thacker solution. */

#define BTH (4.)
#define DELTATH (0.1)
#define ATH (1.)

#define KTH ((K0)*sqrt (pi))
#define CTH (sqrt (0.5*(KTH)*(R0)))
#define WTH (2.*(CTH)/(BTH))
#define TAUTH (1./sqrt (abs ((DELTATH)*sq (WTH))))
#define TTH (2.*pi*(TAUTH))

#define r0thacker(x) ((R0)*(1. + (DELTATH)*(1. - sq ((x)/(BTH)))))
#define pthacker(t,x) (-0.25*(1.)*sq ((ATH))*cos (2.*(t)/(TAUTH)) - (1.)*(ATH)*(x)/(TAUTH)*cos ((t)/(TAUTH)))
#define rthacker(t,x) (1./(KTH)*(pthacker ((t), (x))) + (r0thacker ((x))))
#define athacker(t,x) (pi*sq (rthacker ((t), (x))))
#define uthacker(t) ((ATH)*sin ((t)/(TAUTH)))
#define qthacker(t,x) (uthacker ((t))*athacker ((t),(x)))

int main() {

  /**
  The domain is 2BTH.*/

  L0 = 2.*(BTH);
  size (L0);
  origin (-(L0)/2.);
  
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

We impose the Thacker solution for all variables. */

k[left] = (K0);
zb[left] = (K0)*sqrt (pi)*(r0thacker (x));
eta[left] = (K0)*sqrt (athacker (t, x)) - (K0)*sqrt (pi)*(r0thacker (x));
a[left] = (athacker(t, x)); // Used in elevation.h with TREE
q[left] = (qthacker (t, x));

k[right] = (K0);
zb[right] = (K0)*sqrt (pi)*(r0thacker (x));
eta[right] = (K0)*sqrt (athacker (t, x)) - (K0)*sqrt (pi)*(r0thacker (x));
a[right] = (athacker(t, x)); // Used in elevation.h with TREE
q[right] = (qthacker (t, x));

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
    zb[] = k[]*sqrt (pi)*(r0thacker (x));
    a[] = (athacker (0., x));
    q[] = (qthacker (0., x));
  }
}

/**
## Post-processing

We output the computed fields. */

event field (t = {0., 0.1*0.422653, 0.3*0.422653, 0.6*0.422653, 0.8*0.422653})
{
  if (N == 128) {

    char name[80];
    sprintf (name, "fields-%.1f-pid-%d.dat", t/0.422653, pid());
    FILE * ff = fopen (name, "w");

    foreach()
      fprintf (ff, "%g %g %g %g %g %g %g\n",
	       x, k[]/(K0), sq (zb[]/k[])/(pi*sq ((R0))),
	       (athacker (t,x) - sq (zb[]/k[]))/(pi*sq ((R0))), (a[] - sq (zb[]/k[]))/(pi*sq ((R0))),
	       (qthacker (t,x)), q[]
	       );
  }
}

/**
Next, we compute the spatial error for the flow rate. */

event error (t = 0.6*0.422653)
{

  scalar err_a[], err_q[];
  foreach() {
    err_a[] = fabs (a[] - (athacker (t, x)));
    err_q[] = fabs (q[] - (qthacker (t, x)));
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

event stop_run (t = 1)
{
  return 0;
}

/**
# Results for second order

#### Arterial properties

~~~gnuplot $a_0$ and $k$ for $N=128$
reset
set key top right
set xlabel 'x'
set ylabel 'a_0,k'
plot '< cat fields-0.0-pid-*' u 1:3 w l lw 2 lc rgb 'blue' t 'a_0/a_0(0)', \
     '< cat fields-0.0-pid-*' u 1:2 w l lw 2 lc rgb 'red' t 'k/k(0)'
~~~

#### Cross-sectional area and flow rate

Next, we plot the spatial evolution of the cross-sectional area $a$
and flow rate $q$ at $t={0, 0.1, 0.2, 0.3, 0.4}*0.422653$ for $N=128$.

~~~gnuplot $a/a_0$ for $N=128$. 
reset
set xlabel 'x'
set ylabel '(a - a_0)/a_0'
plot '< cat fields-0.0-pid-*' u 1:4 w l lw 3 lc rgb "black" t 'analytic', \
     '< cat fields-0.1-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.3-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.6-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.8-pid-*' u 1:4 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:5 w l lw 2 lc rgb "blue" t 't=0', \
     '< cat fields-0.1-pid-*' u 1:5 w l lw 2 lc rgb "red" t 't=0.1', \
     '< cat fields-0.3-pid-*' u 1:5 w l lw 2 lc rgb "sea-green" t 't=0.3', \
     '< cat fields-0.6-pid-*' u 1:5 w l lw 2 lc rgb "coral" t 't=0.6', \
     '< cat fields-0.8-pid-*' u 1:5 w l lw 2 lc rgb "dark-violet" t 't=0.8'
~~~

~~~gnuplot $q$ for $N=128$. 
reset
set key bottom right
set xlabel 'x'
set ylabel 'q'
plot '< cat fields-0.0-pid-*' u 1:6 w l lw 3 lc rgb "black" t 'analytic', \
     '< cat fields-0.1-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.3-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.6-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.8-pid-*' u 1:6 w l lw 3 lc rgb "black" notitle, \
     '< cat fields-0.0-pid-*' u 1:7 w l lw 2 lc rgb "blue" t 't=0', \
     '< cat fields-0.1-pid-*' u 1:7 w l lw 2 lc rgb "red" t 't=0.1', \
     '< cat fields-0.3-pid-*' u 1:7 w l lw 2 lc rgb "sea-green" t 't=0.3', \
     '< cat fields-0.6-pid-*' u 1:7 w l lw 2 lc rgb "coral" t 't=0.6', \
     '< cat fields-0.8-pid-*' u 1:7 w l lw 2 lc rgb "dark-violet" t 't=0.8'
~~~

#### Convergence

Finally, we plot the evolution of the error for the cross-sectional
area $a$ and the flow rate $q$ with the number of cells $N$. We
compare the results with those obtained with a [first-order
scheme](thacker.c).

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
fit f1(x) '../thacker/log' u (log($1)):(log($2)) via a1, b1
fit f2(x) '../thacker/log' u (log($1)):(log($3)) via a2, b2
fit f3(x) '../thacker/log' u (log($1)):(log($4)) via a3, b3

f11(x) = a11 + b11*x
f22(x) = a22 + b22*x
f33(x) = a33 + b33*x
fit f11(x) 'log' u (log($1)):(log($2)) via a11, b11
fit f22(x) 'log' u (log($1)):(log($3)) via a22, b22
fit f33(x) 'log' u (log($1)):(log($4)) via a33, b33

plot '../thacker/log' u 1:2 w p pt 6 ps 1.5 lc rgb "blue" t '|a|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:2 w p pt 6 ps 1.5 lc rgb "sea-green" t '|a|_1, '.ftitle(a11, b11), \
     exp (f11(log(x))) ls 1 lc rgb "coral" notitle, \
     '../thacker/log' u 1:3 w p pt 7 ps 1.5 lc rgb "navy" t '|a|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:3 w p pt 7 ps 1.5 lc rgb "dark-green" t '|a|_2, '.ftitle(a22, b22), \
     exp (f22(log(x))) ls 1 lc rgb "coral" notitle, \
     '../thacker/log' u 1:4 w p pt 5 ps 1.5 lc rgb "skyblue" t '|a|_{max}, '.ftitle(a3, b3), \
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
fit f1(x) '../thacker/log' u (log($1)):(log($5)) via a1, b1
fit f2(x) '../thacker/log' u (log($1)):(log($6)) via a2, b2
fit f3(x) '../thacker/log' u (log($1)):(log($7)) via a3, b3

f11(x) = a11 + b11*x
f22(x) = a22 + b22*x
f33(x) = a33 + b33*x
fit f11(x) 'log' u (log($1)):(log($5)) via a11, b11
fit f22(x) 'log' u (log($1)):(log($6)) via a22, b22
fit f33(x) 'log' u (log($1)):(log($7)) via a33, b33

plot '../thacker/log' u 1:5 w p pt 6 ps 1.5 lc rgb "blue" t '|q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:5 w p pt 6 ps 1.5 lc rgb "sea-green" t '|q|_1, '.ftitle(a11, b11), \
     exp (f11(log(x))) ls 1 lc rgb "coral" notitle, \
     '../thacker/log' u 1:6 w p pt 7 ps 1.5 lc rgb "navy" t '|q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:6 w p pt 7 ps 1.5 lc rgb "dark-green" t '|q|_2, '.ftitle(a22, b22), \
     exp (f22(log(x))) ls 1 lc rgb "coral" notitle, \
     '../thacker/log' u 1:7 w p pt 5 ps 1.5 lc rgb "skyblue" t '|q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 lc rgb "red" notitle, \
     'log' u 1:7 w p pt 5 ps 1.5 lc rgb "forest-green" t '|q|_{max}, '.ftitle(a33, b33), \
     exp (f33(log(x))) ls 1 lc rgb "coral" notitle
~~~
*/
 

