/**
# Analytic solution proposed by O. Delestre

We solve the inviscid 1D blood flow equations in a straight
artery. The computed blood flow is compared to an
analytic solution. */

/**
# Analytic solution

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

We finally obtain the following analytic solution for the cross-sectional area $A$ and the flow rate $Q$:

$$
\left\{
\begin{aligned}
& A\left(x,\, t\right) = \frac{C_2}{C_1 + t}\\
& U\left(x,\, t\right) = \frac{C_3 + x}{C_1 + t} \\
& Q\left(x,\, t\right) = AU,
\end{aligned}
\right.
$$

where $C_1$, $C_2$ and $C_3$ are constants chosen through the inlet and outlet boundary conditions. */


/**
# Code
*/

#include "grid/cartesian1D.h"
#include "bloodflow.h"

double R0 = 1.;
double K0 = 1.e4;
double L = 10.;

double C1 = -1.;
double C2 = -pi;
double C3 = 5.;

scalar ea[], eq[];
double ea1, ea2, eamax, eq1, eq2, eqmax;
int ne = 0;

double celerity (double a, double k) {

  return sqrt(0.5*k*sqrt(a));
}

double analyticA (double t, double x,
		  double C1, double C2) {

  return (C2)/(C1 + t);
}

double analyticU (double t, double x,
		  double C1, double C3) {

  return (C3 + x)/(C1 + t);
}

int main() {

  origin (0., 0.);
  L0 = L;

  for (N = 32; N <= 256; N *= 2) {

    ea1 = ea2 = eamax = 0.;
    eq1 = eq2 = eqmax = 0.;
    ne = 0;
    
    run();

    printf ("ea, %d, %g, %g, %g\n", N, ea1/ne,
	    sqrt(ea2/ne), eamax);
    printf ("eq, %d, %g, %g, %g\n", N, eq1/ne,
	    sqrt(eq2/ne), eqmax);
  }
  
  return 0; 
}

q[left] = dirichlet(analyticA (t, 0., C1, C2)*analyticU (t, 0., C1, C3));
a[right] = dirichlet(analyticA (t, L0, C1, C2));

event defaults (i=0) {

  gradient = order1;
  riemann = hll_glu;
}

event init (i=0) {
  
  foreach() {
    k[] = K0;
    a0[] = pi*pow(R0, 2.);
    a[] = analyticA (0., x, C1, C2);
    q[] = analyticA (0., x, C1, C2)*analyticU (0., x, C1, C3);
  }
}

event field (t = {0., 0.1, 0.2, 0.3, 0.4}) {

  if (N == 128) {
    foreach() {
      fprintf (stderr, "%g, %.6f, %.6f\n", x, a[], q[]) ;
    }
    fprintf (stderr, "\n\n") ;
  }
}

event error (i++) {

  ne++;
  
  foreach() {
    ea[] = a[] - analyticA (t, x, C1, C2);
    eq[] = q[] - analyticA (t, x, C1, C2)*
      analyticU (t, x, C1, C3);
  }
  
  norm na = normf (ea);
  ea1 += na.avg;
  ea2 += na.rms*na.rms;
  if (na.max > eamax)
    eamax = na.max;
  
  norm nq = normf (eq);
  eq1 += nq.avg;
  eq2 += nq.rms*nq.rms;
  if (nq.max > eqmax)
    eqmax = nq.max;
}

event end (t = 0.4) {
  printf ("#Inviscid Delestre test case completed\n") ;
}

/**
# Plots
*/

/**
~~~gnuplot Spatial evolution of the cross-sectional area $A$ computed at $t={0, 0.1, 0.2, 0.3, 0.4}$ for $N=128$ 

reset

mycolors = "dark-blue red sea-green dark-violet orange"
mypoints = '1 2 3 4 6'

set datafile separator ','
set style line 1 lt -1 lw 2 lc 'black'

k0 = 1.e4
a0 = pi

c1 = -1.
c2 = -pi
c3 = 5.

a(t, x) = c2/(c1 + t)
u(t, x) = (c3 + x)/(c1 +t)
q(t, x) = a(t, x)*u(t, x)

set output 'A.png'
set xlabel 'x [cm]'
set ylabel 'A [cm^2]'
set key top right

plot a(0., x) w l ls 1 t 'Analytic', \
     a(0.1, x) w l ls 1 notitle, \
     a(0.2, x) w l ls 1 notitle, \
     a(0.3, x) w l ls 1 notitle, \
     a(0.4, x) w l ls 1 notitle, \
     'log' i 0 u 1:2 every 5 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     'log' i 1 u 1:2 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.1', \
     'log' i 2 u 1:2 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.2', \
     'log' i 3 u 1:2 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.3', \
     'log' i 4 u 1:2 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.4'
~~~

~~~gnuplot Spatial evolution of the flow rate $Q$ computed at $t={0, 0.1, 0.2, 0.3, 0.4}$ for $N=128$ 

set output 'Q.png'
set xlabel 'x [cm]'
set ylabel 'Q [cm^3/s]'
set key top right

plot q(0., x) w l ls 1 t 'Analytic', \
     q(0.1, x) w l ls 1 notitle, \
     q(0.2, x) w l ls 1 notitle, \
     q(0.3, x) w l ls 1 notitle, \
     q(0.4, x) w l ls 1 notitle, \
     'log' i 0 u 1:3 every 5 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     'log' i 1 u 1:3 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.1', \
     'log' i 2 u 1:3 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.2', \
     'log' i 3 u 1:3 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.3', \
     'log' i 4 u 1:3 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.4'
~~~

~~~gnuplot Comparison of the evolution of the relative error norms for the cross-sectional area $A$ with the number of cells $N$

set xtics 32,2,256
set logscale

ftitle(a,b) = sprintf('Order %4.2f', -b)
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x

fit f1(x) '< grep ^ea out' u (log($2)):(log($3)) via a1, b1
fit f2(x) '< grep ^ea out' u (log($2)):(log($4)) via a2, b2
fit f3(x) '< grep ^ea out' u (log($2)):(log($5)) via a3, b3

set output 'errA.png'
set xlabel 'N'
set ylabel 'Relative error norms'
set key top right

set cbrange[1e-10:1e10]

plot '< grep ^ea out' u 2:3 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|A|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^ea out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|A|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^ea out' u 2:5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|A|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~

~~~gnuplot Comparison of the evolution of the relative error norms for the flow rate $Q$ with the number of cells $N$

fit f1(x) '< grep ^eq out' u (log($2)):(log($3)) via a1, b1
fit f2(x) '< grep ^eq out' u (log($2)):(log($4)) via a2, b2
fit f3(x) '< grep ^eq out' u (log($2)):(log($5)) via a3, b3

set output 'errQ.png'
set xlabel 'N'
set ylabel 'Relative error norms'
set key top right

plot '< grep ^eq out' u 2:3 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|Q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^eq out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|Q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^eq out' u 2:5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|Q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~
*/
