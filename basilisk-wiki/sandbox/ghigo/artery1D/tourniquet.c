/**
# Inviscid tourniquet

We solve the inviscid 1D blood flow equations in a straight
artery. The artery is initially inflated on its left half, mimicking
the action of a tourniquet. At $t>0$, the vessel relaxes towards its
steady-state at rest. The solution of the flow of blood generated by
this elastic relaxation is obtained numerically and compared to the
solution of the corresponding Riemann problem. */

/**
# Analytic solution of the corresponding Riemann problem

The initial conditions are:

$$
A(t=0,\, x) = 
\left\{
\begin{aligned}
& A_L \quad \mathrm{if} \: x < x_m \\
& A_R \quad \mathrm{if} \: x >= x_m
\end{aligned}
\right.
, \qquad Q(t=0,\, x) = 0,
$$

with $A_L > A_R$ and $x_m = 0$. Using the characteristic method, we obtain the following analytic solution:

$$
\begin{aligned}
& \mathrm{if} \: x < x_A = - c_L t \quad
\left\{
\begin{aligned}
& A(t,\, x) = A_L \\
& U(t,\, x) = 0
\end{aligned}
\right.
\\
& \mathrm{if} \: x_A <= x < x_B = (u_M - c_M)t \quad
\left\{
\begin{aligned}
& U(t,\, x) = \frac{4}{5}\frac{x}{t} + \frac{4}{5}c_L \\
& c(t,\, x) = -\frac{1}{5}\frac{x}{t} + \frac{4}{5}c_L
\end{aligned}
\right. \\
& \mathrm{if} \: x_B <= x < x_C = \frac{A_M U_M}{A_M - A_R} t \quad
\left\{
\begin{aligned}
& A(t,\, x) = A_M  \\
& U(t,\, x) = U_M
\end{aligned}
\right.
\\
& \mathrm{if} \: x_C <= x \quad
\left\{
\begin{aligned}
& A(t,\, x) = A_R  \\
& U(t,\, x) = 0
\end{aligned}
\right.
.
\end{aligned}
$$

The variables $A_M$ and $U_M$ are obtain by solving the following system iteratively:

$$
\left\{
\begin{aligned}
& U_M + 4 c_M = 4 c_L\\
& Q_M = \frac{A_M U_M}{A_M - A_R} \left[A_M - A_R\right]\\
& \left[ \frac{Q_M^2}{A_M} + \frac{K}{3\rho} A_M^{\frac{3}{2}} \right] - \frac{K}{3 \rho} A_R^{\frac{3}{2}} = \frac{A_M U_M}{A_M - A_R} Q_M .
\end{aligned}
\right.
$$
*/

/**
# Code
*/

#include "grid/cartesian1D.h"
#include "bloodflow.h"

double R0 = 1.;
double K0 = 1.e4;
double L = 10.;

double XM = 0.;
double DR = 1.e-1;

scalar ea[], eq[];
double ea1, ea2, eamax, eq1, eq2, eqmax;
int ne = 0;

double celerity (double a, double k) {

  return sqrt(0.5*k*sqrt(a));
}

double analyticA (double t, double x, double dr,
		  double a0, double k0) {

  double al = a0*pow(1 + dr, 2.);
  double cl = celerity(al, k0);

  double ar = a0;

  double am = 3.459578046858399;
  double cm = celerity(am, k0);
  double um = 9.192473939896399;
  double s = am*um/(am-ar);

  double xa = -cl*t;
  double xb = (um - cm)*t;
  double xc = s*t;

  if (x <= xa) return al;
  else if (xa < x && x <= xb) {
    return pow(2./k0, 2.)*pow( -1./5.*x/t + 4./5.*cl, 4.);
  }
  else if (xb < x && x <= xc) return am;
  else return ar;	     
}

double analyticU (double t, double x, double dr,
		  double a0, double k0) {

  double al = a0*pow(1 + dr, 2.);
  double cl = celerity(al, k0);

  double ar = a0;
  
  double am = 3.459578046858399;
  double cm = celerity(am, k0);
  double um = 9.192473939896399;
  double s = am*um/(am-ar);

  double xa = -cl*t;
  double xb = (um - cm)*t;
  double xc = s*t;

  if (x <= xa) return 0.;
  else if (xa < x && x <= xb) return  4./5.*x/t + 4./5.*cl;
  else if (xb < x && x <= xc) return um;
  else return 0.;	     
}

int main() {

  origin (-L/2., 0.);
  L0 = L;

  for (N = 128; N <= 1024; N *= 2) {

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

q[left] = neumann(0.);
q[right] = neumann(0.);

event defaults (i=0) {

  gradient = order1;
  riemann = hll_glu;
}

event init (i=0) {
  
  foreach() {
    k[] = K0;
    a0[] = pi*pow(R0, 2.);
    a[] = a0[]*(x < XM ? pow(1 + DR, 2.) : 1);
    q[] = 0.;
  }
}

event field (t = {0., 0.01, 0.02, 0.03, 0.04}) {

  if (N == 1024) {
    foreach() {
      fprintf (stderr, "%g, %.6f, %.6f\n", x, a[], q[]) ;
    }
    fprintf (stderr, "\n\n") ;
  }
}

event error (i++) {

  ne++;
  
  foreach() {
    ea[] = a[] - analyticA (t, x, DR, a0[], k[]);
    eq[] = q[] - analyticU (t, x, DR, a0[], k[])*
      analyticA (t, x, DR, a0[], k[]);
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

event end (t = 4.e-2) {
  printf ("#Tourniquet test case completed\n") ;
}

/**
# Plots
*/

/**
~~~gnuplot Spatial evolution of the cross-sectional area $A$ computed at $t={0, 0.01, 0.02, 0.03, 0.04}$ for $N=1024$ 

reset

mycolors = "dark-blue red sea-green dark-violet orange"
mypoints = '1 2 3 4 6'

set datafile separator ','
set style line 1 lt -1 lw 2 lc 'black'

k0 = 1.e4
a0 = pi
dr = 0.1

al = a0*(1 + dr)**2.
cl = sqrt(0.5*k0*sqrt(al))
ar = a0
am = 3.459578046858399;
cm = sqrt(0.5*k0*sqrt(am))
um = 9.192473939896399;
s = 100.01113797047884;
xa(t) = -cl*t;
xb(t) = (um - cm)*t;
xc(t) = s*t;
a(t, x) = x <= xa(t) ? al : \
     	  xa(t) < x && x <= xb(t) ? ((2./k0)**2.)*(-1./5.*x/t + 4./5.*cl)**4. : \
	  xb(t) < x && x <= xc(t) ? am : ar
u(t, x) = x <= xa(t) ? 0. : \
     	  xa(t) < x && x <= xb(t) ? 4./5.*x/t + 4./5.*cl : \
	  xb(t) < x && x <= xc(t) ? um : 0.
q(t, x) = a(t, x)*u(t, x)

set output 'A.png'
set xlabel 'x [cm]'
set ylabel 'A [cm^2]'
set key top right

plot a(0., x) w l ls 1 t 'Analytic', \
     a(0.01, x) w l ls 1 notitle, \
     a(0.02, x) w l ls 1 notitle, \
     a(0.03, x) w l ls 1 notitle, \
     a(0.04, x) w l ls 1 notitle, \
     'log' i 0 u 1:2 every 5 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     'log' i 1 u 1:2 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.01', \
     'log' i 2 u 1:2 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.02', \
     'log' i 3 u 1:2 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.03', \
     'log' i 4 u 1:2 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.04'
~~~

~~~gnuplot Spatial evolution of the flow rate $Q$ computed at $t={0, 0.01, 0.02, 0.03, 0.04}$ for $N=1024$ 

set output 'Q.png'
set xlabel 'x [cm]'
set ylabel 'Q [cm^3/s]'
set key center center

plot q(0., x) w l ls 1 t 'Analytic', \
     q(0.01, x) w l ls 1 notitle, \
     q(0.02, x) w l ls 1 notitle, \
     q(0.03, x) w l ls 1 notitle, \
     q(0.04, x) w l ls 1 notitle, \
     'log' i 0 u 1:3 every 5 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     'log' i 1 u 1:3 every 5 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.01', \
     'log' i 2 u 1:3 every 5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.02', \
     'log' i 3 u 1:3 every 5 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.03', \
     'log' i 4 u 1:3 every 5 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.04'
~~~

~~~gnuplot Comparison of the evolution of the error norms for the cross-sectional area $A$ with the number of cells $N$

set xtics 128,2,1024
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
set ylabel 'Error norms'
set key top right

set cbrange[1e-10:1e10]

plot '< grep ^ea out' u 2:3 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|A|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^ea out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|A|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^ea out' u 2:5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|A|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~

~~~gnuplot Comparison of the evolution of the error norms for the flow rate $Q$ with the number of cells $N$

fit f1(x) '< grep ^eq out' u (log($2)):(log($3)) via a1, b1
fit f2(x) '< grep ^eq out' u (log($2)):(log($4)) via a2, b2
fit f3(x) '< grep ^eq out' u (log($2)):(log($5)) via a3, b3

set output 'errQ.png'
set xlabel 'N'
set ylabel 'Error norms'
set key top right

plot '< grep ^eq out' u 2:3 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|Q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^eq out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|Q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^eq out' u 2:5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|Q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~
*/
