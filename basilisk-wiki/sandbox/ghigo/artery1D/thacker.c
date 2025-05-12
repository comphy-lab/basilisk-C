/**
# Thacker solution

We solve, using the 1D blood flow equations, the inviscid pressure
oscillations in a parabolic aneurysm. This test case is greatly
inspired by the solutions of Thacker and Sampson, well-known in the
shallow water community. */

/**
# Analytic solution

We choose the following spatial variations of the neutral cross-sectional area $A_0$ and the arterial wall rigidity $K^\prime$ in the domain $x\in \left[-a,\, a \right]$:

$$
\left\{
\begin{aligned}
& R_0 = \bar{R_0}\left[1 + \Delta \mathcal{G} \left[1 - \frac{x^2}{a^2} \right] \right] && \text{ with } \bar{R_0} >0 ,\, a >0 ,\, \Delta \mathcal{G} > -1 \\
& K^\prime = \mathrm{cst} && \text{ with } K^\prime >0,
\end{aligned}
\right.
$$

where $\Delta \mathcal{G} > 0$. Simple analysis allows us to obtain the following analytic solution:

$$
\left\{
\begin{aligned}
& U_0 =  U \sin \left( \frac{t}{\tau}  \right)  && \text{ with } U = \mathrm{cst}\\
& p = - \frac{1}{4} \rho U^2 \cos \left( 2 \frac{t}{\tau} \right) 
 - \rho \frac{x}{\tau} U \cos \left( \frac{t}{\tau} \right),
\end{aligned}
\right.
$$

where: 

$$
\tau = \lvert \Delta \mathcal{G} \omega^2 \rvert^{-\frac{1}{2}} .
$$

Finally, we choose $U$ with respect to the following nonlinear stability arguments, namely that the radius $R$ remains positive and that the flow remains subcritical. */

/**
# Code
*/

#include "grid/cartesian1D.h"
#include "bloodflow.h"

double R0 = 1.;
double K0 = 1.e4;

double BTH = 4.;
double DELTATH = 0.1;
double ATH = 1.;

double KTH = 0.;
double CTH = 0.;
double WTH = 0.;
double TAUTH = 0.;
double TTH = 0.;

scalar eq[];
double eq1, eq2, eqmax;
int ne;

int ihr = 0;

double celerity (double a, double k) {

  return sqrt(0.5*k*sqrt(a));
}

double r0thacker(double x, double r0, double b, double delta) {
  
  return r0*(1. + delta*(1. - pow(x/b, 2.)));
}

double pthacker(double x, double t, double r0, double rho,
		double tau, double ath) {
  
  return -0.25*rho*pow(ath, 2.)*cos(2.*t/tau) - rho*ath*x/tau*cos(t/tau);
}

double rthacker(double x, double t, double r0, double b,
		double delta, double k, double rho, double tau, double ath) {
  
  return 1./k*pthacker(x, t, r0, rho, tau, ath) + r0thacker(x, r0, b, delta);
}

double athacker(double x, double t, double r0, double b,
		double delta, double k, double rho, double tau, double ath) {
  
  return pi*pow(rthacker(x, t, r0, b, delta, k, rho, tau, ath), 2.);
}

double uthacker(double t, double tau, double ath) {

  return ath*sin(t/tau);
}

double qthacker(double x, double t, double r0, double b,
		double delta, double k, double rho, double tau, double ath) {
  
  return uthacker(t, tau, ath)*athacker(x, t, r0, b, delta, k, rho, tau, ath);
}

int main() {

  origin (-BTH, 0.);
  L0 = 2.*BTH;

  KTH = K0*sqrt(pi);
  CTH = sqrt(0.5*KTH*R0);
  WTH = 2.*CTH/BTH;
  TAUTH = 1./sqrt(abs(DELTATH*pow(WTH, 2.)));
  TTH = 2.*pi/sqrt(abs(DELTATH*pow(WTH, 2.)));

  ihr = 0;
  for (ihr = 0; ihr <= 2; ihr ++) {
    for (N = 32; N <= 256; N *= 2) {

      eq1 = eq2 = eqmax = 0.;
      ne = 0;
      
      run();

      printf ("eq%d, %d, %g, %g, %g\n", ihr, N, eq1/ne,
	      sqrt(eq2/ne), eqmax);
    }
  }
  
  return 0; 
}

q[left] = dirichlet (qthacker (-BTH, t, R0, BTH, DELTATH, KTH, 1., TAUTH, ATH));
q[right] = dirichlet (qthacker (BTH, t, R0, BTH, DELTATH, KTH, 1., TAUTH, ATH));

event defaults (i=0) {

  gradient = order1;
  if (ihr == 0) riemann = hll_hr;
  else if (ihr == 1) riemann = hll_hrls;
  else if (ihr == 2) riemann = hll_glu;
}

event init (i=0) {
  
  foreach() {
    k[] = K0;
    a0[] = pi*pow(r0thacker (x, R0, BTH, DELTATH), 2.);
    a[] = athacker (x, 0., R0, BTH, DELTATH, KTH, 1., TAUTH, ATH);
    q[] = qthacker (x, 0., R0, BTH, DELTATH, KTH, 1., TAUTH, ATH);
  }
}

event field (t = {0., 0.1*0.422653, 0.3*0.422653, 0.6*0.422653, 0.8*0.422653}) {

  if (N == 128) {
    foreach() {
      fprintf (stderr, "%d, %g, %.6f, %.6f, %.6f, %.6f\n", ihr, x,
	       a0[]/(pi*R0*R0), k[]/K0, a[]-a0[], q[]);
    }
    fprintf (stderr, "\n\n");
  }
}

event error (i++) {

  ne++;
  
  foreach() {
    eq[] = q[] - qthacker (x, t, R0, BTH, DELTATH,
				KTH, 1., TAUTH, ATH);
  }
    
  norm nq = normf (eq);
  eq1 += nq.avg;
  eq2 += nq.rms*nq.rms;
  if (nq.max > eqmax)
    eqmax = nq.max;
}

event end (t = 1.) {
  printf ("#Thacker test case completed\n");
}

/**
# Plots
*/

/**
~~~gnuplot Spatial evolution of the cross-sectional area at rest $A_0$ and arterial wall rigidity $K$ 

reset

mycolors = "dark-blue red sea-green dark-violet orange brown4"
mypoints = '1 2 3 4 6'

set datafile separator ','
set style line 1 lt -1 lw 2 lc 'black'

set output 'A0K.png'
set xlabel 'x [cm]'
set ylabel 'Dimensionless arterial properties'
set key top right

plot '< grep -e "^0" -e "^$" log' i 0 u 2:3 w l ls 1 lc rgb word(mycolors, 1) t 'A_0/A_0(0)', \
     '< grep -e "^0" -e "^$" log' i 0 u 2:4 w l ls 1 lc rgb word(mycolors, 2) t 'K/K(0)'
~~~

~~~gnuplot Spatial evolution of the cross-sectional area $A-A_0$ computed at $t=\{0,\, 0.1,\, 0.3,\, 0.6,\, 0.8\}T$ for $N=128$ using the HR well-balanced scheme

r0 = 1.
k0 = 1.e4

bth = 4.
deltath = 0.1
ath = 1.

kth = k0*sqrt(pi)
cth = sqrt(0.5*kth*r0)
wth = 2.*cth/bth
tauth = 1./sqrt(abs(deltath*(wth)**2.))
tth = 2.*pi/sqrt(abs(deltath*(wth)**2.))

r0thacker(x) = r0*(1. + deltath*(1. - (x/bth)**2.))
a0thacker(x) = pi*r0thacker(x)**2.
pthacker(x, t) = -0.25*(ath)**2.*cos(2.*t/tauth) - ath*x/tauth*cos(t/tauth)
rthacker(x, t) = 1./kth*pthacker(x, t) + r0thacker(x)
athacker(x, t) = pi*rthacker(x, t)**2.

uthacker(t) = ath*sin(t/tauth)
qthacker(x, t) = uthacker(t)*athacker(x, t)

set output 'AmA0.png'
set xlabel 'x [cm]'
set ylabel 'A-A_0 [cm^2]'
set key top center

plot athacker(x, 0) - a0thacker(x) w l ls 1 t 'Analytic', \
     athacker(x, 0.1*0.422653) - a0thacker(x) w l ls 1 notitle, \
     athacker(x, 0.3*0.422653) - a0thacker(x) w l ls 1 notitle, \
     athacker(x, 0.6*0.422653) - a0thacker(x) w l ls 1 notitle, \
     athacker(x, 0.8*0.422653) - a0thacker(x) w l ls 1 notitle, \
     '< grep -e "^0" -e "^$" log' i 0 u 2:5 every 2 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     '< grep -e "^0" -e "^$" log' i 1 u 2:5 every 2 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.1T', \
     '< grep -e "^0" -e "^$" log' i 2 u 2:5 every 2 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.3T', \
     '< grep -e "^0" -e "^$" log' i 3 u 2:5 every 2 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.6T', \
     '< grep -e "^0" -e "^$" log' i 4 u 2:5 every 2 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.8T'
~~~

~~~gnuplot Spatial evolution of the flow rate $Q$ computed at $t=\{0,\, 0.1,\, 0.3,\, 0.6,\, 0.8\}T$ for $N=128$ using the HR well-balanced scheme

set output 'Q.png'
set xlabel 'x [cm]'
set ylabel 'Q [cm^3/s]'
set key center center

plot qthacker(x, 0) w l ls 1 t 'Analytic', \
     qthacker(x, 0.1*0.422653) w l ls 1 notitle, \
     qthacker(x, 0.3*0.422653) w l ls 1 notitle, \
     qthacker(x, 0.6*0.422653) w l ls 1 notitle, \
     qthacker(x, 0.8*0.422653) w l ls 1 notitle, \
     '< grep -e "^0" -e "^$" log' i 0 u 2:6 every 2 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0', \
     '< grep -e "^0" -e "^$" log' i 1 u 2:6 every 2 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=0.1T', \
     '< grep -e "^0" -e "^$" log' i 2 u 2:6 every 2 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=0.3T', \
     '< grep -e "^0" -e "^$" log' i 3 u 2:6 every 2 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=0.6T', \
     '< grep -e "^0" -e "^$" log' i 4 u 2:6 every 2 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=0.8T'
~~~

~~~gnuplot Comparison of the evolution of the $L_2$ relative error for the flow rate $|Q|_2$ with the number of cells $N$ computed using the HR, HR-LS and GLU well-balanced schemes.
 
set xtics 32,2,256
set logscale

ftitle(a,b) = sprintf('Order %4.2f', -b)
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x

fit f1(x) '< grep ^eq0 out' u (log($2)):(log($4)) via a1, b1
fit f2(x) '< grep ^eq1 out' u (log($2)):(log($4)) via a2, b2
fit f3(x) '< grep ^eq2 out' u (log($2)):(log($4)) via a3, b3

set output 'errQL2.png'
set xlabel 'N'
set ylabel '|Q|_2'
set key top right

set cbrange[1e-10:1e10]

plot '< grep ^eq0 out' u 2:4 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t 'HR, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^eq1 out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t 'HR LS, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^eq2 out' u 2:4 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t 'GLU, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~

~~~gnuplot Evolution of the relative error norms for the flow rate $Q$ with the number of cells $N$ computed using the HR well-balanced scheme

fit f1(x) '< grep ^eq0 out' u (log($2)):(log($3)) via a1, b1
fit f2(x) '< grep ^eq0 out' u (log($2)):(log($4)) via a2, b2
fit f3(x) '< grep ^eq0 out' u (log($2)):(log($5)) via a3, b3

set output 'errQHR.png'
set xlabel 'N'
set ylabel 'Relative error norms for HR'
set key top right

plot '< grep ^eq0 out' u 2:3 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|Q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 notitle, \
     '< grep ^eq0 out' u 2:4 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|Q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 notitle, \
     '< grep ^eq0 out' u 2:5 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|Q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 notitle
~~~
*/


