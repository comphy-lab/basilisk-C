/**
# Inviscid wave propagation in a tapered artery

We solve the propgation of an inviscid wave in a tapered artery using
the inviscid 1D blood flow equations. */

/**
# Analytic solution

To capture the perturbations induced by vessel tapering, we introduce the following non-dimensional variables:

$$
t = T \bar{t},\, x = X \bar{x},\,R_0 = R_0 \bar{f},\,R = R_0 \left[ \bar{f} + \Delta_R \bar{R} \right],\,K = K_0 \bar{g},\,Q = Q \bar{Q},\,p = p_0,\,\Pi \tilde{p}.
$$

Injecting these non-dimensional variables in the 1D blood flow equations which we then linearize, we obtain the following simplified equation:

$$
\frac{1}{\bar{c}^2} \frac{\partial^2 \bar{Q}}{\partial \bar{t}^2} - \frac{\partial^2 \bar{Q} }{\partial \bar{x}^2 } = \left[ \frac{1}{\bar{g}}\frac{\partial \bar{g} }{\partial \bar{x} } - \frac{1}{\bar{f}}\frac{\partial \bar{f} }{\partial \bar{x} } \right] \frac{\partial \bar{Q} }{\partial \bar{x} },
$$
where $\bar{c} = \sqrt{\bar{f}\bar{g}}$ is the dimensionless speed.

We then search for a solution of the form:

$$
\bar{Q} = \tilde{Q}\left( \bar{x} \right)\exp\left( i \omega \bar{t} \right), \quad \omega \in \mathbb{R},
$$

and therefore rewrite the previous equation as:

$$
\frac{\mathrm{d}^2 \tilde{Q} }{\mathrm{d} \bar{x}^2 } + \frac{\omega^2}{\bar{c}^2} \tilde{Q} = - \left[ \frac{1}{\bar{g}}\frac{\mathrm{d} \bar{g} }{\mathrm{d} \bar{x} } - \frac{1}{\bar{f}}\frac{\mathrm{d} \bar{f} }{\mathrm{d} \bar{x} } \right] \frac{\mathrm{d} \tilde{Q} }{\mathrm{d} \bar{x} }.
$$

To keep track of the slowly varying neutral radius $R_0$ and arterial wall rigidity $K$, we use the following change of variables, to place ourselves at long $x$ while keeping track of local variations of the wave speed:

$$
\frac{d \xi}{d \bar{x}} = \Phi^{\prime}\left( X \right) , \quad X = \epsilon \bar{x},
$$

where $\epsilon$ is the small parameter characterizing the slow variations of the neutral radius $R_0$ and the arterial wall rigidity $K$. The function $\Phi^\prime$ represents the wave distortion. Using this change of variables, we have obtain the following solution at first order:

$$
\tilde{Q_0} = \frac{B}{\sqrt{\omega}} \frac{ \bar{f}^{\frac{3}{4}}}{\bar{g}^{\frac{1}{4}}} \exp \left( i\omega \left[ \bar{t} - \frac{1}{\epsilon} \int_{0}^{X} \frac{1}{c} \mathrm{d}X \right] \right) + C.C., \qquad B = \mathrm{cst}.
$$

*/
/**
# Code
*/

#include "grid/cartesian1D.h"
#include "bloodflow.h"

double R0 = 1.;
double K0 = 1.e4;
double L = 200.;

double T0 = 1.;

double DR0 = -0.1;
double DK0 = 0.1;
double DR = 0., DK = 0.;

double SH = 1.e-2;
double AIN = 0., QIN = 0.;

/* scalar eq[]; */
/* double eq1, eq2, eqmax; */
/* int ne = 0; */

int ihr = 0;

double celerity (double a, double k) {

  return sqrt(0.5*k*sqrt(a));
}

double shape (double x, double dr) {

  return 1 + dr*x;
}

double inlet (double t, double t0) {

  if (t < t0) return max(0., 0.5*(1. + cos(pi + 2.*pi/T0*t)));
  else return 0.;
}

int main() {

  origin (0., 0.);
  L0 = L;

  DR = DR0/T0/celerity(pi*R0*R0, K0);
  DK = DK0/T0/celerity(pi*R0*R0, K0);
  
  AIN = pi*pow(R0*(1 + SH), 2.);
  QIN = SH*AIN*celerity(AIN, K0);

  ihr = 0;
  for (ihr = 0; ihr <= 2; ihr ++) {
    for (N = 128; N <= 128; N *= 2) {

      /* eq1 = eq2 = eqmax = 0.; */
      /* ne = 0; */
    
      run();

      /* printf ("eq%d, %d, %g, %g, %g\n", ihr, N, eq1/QIN/ne, */
      /* 	      sqrt(eq2)/QIN/ne, eqmax/QIN); */
    }
  }
  
  return 0; 
}

q[left] = dirichlet (QIN*inlet (t, T0));

event defaults (i=0) {

  gradient = order1;
  if (ihr == 0) riemann = hll_hr;
  else if (ihr == 1) riemann = hll_hrls;
  else if (ihr == 2) riemann = hll_glu;
}

event init (i=0) {
  
  foreach() {
    k[] = K0*shape (x, DK);
    a0[] = pi*pow(R0*shape (x, DR), 2.);
    a[] = a0[];
    q[] = 0.;
  }
}

event field (t = {0.5, 1., 1.5, 2., 2.5}) {

  if (N == 128) {
    foreach() {
      fprintf (stderr, "%d, %g, %.6f, %.6f, %.6f, %.6f\n", ihr, x,
	       a0[]/(pi*R0*R0), k[]/K0, a[] - a0[], q[]) ;
    }
    fprintf (stderr, "\n\n") ;
  }
}

/* event error (i++) { */

/*   ne++; */
  
/*   foreach() { */
/*     eq[] = q[] - QIN; */
/*   } */
    
/*   norm nq = normf (eq); */
/*   eq1 += nq.avg; */
/*   eq2 += nq.rms*nq.rms; */
/*   if (nq.max > eqmax) */
/*     eqmax = nq.max; */
/* } */

event end (t = 2.5) {
  printf ("#Taper test case completed\n") ;
}

/**
# Plots
*/

/**
~~~gnuplot Spatial evolution of the cross-sectional area at rest $A_0$ and arterial wall rigidity $K$ 

reset

mycolors = "dark-blue red sea-green dark-violet orange"
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

~~~gnuplot Spatial evolution of the cross-sectional area $A-A_0$ computed at $t={0.5, 1, 1.5, 2, 2.5}$ for $N=128$ using the HR well-balanced scheme

set output 'AmA0.png'
set xlabel 'x [cm]'
set ylabel 'A-A_0 [cm^2]'
set key top right

plot '< grep -e "^0" -e "^$" log' i 0 u 2:5 every 2 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0.5', \
     '< grep -e "^0" -e "^$" log' i 1 u 2:5 every 2 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=1', \
     '< grep -e "^0" -e "^$" log' i 2 u 2:5 every 2 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=1.5', \
     '< grep -e "^0" -e "^$" log' i 3 u 2:5 every 2 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=2', \
     '< grep -e "^0" -e "^$" log' i 4 u 2:5 every 2 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=2.5'
~~~

~~~gnuplot Spatial evolution of the flow rate $Q$ computed at $t={0.5, 1, 1.5, 2, 2.5}$ for $N=128$ using the HR well-balanced scheme

set output 'Q.png'
set xlabel 'x [cm]'
set ylabel 'Q [cm^3/s]'
set key top right

plot '< grep -e "^0" -e "^$" log' i 0 u 2:6 every 2 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 2 lw 1 t 't=0.5', \
     '< grep -e "^0" -e "^$" log' i 1 u 2:6 every 2 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 2 lw 1 t 't=1', \
     '< grep -e "^0" -e "^$" log' i 2 u 2:6 every 2 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 2 lw 1 t 't=1.5', \
     '< grep -e "^0" -e "^$" log' i 3 u 2:6 every 2 w p lt word(mypoints, 4) lc rgb word(mycolors, 4) ps 2 lw 1 t 't=2', \
     '< grep -e "^0" -e "^$" log' i 4 u 2:6 every 2 w p lt word(mypoints, 5) lc rgb word(mycolors, 5) ps 2 lw 1 t 't=2.5'
~~~
*/
