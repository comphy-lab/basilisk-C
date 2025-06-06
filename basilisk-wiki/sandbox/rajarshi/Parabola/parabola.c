/**
#Parabola - Saint Venant
*/

#include "grid/multigrid1D.h"
#include "saint-venant-implicit-raj.h"

double e1 = 0., e2 = 0., emax = 0.;
int ne = 0;

double A = 3000.;
double h0 = 10.;
double tau = 1e-3;
double Bf = 5.;

int main()
{
  origin (-5000);
  size (10000);
  G = 9.81;
  DT = 10.;
  dry = 1e-6;
  CFLa = 0.25;

  for (N = 32; N <= 256; N *= 2) {
    e1 = e2 = emax = 0.;
    ne = 0;
    run();
    fprintf (stderr, "%d %g %g %g\n", N, e1/ne/h0, sqrt(e2/ne)/h0, emax/h0);
  }
}

double Psi (double x, double t)
{
  // Analytical solution, see Sampson, Easton, Singh, 2006
  double p = sqrt (8.*G*h0)/A;
  double s = sqrt (p*p - tau*tau)/2.;
  return A*A*Bf*Bf*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + 
	    (tau*tau/4. - s*s)*cos (2.*s*t)) - Bf*Bf*exp(-tau*t)/(4.*G) -
    exp (-tau*t/2.)/G*(Bf*s*cos (s*t) + tau*Bf/2.*sin (s*t))*x;
}

event init (i = 0)
{
  foreach() {
    zb[] = h0*sq(x/A);
    h[] = max(h0 + Psi(x,0) - zb[], 0.);
  }
}


event friction (i++) {

  foreach()
    q.x[] /= 1. + tau*dt;
  boundary ({q.x});

}


scalar e[];

event error (i++) {
  foreach()
    e[] = h[] - max(h0 + Psi(x,t) - zb[], 0.);
  norm n = normf (e);
  e1 += n.avg;
  e2 += n.rms*n.rms;
  ne++;
  if (n.max > emax)
    emax = n.max;
  printf ("e %g %g %g %g %g\n", t, n.avg, n.rms, n.max, dt);
}

event field (t = 1500) {
  if (N == 64) {
    foreach()
      printf ("p %g %g %g %g %g\n", x, h[], q.x[], zb[], e[]);
    printf ("p\n");
  }
}

event umean (t += 50; t <= 6000) {
  if (N == 128) {
    double sq = 0., sh = 0.;
    foreach() {
      sq += Delta*q.x[];
      sh += Delta*h[];
    }
    printf ("s %g %g %f\n", t, sq/sh, sh);
  }
}

/**
~~~gnuplot Free surface and topography at $t=1500$ for $N=64$ grid points.
h0 = 10.
a = 3000.
tau = 1e-3
B = 5.
G = 9.81
p = sqrt (8.*G*h0)/a
s = sqrt (p*p - tau*tau)/2.
u0(t) = B*exp (-tau*t/2.)*sin (s*t)

set xlabel 'x (m)'
set ylabel 'z (m)'
t = 1500
psi(x) = a*a*B*B*exp (-tau*t)/(8.*G*G*h0)*(- s*tau*sin (2.*s*t) + \
      (tau*tau/4. - s*s)*cos (2.*s*t)) - B*B*exp(-tau*t)/(4.*G) - \
      exp (-tau*t/2.)/G*(B*s*cos (s*t) + tau*B/2.*sin (s*t))*x + h0
bed(x) = h0*(x/a)**2
set key top center
plot [-5000:5000] \
      '< grep ^p out' u 2:5:($5+$3) w filledcu lc 3 t 'Implicit', \
      psi(x) > bed(x) ? psi(x) : bed(x) lc 2 t 'Analytical', \
      bed(x) lw 3 lc 1 lt 1 t 'Bed profile'
~~~

~~~gnuplot Evolution of the axial velocity with time
reset
set key top right
set ylabel 'u0'
set xlabel 'Time'
plot u0(x) t 'analytical', \
     '< grep ^s out' u 2:3 every 2 w p t 'Implicit', \
     '< grep ^s ../parabola-explicit/out' u 2:3 every 2 w p t 'Explicit'
~~~

~~~gnuplot Convergence of the error on the free surface position
reset
set xlabel 'Resolution'
set ylabel 'Relative error norms'
set key bottom left
set logscale
set cbrange [1:2]
set xtics 32,2,512
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f1(x)=a1+b1*x
fit f1(x) 'log' u (log($1)):(log($2)) via a1,b1
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' u (log($1)):(log($4)) via am,bm
plot exp (f1(log(x))) t ftitle(a1,b1), \
     exp (f2(log(x))) t ftitle(a2,b2), \
     exp (fm(log(x))) t ftitle(am,bm),  \
     'log' u 1:2 t '|h|_1' ps 1.5, \
     'log' u 1:3 t '|h|_2' ps 1.5, \
     'log' u 1:4 t '|h|_{max}' ps 1.5 lc 0, \
     '../parabola-explicit/log' u 1:2 t '|h|_1 (explicit)' ps 1.5, \
     '../parabola-explicit/log' u 1:3 t '|h|_2 (explicit)' ps 1.5, \
     '../parabola-explicit/log' u 1:4 t '|h|_{max} (explicit)' ps 1.5
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/parabola.html)
*/
