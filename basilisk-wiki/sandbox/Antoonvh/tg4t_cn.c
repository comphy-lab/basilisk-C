/**
# Taylor-Green vortices

Vortex stability test with CN time-advancement

![vorticity](tg4t_cn/o.mp4)


~~~gnuplot Error
set xlabel 't'
set ylabel 'Error'
set logscale y
set yrange [1e-7:10]
set key top left
plot 'log' u 1:3 w l t 'N = 64',  \
     2.6e-11*exp(0.82*x),         \
     'out' u 1:3 w l t 'N = 128', \
     1e-13*exp(0.88*x)
~~~
 */
#include "nsf4t_cn.h"

double muv = 0;

double tend = 100;

int main(int argc, char ** argv) {
  if (argc > 1)
    tend = atof (argv[1]);
  X0 = Y0 = -L0/2;
  foreach_dimension()
    periodic (right);
  // for (N = 8; N <= 64; N *= 2)
  N = 64;
  run();

  N = 128;
  run();

  muv = 0.005;
  const scalar muc[] = muv;
  nu = kappa = muc;
  for (N = 8; N <= 64; N *= 2)
    //  run();
    ;
}

double u_x (double x, double y) {
  return -cos(2.*pi*x)*sin(2.*pi*y)*exp(-2.*muv*sq(2.*pi)*t) ;
}

double u_y (double x, double y) {
  return  sin(2.*pi*x)*cos(2.*pi*y)*exp(-2.*muv*sq(2.*pi)*t) ;
}

double s_a (double x, double y) {
  return  cos(2.*pi*x)*cos(2.*pi*y)*exp(-2*muv*sq(2.*pi)*t);
}

event init (t = 0) {
  TOLERANCE = 1e-5;
  foreach_face() 
    u.x[] = Gauss6_x (x, y, Delta, u_x);
  
}

event mov (i += 5) {
  if (N == 64) {
    scalar omg[];
    vorticityf (u, omg);
    output_ppm (omg, file = "o.mp4", n = 300, min = -10, max = 10);
  }
}

event errors (t += 1) {
  if (i % 5 != 0)
    event ("mov");
  double e = 0, es = 0, E = 0;
  double em = -1, esm = -1;
  foreach_face (serial) {
    E  += sq(Delta)*sq(u.x[]);
    double el = fabs(Gauss6_x (x, y, Delta, u_x) - u.x[]);
    e  += sq(Delta)*el;
    if (el > em)
      em = el;
  }
  fprintf (N == 64 ? stderr : stdout,
	   "%g %g %g %g %g %g %g\n", t, sqrt(grid->tn), e, em, E, es, esm);
  fflush (stdout);
  if (e > 1) 
    return 1;
  if (N == 64 && muv == 0) {
    scalar omg[];
    vorticityf (u, omg);
    output_ppm (omg, n = 400, file = "o.png", min = -4, max = 4);
  }
}

event stop (t = tend) {
  return 1;
}
