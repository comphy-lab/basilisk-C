/**
# Finite difference Advection in 2D

![Test case](vad2D/s.mp4)

![Error distribution on equidistant grid](vad2D/er.png)

~~~gnuplot 
set xr[30: 300]
set yr[0.04:10]
set grid
set size square
set ylabel 'L_1 Error'
set xlabel 'N_{eff}'
set logscale xy 2
plot 'log' t 'Equidistant', 1e6*x**(-3) t '3rd order',\
 'out' t 'Adaptive', 3e6*x**(-4) t '4th order'
~~~
 */
#define RKORDER 4
#include "lsrk.h"
#include "my_vertex.h"

vector u[];
#define SOL(a) (exp(-sq(x - (a)) - sq(y)))
#define FUNC (SOL(2) - SOL(-4))

foreach_dimension() {
#define duc_x(uc) ((u.x[] > 0 ? (2*uc[1]     + 3*uc[] -6*uc[-1]     + uc[-2])/6.     : \
		   (-2.*uc[-1]     - 3.*uc[] + 6.*uc[1]     - uc[2])/6.)/Delta)
}
void adv (scalar * sl, scalar * dsl) {
  boundary (sl);
  foreach_vert() {
    scalar s, ds;
    for (s, ds in sl, dsl) {
      ds[] = 0;
      foreach_dimension()
	ds[] += -u.x[]*duc_x(s);
    }
  }
}

int maxlevel = 15;
vertex scalar s[];
double es = 0.001;
int j = 0;
int main() {
  foreach_dimension()
    periodic (left);
  L0 = 25;
  X0 = Y0 = -L0/2;
  N = 1 << 11;
  for (es = 0.01; es > 1e-4; es /= 2)
    run();
  j = 1;
  for (N = 32; N < 512; N *= 2)
    run();
}

event init (t = 0) {
  s.prolongation = s.refine = refine_vert5;
  s.coarsen = s.restriction = restriction_vert;
  foreach_vert() 
    s[] = FUNC;
  boundary ({s});
}

event advance (i++, last) {
  DT = HUGE;
  foreach_vert() {
    u.x[] = y;
    u.y[] = -x;
    double um = max(fabs(x), fabs(y));
    double dts = um > 0 ? Delta/um : HUGE;
    DT = min (DT, dts);
  }
  dt = dtnext (DT);
  A_Time_Step ({s}, dt, adv);
}

scalar w[];
event mov (t += pi/20) {
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (s, file = "s.mp4", n = 300);
  output_ppm (lev, file = "l.mp4", n = 300, min = 1, max = maxlevel);
  output_ppm (w, file = "w.mp4", n = 500);
}

#include "adapt_field.h"

event adapt (i++) {
  restriction ({s});
  boundary ({s});
  for (int l = 1; l <= depth(); l++) {
    foreach_vert_level(l) {
      w[] = 0;
      foreach_dimension()
	w[] += sq(fabs(u.x[])*fabs(s[] - interp4_x(point, s)));
      w[] = sqrt(w[]);
    }
  }
  if (!j)
    adapt_field (w, es, es/1.5, maxlevel, 5);
}

event stop (t = 2*pi) {
  double e = 0;
  scalar er[];
  foreach_vert() {
    er[] = fabs(s[] - FUNC);
    e += sq(Delta)*fabs(s[] - FUNC);
  }
  output_ppm (er, file = "er.png", n = 500);
  if (!j)
    printf ("%g %g\n", sqrt(grid->tn), e);
  else
    fprintf (stderr, "%g %g\n", sqrt(grid->tn), e);
    
}
