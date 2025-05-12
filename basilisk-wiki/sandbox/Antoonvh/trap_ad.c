/**
# Advection and the trapezoidal rule 

We test an implicit time integration scheme for the advection
equation, using simple 2nd-order accurate and colocated central finite
differences.

![Convergence test: Evolution of the tracer field](trap_ad/s.mp4)

~~~gnuplot Convergence rate
set logscale x 2
set logscale y 4
set grid
set xr [8:2048]
set yr [2e-2:50]
set size square
set xlabel 'N'
set ylabel 'Error'
plot 'out' u 1:3 t 'L_1 error', 1e4*x**(-2) t 'Second order'
~~~

We also plot the linear-solver statistics of the last
step. Convergence rates seem OK, even not to bad for the noisy
solutions.

~~~gnuplot linear-system solve convergence rate
reset
set logscale x 2
set yr [0:6]
set xr [8:2048]
set grid
set xlabel 'N'
plot 'out' u 1:4 t 'mg.i', 'out' u 1:5 t 'mg.nrelax'
~~~
*/

#include "run.h"
#include "solve.h"

scalar s[];
vector u[];

double x0 = 3;

int main() {
  L0 = x0*10;
  X0 = Y0 = -L0/2;
  CFL = 1.4; 
  TOLERANCE = 1e-6; 
  for (N = 16; N <= 1024; N *= 2) {
    DT = CFL* (L0/N) / (x0 + 3);
    run();
  }
}

#define Rad (sqrt(sq(x) + sq(y)))

event init (t = 0) {
  foreach() {
    s[] = exp(-sq(x - x0) - sq(y));
    u.x[] = -y*(Rad < (x0 + 3));
    u.y[] =  x*(Rad < (x0 + 3));
  } 
}
/**
## Implicit scheme

The 2nd-order accurate Adams-Moulton scheme is the trapezoidal rule, 

$$s_{n+1} = s_n + \frac{1}{2} \mathrm{d}t \left( f(s_n) + f(s_{n+1})\right).$$

To `solve()`:

$$s_{n+1} - \frac{1}{2}\mathrm{dt}f(s_{n+1}) = s_n + \frac{1}{2} \mathrm{d}t f(s_n).$$

with,

$$f(s) = -\sum_{\mathrm{dim}}\mathrm{u.x[\ ]} \frac{\mathrm{s[1]} - \mathrm{s[-1]}}{2\Delta}$$
*/

mgstats solve_stats;

event advection (i++, last) {
  dt = dtnext (DT);
  scalar st[];
  foreach() {
    st[] = s[];
    foreach_dimension() 
      st[] -= dt*u.x[]*(s[1] - s[-1])/(4*Delta); //rhs
  }
  foreach()
    s[] += 2*(st[] - s[]);                       //Guess
  solve_stats = 
    solve (s, s[] + 0.25*dt/Delta*(+ u.x[]*(s[1]   - s[-1])
                                   + u.y[]*(s[0,1] - s[0,-1])) , st[]);
}

event mov (i += 5) {
  output_ppm (s, file = "s.mp4", n = 300, min = -1, max = 1);
}

event stop (t = 2*pi) {
  double err = 0;
  foreach (reduction (+:err)) 
    err += sq(Delta)*fabs(s[] - exp(-sq(x - x0) - sq(y)));
  printf ("%d %d %g %d %d\n", N, i, err,  solve_stats.i, solve_stats.nrelax);
  return 1;
}

