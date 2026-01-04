/**
# Couette flow with Robin boundary condition

A robin boundary condition is implemented following the [robin boundary condition from Antoon](/sandbox/Antoonvh/robin.c). 

Consider a simple couette flow in a channel of height 2H and periodic in the $x$ direction, with the top boundary moving with velocity $u_x|_{y=H} = U_t$ and the bottom boundary being a robin boundary condition of the form 

$$u_x|_{y=-H} = U_b + \lambda\frac{\partial{u_x}}{\partial y}.$$

The vertical flow profile $u_x(y)$ can be analytically computed with these boundary conditions and is given by, 

$$u_x(y) = \frac{(U_t - U_b)(y - H)}{2H + \lambda} + U_t$$.

In the limit $\lambda \rightarrow 0$, $u_x(-H) = U_b$ which is the expect result for a no-slip boundary condition. If $k \rightarrow \infty$, $u_x(-H) = U_t$ which is the bottom wall velocity. As a sanity check, the solution is bounded between $U_b$ and $U_t$ which is expected.

# Numerical implementation of robin boundary condition 

see https://basilisk.fr/sandbox/Antoonvh/robin.c

*/


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define robin(slip, vel) ((dirichlet ((vel)*Delta/(2*slip + Delta))) + ((neumann (0))*((2*(slip) - Delta)/(2*(slip) + Delta) + 1.)))

double slip_length, vel_bottom;

u.t[top] = dirichlet(1.);
u.n[top] = dirichlet(0.);

u.t[bottom] = robin(slip_length, vel_bottom);

int main()
{  
    
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);
  
  stokes = true;
  mu = fm;

  TOLERANCE = 1e-4;

  periodic(right);

  // N = 512;
  /** We Vary the domain size and slip length */
  for (N = 32; N <= 512; N *= 2){
    for (L0 = 0.5; L0 <= 4.; L0 += 0.5){
      for (double vel = 0.2; vel <= 0.4; vel += 0.2){
        slip_length = 0.4;
        vel_bottom = vel;
        run();
      }
    }
  }
}

/**
We look for a stationary solution. */

scalar un[];
event logfile (i++; i <= 1000) {
  double du = change (u.x, un);
  fprintf(fout, "change %g %d %g \n", t, i, du);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/** We store the minimum velocity which is the velocity at the bottom boundary (Slip velocity) as a function of domain size */
#define robin_vel(slip, L0, vel_bottom) ((slip + vel_bottom * L0) /(L0 + slip))

event profile (t = end)
{
  
  if (N ==128){
    char file_slip[80];
    sprintf(file_slip, "vel_%g.dat", vel_bottom);
    FILE *fp_slip = fopen(file_slip, "a");
    fprintf(fp_slip, "%g %g\n", L0, statsf(u.x).min);
  }

  if (vel_bottom == 0.4 && L0 == 4.){
    double emax = 0.;
    double u_anal = robin_vel (slip_length, L0, vel_bottom);
    double e = statsf(u.x).min - u_anal ;
    fprintf(ferr, "e %g uanal %g\n", e, u_anal);
    emax = max(e, emax);
    FILE *fp_err = fopen("error.dat", "a");
    fprintf(fp_err, "%d %g\n", N, emax);
  }
}


/**
# Results
~~~gnuplot Slip velocity as a function of Domain size for $\lambda = 0.4, U_b = 0.2$ and $U_t = 1$ 
set xlabel 'L0'
set ylabel 'u_x|y = 0'
set grid
f(x) = (0.4 + 0.2*x) / (x + 0.4)
plot 'vel_0.2.dat' u 1:2 ps 2 pt 4 title 'Simulation Data', \
    f(x) with lines lw 2 lc rgb "black" title 'Analytical slip = 0.4, U_b = 0.2'
~~~

~~~gnuplot Slip velocity as a function of Domain size for $\lambda = 0.4, U_b = 0.4$ and $U_t = 1$ 
set xlabel 'L0'
set ylabel 'u_x|y = 0'
set grid
f(x) = (0.4 + 0.4*x) / (x + 0.4)
plot 'vel_0.4.dat' u 1:2 ps 2 pt 4 title 'Simulation Data', \
    f(x) with lines lw 2 lc rgb "black" title 'Analytical slip = 0.4, U_b = 0.4'
~~~

~~~gnuplot Error as a function of Resolution for $\lambda = 0.4, U_b = 0.4, U_t = 1$ and $L_0 = 4$  
unset arrow
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) 'error.dat' u (log($1)):(log($2)) via a,b
set xlabel 'Resolution'
set logscale
set xtics 8,2,1024
set ytics format "% .0e"
set grid ytics
set cbrange [1:2]
set xrange [16:1024]
set ylabel 'Error'
set yrange [*:*]
set key top right
plot '' u 1:2 pt 6 t 'max', exp(f(log(x))) t ftitle(a,b)
~~~
*/