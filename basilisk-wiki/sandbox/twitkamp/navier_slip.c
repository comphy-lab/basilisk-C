
/**
# Couette flow with Navier slip boundary condition

A Navier slip boundary condition is implemented following (Afkhami 2009).^[Afkhami, S., Zaleski, S., & Bussmann, M. (2009). A mesh-dependent model for applying dynamic contact angles to VOF simulations. Journal of computational physics, 228(15), 5370-5389.] and the [robin boundary condition from Antoon](/sandbox/Antoonvh/robin.c). 

Consider a simple couette flow in a channel of height 2H and periodic in the $x$ direction, with the top boundary moving with velocity $u_x|_{y=H} = U$ and the bottom boundary being a Navier slip boundary condition (NSBC) of the form 

$$u_x|_{y=-H} = \lambda\frac{\partial{u_x}}{\partial y}.$$

The vertical flow profile $u_x(y)$ can be analytically computed with these boundary conditions and is given by, 

$$u_x(y) = \frac{U(y - H)}{2H + \lambda} + U$$.

In the limit $\lambda \rightarrow 0$, $u_x(-H) = 0$ which is the expect result for a no-slip boundary condition. If $k \rightarrow \infty$, $u_x(-H) = U$ which is the free slip. As a sanity check, the solution is bounded between 0 and $U$ which is expected.

# Numerical implementation of NSBC 

To enforce the Navier-slip boundary condition, the ghost cell velocity is computed using a linear interpolation, knowing the velocity $v_1$ at the first grid cell ($\Delta /2$ away from the wall) and the slip length $\lambda$. 

![Free-slip, Navier-slip and No-slip boundary conditions (Taken from Afkhami 2009)](NavierSlipAfkhami2009.png)

Let $v$ be the velocity at the wall $x = 0$, then the velocity $v(x)$ is given by

$$v(x) = v_a + \frac{v_b - v_a}{x_b - x_a} (x - x_a).$$ 

Two choices are possible namely, 

$$\{v_a, v_b, x_a, x_b\} = \{v, v_1, 0, \Delta\} \hspace{1cm} \text{or} \hspace{1cm} \{v_g, v, -\Delta, 0\},$$

giving 

$$v(x) = v - \frac{v_1 - v}{\Delta/2}x \hspace{1cm} \text{or} \hspace{1cm} v(x) = v_g + \frac{v - v_g}{\Delta/2}(x - \Delta/2).$$

Knowing that at $x = -\lambda$, $v(-\lambda) = 0$, $v_g$ can be solved as a function of $v_1$ as 
$$v_g = \frac{2\lambda - \Delta}{2\lambda + \Delta} v_1.$$

As a sanity check, $\lambda \rightarrow 0$ gives $v_g = -v_1$, i.e. Dirichlet boundary condition. For $\lambda \rightarrow \infty$, $v_g = v_1$ which is the neumann condition. 

Finally, numerically, the ghost cell is set to neumann(0.) $\frac{2\lambda - \Delta}{2\lambda + \Delta}$ since neumann(0.) sets $v_g = v_1$.
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

#define navier(slip) neumann(0.) * (2*slip - Delta)/(2*slip + Delta)
double slip_length;


u.t[top] = dirichlet(1.);
u.n[top] = dirichlet(0.);

u.t[bottom] = navier(slip_length);

int main()
{  
    
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);
  
  stokes = true;
  mu = fm;

  TOLERANCE = 1e-4;

  periodic(right);

  /** We Vary the domain size and slip length */
  for (N = 32; N <= 512; N *= 2){
    for (L0 = 0.5; L0 <= 4.; L0 += 0.5){
      for (double slip = 0.2; slip <= 0.4; slip += 0.2){
        slip_length = slip;
        run();
      }
    }
  }
}

/**
We look for a stationary solution. */

scalar un[];
event logfile (i++; i <= 20) {
  double du = change (u.x, un);
  fprintf(fout, "change %g %d %g \n", t, i, du);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/** We store the minimum velocity which is the velocity at the bottom boundary (Slip velocity) as a function of domain size 

The theoretical expression for the slip velocity, given $U = 1$ and $2H = L0$, is $u_x(y=-L0/2) = \frac{\lambda}{L0 + \lambda}$
*/
#define slip_vel(slip, L0) (slip/(L0 + slip))

event profile (t = end)
{
  
  if (N ==512){
    char file_slip[80];
    sprintf(file_slip, "slip_%g.dat", slip_length);
    FILE *fp_slip = fopen(file_slip, "a");
    fprintf(fp_slip, "%g %g\n", L0, statsf(u.x).min);
  }

  if (slip_length == 0.4 && L0 == 4.){
    double emax = 0.;
    double e = statsf(u.x).min - slip_vel (slip_length, L0);
    emax = max(e, emax);
    FILE *fp_err = fopen("error.dat", "a");
    fprintf(fp_err, "%d %g\n", N, emax);
  }
}


/**
## Results
~~~gnuplot Slip velocity as a function of Domain size
set xlabel 'L0'
set ylabel 'u_x|y = 0'
set grid
f(x) = 0.2 / (x + 0.2)
plot 'slip_0.2.dat' u 1:2 ps 2 pt 4 title 'Sim. data slip length = 0.2', \
    f(x) with lines lw 2 lc rgb "black" title 'Analytical slip = 0.2'
~~~
~~~gnuplot Slip velocity as a function of Domain size
set xlabel 'L0'
set ylabel 'u_x|y = 0'
set grid
f(x) = 0.4 / (x + 0.4)
plot 'slip_0.4.dat' u 1:2 ps 2 pt 4 title 'Sim. data slip length = 0.4', \
    f(x) with lines lw 2 lc rgb "black" title 'Analytical slip = 0.4'
~~~

~~~gnuplot Error as a function of gridsize for L0 = 4. and slip_length = 0.4
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
