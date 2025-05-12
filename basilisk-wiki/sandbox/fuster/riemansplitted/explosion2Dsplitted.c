/**
# Two-dimensional explosion with splitted Riemann solver

We solve the Euler equations for a compressible gas.*/

#include "vofsplitriemannsolver.h"

#if dimension == 2
# define LEVEL 7
#else // 3D
# define LEVEL 6
#endif

double gammao = 1.4;
double dtmax = 0.005;

foreach_dimension () 
void flux_x (const double * s, double * f, double e[2])
{

  /**
  We first recover each value ($\rho$, $E$, $w_x$ and $w_y$) and then
  compute the corresponding fluxes (`f[0]`, `f[1]`, `f[2]` and
  `f[3]`). */

  double rho = s[0], E = s[1], wn = s[2 + sweep_dir], w2 = 0.;
  for (int i = 2; i < 2 + dimension; i++)
    w2 += sq(s[i]);
  double un = wn/rho;

  double p = (E - 0.5*w2/rho)*(gammao - 1.);
  
  f[0] = wn;
  f[1] = un*(E + p);
  for (int i = 2; i <= 2 + dimension; i++)
    f[i] = un*s[i];
  f[2+sweep_dir] += p;

  /**
  The minimum and maximum eigenvalues for the Euler system are the
  characteristic speeds $u \pm \sqrt(\gamma p / \rho)$. */

  double c = sqrt(gammao*p/rho);
  e[0] = un - c; // min
  e[1] = un + c; // max
}


int main() {

  /**
  We make boundary conditions free outflow. */

  foreach_dimension() {
    q.n[right] = neumann(0);
    q.n[left]  = neumann(0);
  }
  
  /**
  The domain spans $[-1:1]\times[-1:1]\times[-1:1]$. */

  origin (-1, -1, -1);
  size (2.);
  init_grid (1 << LEVEL);
  run(); 
}

/**
Initial conditions come from Toro's book (Riemann Solvers and
Numerical Methods for Fluid Dynamics, 3rd Edition, Springer Ed.)
Chapter 17 section 17.1.1 are given in terms of density ($\rho$),
pression ($p$), velocity ($u$) both at the left and right side of the
discontinuity placed at $R=0.4$. */

event init (t = 0)
{
  double R = 0.4 ;
  double rhoL = 1., rhoR = 0.125 ;
  double pL = 1.0,  pR = 0.1 ;
  
  /**
  To define an accurate (i.e. symmetrical) initial sphere of rayon
  $R$, we compute the volume fraction corresponding to the spherical
  interface. */

  fraction (c, sq(x) + sq(y) + sq(z) - sq(R));
  
  /**
  Left and right initial states for $\rho$, $\mathbf{w}$ and energy
  $E = \rho \mathbf{u}^2/2 + p/(\gamma-1)$. */
  
  foreach() {
    rho[] = rhoR*c[] + rhoL*(1. - c[]);
    foreach_dimension()
      q.x[] = 0.;
    Etot[] = (pR*c[] + pL*(1. - c[]))/(gammao - 1.);
  }
}

event print (t = 0.25)
{

  /**
  At $t=0.25$ we output the values of $\rho$ and the normal velocity
  $\mathbf{u}_n$ as functions of the radial coordinate. */

  foreach() {
    double r = sqrt(sq(x) + sq(y) + sq(z));
    double wn = (q.x[]*x + q.y[]*y + q.z[]*z)/r;
    printf ("%g %g %g\n", r, rho[], wn/rho[]);
  }

}

/**
## Results

Results are presented in terms of $\rho$ and normal velocity $u_n$ for
the splitted and non-splitted solver. 

~~~gnuplot Radial density profile
set xrange [0:1]
set xlabel 'r'

set output 'rho.png'
set ylabel 'rho'
plot './out' u 1:2 w p pt 7 ps 0.2 t 'splitted', \
     './../../../../src/test/explosion/cout' u 1:2 w p pt 7 ps 0.2 t 'unsplitted'
~~~

~~~gnuplot Normal velocity
set output 'velocity.png'
set ylabel 'Normal velocity'
plot './out' u 1:3 w p pt 7 ps 0.2 t 'splitted',		  \
     './../../../../src/test/explosion/cout' u 1:3 w p pt 7 ps 0.2 t 'unsplitted'
~~~

*/