/**
# Bug: Poisson solver does not converge if it is compiled with mpi and run in serial

The error is triggered if you do

CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -o bug bugpoisMPI.c -lm

./bug

*/

#include "viscosity.h"
#include "fractions.h"

double pr = 1.001;
double pL0, pb0, dt = 0.001;
scalar p[], ps[], rhoc2[], f[];
double Ldomain = 20;
double gamma1 = 1.4, gamma2 = 7.14;
double PI2;
double Ma = 0.01;

face vector alphav[];
(const) face vector alpha = unityf;
(const) scalar rho = unity;

int main()
{  

  /** Size of the domain: */
  X0 = Y0 = -Ldomain/2.;
  size (Ldomain);

  init_grid (256);
 
  TOLERANCE = 1.e-6;
  double PI2 = pr/(1.-pr) + 1./(gamma2*sq(Ma));

  /** Initial pressures */
  pb0 = 1./(pr - 1.);
  pL0 = pr/(pr - 1.);

  alpha = alphav;

  fraction (f, (sq(1.) - sq(x) - sq(y)));
  
  foreach () {
    double pL = pL0 + (pb0 - pL0)/sqrt(sq(x) + sq(y));
    p[] = f[]*pb0 + (1. - f[])*pL;
    ps[] = p[];

    double fc = clamp (f[],0.,1.);
    double invgammaavg = fc/(gamma1 - 1.) + (1. - fc)/(gamma2 - 1.);
    double PIGAMMAavg = (1. - fc)*PI2*gamma2/(gamma2 - 1.);
    rhoc2[] = (p[]*(invgammaavg + 1.) + PIGAMMAavg)/invgammaavg;
  }

  scalar rhs = ps;
  scalar lambda = rhoc2;
  
  foreach() {
      lambda[] = - cm[]/(sq(dt)*rhoc2[]);
      rhs[] = lambda[]*ps[];
  }
  boundary ({lambda, rhs, p});
 
  mgstats mgp = poisson (p, rhs, alpha, lambda, tolerance = TOLERANCE/sq(dt));

  printf("# %i \n", mgp.i);
  foreach ()
    printf("%g %g %g \n", x, y, p[]);

}

/**

~~~gnuplot solution
p "out" u 1:2:3 w p palette
~~~

*/
