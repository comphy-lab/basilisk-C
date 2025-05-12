/**
# Propagation of an acoustic disturbance in a tube

This test is the axisymmetric variant of the results presented in
Section 4.1 in [Fuster and Popinet, 2018](#fuster2018). We quantify
the dissipation properties of the all-Mach solver in the acoustic
limit by simulating the propagation of a gaussian pulse of small
amplitude.  The results can be compared with the results obtained with
a classical [Riemann
solver](http://basilisk.fr/sandbox/fuster/RiemannSolverExamples/gaussian.c). */

#include "grid/multigrid.h"

/** 
We use the two-phase flow formulation. */

#include "axi.h"
#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/**
Parameters of the problem. */

double tend = 3.;

int main()
{  
  
  /**
  We make everything dimensionless but this should be improved. */

  size (20. [0]);
  DT = HUGE [0];
  
  X0 = -L0/2.;

  /** 
  The EOS for an adiabatic perfect gas is defined by its polytropic
  coefficient $\Gamma = \gamma = 1.4$. */
  
  gamma1 = 1.4;
  
  /** 
  We perform a convergence study. */
  
  N = 128;
  for (CFLac = 0.01; CFLac <= 100; CFLac *= 5)
    run();
}

event init (i = 0)
{   
  double cson = sqrt(gamma1);
  foreach() {
    f[] = 1.;
    p[] = (1. + 1.e-3*exp(-x*x));
    frho1[] = (1. + (p[] - 1.)/sq(cson));
    q.x[] = 0.;
    q.y[] = 0.;
    fE1[] = p[]/(gamma1 - 1.) + 0.5*sq(q.x[])/frho1[];
  }
}

event endprint (t = tend)
{  
  foreach ()
    if (y < Delta) {
      double xref = fabs(x) - tend*sqrt(gamma1);
      fprintf (stderr, "%g %g %g \n", CFLac, xref, (p[] - 1.)/1.e-3);
    }
  fprintf (stderr, "\n");
}

/**
~~~gnuplot Pressure profile
set xrange[-3:3]
set cblabel 'log10(CFL)'
p "log" u 2:3:(log10($1)) not w lp pt 7 palette, 0.5*exp(-x*x) not w l lw 2 lc 0
~~~ 

## References

~~~bib
@hal{fuster2018, hal-01845218}
~~~
*/
