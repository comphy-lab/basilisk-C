/**
# Zero reflection of a wave propagating across an interface between two fluids with impedance matching

In this test proposed by [Denner et al, 2018](#denner2018) a linear wave
propagating in an ideal gas is completely transmitted to another ideal
gas with the same acoustic impedance. */

#include "grid/multigrid1D.h"
#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/** 
Parameters of the problem. */

double tend = 0.6;
double cflac = 0.1;
double uper = 1e-4;
double freq = 15.;
double p0, rho20;

/**
Fixme: `cflac` should be a parameter of two-phase-compressible.h. */

event stability (i++)
{
  double Delta_min = HUGE;
  foreach (reduction(min:Delta_min))
    if (Delta < Delta_min)
      Delta_min = Delta;
  dtmax = Delta_min*cflac;
  DT = dtmax;
}

int main()
{

  /** 
  The EOS for an adiabatic perfect gas is defined by its polytropic
  coefficient $\Gamma = \gamma$. */
  
  gamma1 = 9.872;
  gamma2 = 2.468;
  rho20  = 1./0.25;

  p0 = 1./gamma1;
  
  N = 512;
  run();
}

event init (i = 0)
{
  foreach() {
    double perturb = uper*exp(- sq ((x - 0.3)*freq));
    f[] = (x < 0.5);
    p[] = p0 + perturb;
    frho1[] = f[]*(1. + perturb);
    frho2[] = (1. - f[])*rho20;
    q.x[] = (frho1[] + frho2[])*perturb;
    fE1[] = f[]*p[]/(gamma1 - 1.) + 0.5*sq(q.x[]/(frho1[] + frho2[]))*frho1[];
    fE2[] = (1. - f[])*p[]/(gamma2 - 1.) + 0.5*sq(q.x[]/(frho1[] + frho2[]))*frho2[];
  }
}

event endprint (t = tend)
{  
  scalar perr[];
  foreach () {
    perr[] = fabs((p[] - p0)/uper - exp(- sq((x - 0.6)*freq*sqrt(gamma1*rho20/gamma2))));
    fprintf (stderr, "%g %g %g \n", t, x, p[] - p0);
  }
  fprintf (stderr, "error %g\n", statsf(perr).sum);
}

/**
~~~gnuplot Perfect transmission
set ylabel '{/Symbol D}p/{/Symbol D}p_0'
set xlabel 'x'
set samples 1000
set arrow from 0.5,0 to 0.5,1 nohead lc 0
set label 'fluid 1' at 0.4,0.8
set label 'fluid 2' at 0.51,0.8
p "log" u 2:($3/0.0001) t 'tend' w l lc 1,				        \
   exp(-((x - 0.3)*15)**2) t 't = 0' w l,					\
   exp(-((x - 0.6)*15*sqrt(9.872/2.468/0.25))**2) t 'tend (theory)' w l
~~~

## References

~~~bib
@article{denner2018,
title = {Pressure-based algorithm for compressible interfacial flows
with acoustically-conservative interface discretisation},
journal = {Journal of Computational Physics},
volume = {367},
pages = {192-234},
year = {2018},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2018.04.028},
url = {https://www.sciencedirect.com/science/article/pii/S0021999118302535},
author = {Fabian Denner and Cheng-Nian Xiao and Berend G.M. {van Wachem}}
}
~~~
*/
