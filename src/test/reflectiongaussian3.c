/**
# Transmission/reflection of a wave propagating across an interface between two fluids 

In this test proposed by [Denner et al. 2018](#denner2018) a linear
wave propagating in an ideal gas is partially transmitted to another
ideal gas with a different acoustic impedance. */

#include "grid/multigrid1D.h"
#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/** 
Parameters of the problem. */

double tend = 0.1;
double uper = 0.0001;
double freq = 4000.;
double uref = 347.8;
double p0, rho20, rho10;

int main()
{  
  
  CFLac = 0.25;
  /** 
  The EOS for an adiabatic perfect gas is defined by its polytropic
  coefficient $\Gamma = \gamma = 1.4$. */

  gamma1 = 1.4;
  gamma2 = 1.667;
  rho20  = 0.164/1.157;
  rho10 = 1;;

  p0 = 1./gamma1;
  freq *= sqrt(gamma2/gamma1/rho20)/uref;
  
  /**
  We perform a convergence study. */

  for (N = 256; N <= 1024; N *= 2)
    run();
}

event init (i = 0)
{
  foreach() {
    double perturb = uper*exp(- sq((x - 0.4)*freq));
    f[] = (x > 0.5);
    p[] = p0 + perturb;
    frho1[] = f[]*rho10*(1. + perturb);
    frho2[] = (1. - f[])*rho20;
    q.x[] = (frho1[] + frho2[])*perturb*sqrt(gamma2*p[]/rho20);
    fE1[] = f[]*p[]/(gamma1 - 1.) + 0.5*sq(q.x[]/(frho1[] + frho2[]))*frho1[];
    fE2[] = (1. - f[])*p[]/(gamma2 - 1.) + 0.5*sq(q.x[]/(frho1[] + frho2[]))*frho2[];
  }
}

event endprint (t = tend) 
{
  foreach()
    fprintf (stderr, "%i %f %f\n", N, x, (p[] - p0)/1e-4);
}

/**
~~~gnuplot Reflected wave
ZR = 1.
ZL = 0.164/1.157*sqrt(1.667*1.157/0.164/1.4)
set ylabel '{/Symbol D}p/{/Symbol D}p_0'
set xlabel 'x'
set cblabel '{/Symbol s}/{/Symbol D}x'
set xrange[0.2:0.4]
p "log" u 2:3:(0.1*$1) not w l palette, (ZR-ZL)/(ZR+ZL) t 'theory' w l lc 0
~~~ 

~~~gnuplot Transmitted wave
set xrange[0.5:0.65]
p "log" u 2:3:(0.1*$1) not w l palette, 1./(ZR+ZL)*2*ZR t 'theory' w l lc 0
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
