/**
# Advection of two fluids at different pressures

This problem, initially proposed by [Johnsen and Colonius,
2006](#johnsen2006) tests the solver for the advection of an interface
between two different ideal gases at uniform velocity and pressure
quantifying the amplitude of spurious pressure and velocity
oscillations induced by the method when advecting an interface with
different material properties. See also Figure 4 in [Fuster and
Popinet, 2018](#fuster2018).

$$(\rho,u,p,\gamma)^T_L = (1, 0.5, 1/1.4, 1.2)^T$$
$$(\rho,u,p,\gamma)^T_R = (10, 0.5, 1/1.4, 1.4)^T$$

By construction the method should keep the pressure and velocity
uniform. During the advection step the energy is advected avoiding any
diffusion and therefore both the provisional pressure and velocity are
uniform and do not need to be corrected during the projection step. */

#include "grid/multigrid1D.h"
#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"

/**
We set the problem parameters, size domain and boundary conditions. */

double rhoL = 1., rhoR = 10.;
double pL, pR;
double tend = 8;

int main()
{
  L0 = 2. [1];
  X0 = - L0/2.;
  
  periodic (right);
  
  pL = pR = 1./1.4;
  gamma1 = 1.2;
  gamma2 = 1.4;
 
  N = 128;
  run();
}

/** 
The initial conditions are: */

event init (i = 0)
{  
  double u0 = 0.5;
  foreach() {
    f[] = (x < 0.);
    p[] = f[]*pL + (1. - f[])*pR;
    frho1[] = f[]*rhoL;
    frho2[] = (1. - f[])*rhoR;
    fE1[] = f[]*pL/(gamma1 - 1.) +  0.5*frho1[]*sq(u0);
    fE2[] = (1. - f[])*pR/(gamma2 - 1.) +  0.5*frho2[]*sq(u0);
    q.x[] = (frho1[] + frho2[])*u0;
  }
}

/** 
We output the field variables at the end of the simulation. */

event outputdata (t = tend)
{  
  scalar perr[], uerr[];
  
  foreach () {
    perr[] = fabs(p[] - 1./1.4);
    uerr[] = fabs(q.x[]/rho[] - 0.5);

    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);

    fprintf (stderr, "%g %g %g %g %g %g %g \n", x, t, p[], rho[], q.x[]/rho[], f[],
	     (fE1[] + fE2[] - 0.5*Ek/rho[])/(f[]/(gamma1 - 1.) + (1. - f[])/(gamma2 - 1.)));
  }

  stats sp = statsf(perr), su = statsf(uerr);
  fprintf (stdout, "%g %g\n", sp.sum/sp.volume, su.sum/su.volume);
  fflush (stdout);
  assert (sp.sum/sp.volume < 2e-9 && su.sum/su.volume < 2e-9);
}

/**
The results below are those displayed in Figure 4 of [Fuster and
Popinet, 2018](#fuster2018).

~~~gnuplot Pressure profile
set xlabel 'x'
set ylabel 'p'
unset key
p "log" u 1:3 w lp
~~~ 
 
~~~gnuplot Density profile
set xlabel 'x'
set ylabel '{/Symbol r}'
p "log" u 1:4 w lp
~~~ 

~~~gnuplot Velocity profile
set xlabel 'x'
set ylabel 'u'
p "log" u 1:5 w lp
~~~ 

~~~gnuplot VOF profile
set xlabel 'x'
set ylabel 'VOF function'
p "log" u 1:6 w lp
~~~

### References

~~~bib
@article{johnsen2006,
title = {Implementation of {WENO} schemes in compressible multicomponent flow problems},
journal = {Journal of Computational Physics},
volume = {219},
number = {2},
pages = {715-732},
year = {2006},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2006.04.018},
url = {https://www.sciencedirect.com/science/article/pii/S0021999106002014},
author = {Eric Johnsen and Tim Colonius}
}

@hal{fuster2018, hal-01845218}
~~~
*/
