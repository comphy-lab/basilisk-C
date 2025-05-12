/**
# The Noble-Abel Stiffened-Gas (NASG) Equation Of State

This EOS is typically used in combination with the [two-phase
compressible solver with thermal effects](/src/compressible/thermal.h).

The general form of the NASG EOS ([Le Métayer & Saurel, 2016](#metayer2016)) is
$$
\rho_i e_i = \frac{p_i + \Gamma_i \Pi_i}{\Gamma_i - 1}(1 - \rho_i b_i) + \rho_i q_i
$$
with $\rho_i$, $e_i$ and $p_i$ the densities, internal energies and
pressures of each phase.

These are the coefficients of the NASG EOS for each phase. */

double gamma1 = 1.4 [0], gamma2 = 1.4 [0], PI1 = 0., PI2 = 0.;
double b1 = 0., b2 = 0.;
double q1 = 0., q2 = 0.;
double cv1 = 0., cv2 = 0.;

/**
## Sound speed

In mixture cells, this function returns the maximum between the speeds
in both phases. */

double sound_speed (Point point)
{
  double fc = clamp (f[],0.,1.);
  double c2speed1 = 0., c2speed2 = 0.;

  double Ek = 0.;
  foreach_dimension()
    Ek += sq(q.x[]);
  Ek /= 2.*(frho1[] + frho2[]);
  
  if (fc > 0.00001) {
    double fe1 = fE1[] - fc*Ek;
    double p  = fe1/fc*(gamma1 - 1.) - gamma1*PI1;
    c2speed1 = fc*gamma1*(p + PI1)/frho1[];
  }
  
  if (fc < 0.99999) {
    double fe2 = fE2[] - (1. - fc)*Ek;
    double p  = fe2/(1. - fc)*(gamma2 - 1.) - gamma2*PI2;
    c2speed2 = (1. - fc)*gamma2*(p + PI2)/frho2[];
  }

  return sqrt (max (c2speed1, c2speed2));
}

/**
## Average pressure */

#define PIGAMMA	 double invgammaavg = (fc - frho1[]*b1)/(gamma1 - 1.) +     \
    (1. - fc - frho2[]*b2)/(gamma2 - 1.),				    \
    PIGAMMAavg = (fc - frho1[]*b1)*PI1*gamma1/(gamma1 - 1.) + frho1[]*q1 + \
    (1. - fc - frho2[]*b2)*PI2*gamma2/(gamma2 - 1.) + frho2[]*q2

double average_pressure (Point point)
{
  double fc = clamp (f[],0.,1.);
  PIGAMMA;
  double Ek = 0.;
  foreach_dimension()
    Ek += sq(q.x[]);
  Ek /= 2.*(frho1[] + frho2[]);
  return (fE1[] + fE2[] - Ek - PIGAMMAavg)/invgammaavg;
}

/**
## Bulk compressibility of the mixture

i.e. $\rho c^2$. */

double bulk_compressibility (Point point)
{
  double fc = clamp (f[],0.,1.);
  // Arithmetic mean of the mixture compressibility
  double rhoc2v1 = fc ? gamma1*(p[] + PI1)/(1. - frho1[]*b1/fc) : 1.;
  double rhoc2v2 = (1. - fc) ? gamma2*(p[] + PI2)/(1. - frho2[]*b2/(1. - fc)) : 1.;
  return fc*rhoc2v1 + (1. - fc)*rhoc2v2;
}

/**
## Internal energy */

double internal_energy (Point point, double fc)
{
  PIGAMMA;
  return p[]*invgammaavg + PIGAMMAavg;
}

/**
## Average temperature */

double average_temperature (Point point, double p)
{
  double fc = clamp (f[],0.,1.);
  double rhocpmcvavg = (cp1 - cv1)*frho1[] + (cp2 - cv2)*frho2[];
  double const1 = (fc - frho1[]*b1) + (1. - fc - frho2[]*b2);
  double const2 = (fc - frho1[]*b1)*PI1 + (1. - fc - frho2[]*b2)*PI2;  
  return (const1*p + const2)/rhocpmcvavg;
}

/**
## Thermal expansion coefficient */

double thermal_expansion (Point point)
{
  double fc = clamp (f[],0.,1.);
  return (1. - fc)*(gamma2 - 1.)*cv2/((gamma2 - 1.)*cv2*Ts[] + b2*(ps[] + PI2));
}

/**
# See also

* [The Mie-Gruneisen Equation Of State](Mie-Gruneisen.h)

# References

~~~bib
@article{metayer2016,
  title={The {N}oble--{A}bel stiffened-gas equation of state},
  author={Le M{\'e}tayer, Olivier and Saurel, Richard},
  journal={Physics of Fluids},
  volume={28},
  number={4},
  year={2016},
  publisher={AIP Publishing}
}
~~~
*/
