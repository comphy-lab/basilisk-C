/**
# The Mie--Gruneisen Equation of State

This EOS is typically used in combination with the
[two-phase compressible solver](/src/compressible/two-phase.h).

The general form of the
[Mie--Gruneisen](https://en.wikipedia.org/wiki/Mie%E2%80%93Gr%C3%BCneisen_equation_of_state)
EOS can be written
$$
\rho_i e_i = \frac{p_i + \Gamma_i \Pi_i}{\Gamma_i - 1}
$$
with $\rho_i$, $e_i$ and $p_i$ the densities, internal energies and
pressures of each phase.

These are the coefficients of the Mie-Gruneisen EOS for each phase. */

double gamma1 = 1.4 [0], gamma2 = 1.4 [0], PI1 = 0., PI2 = 0.;

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

#define PIGAMMA	 double invgammaavg = fc/(gamma1 - 1.) + (1. - fc)/(gamma2 - 1.), \
    PIGAMMAavg = fc*PI1*gamma1/(gamma1 - 1.) + (1. - fc)*PI2*gamma2/(gamma2 - 1.)

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
  PIGAMMA;
  return (p[]*(invgammaavg + 1.) + PIGAMMAavg)/invgammaavg;
}

/**
## Internal energy */

double internal_energy (Point point, double fc)
{
  PIGAMMA;
  return p[]*invgammaavg + PIGAMMAavg;
}
