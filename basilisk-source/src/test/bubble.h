/**
# A generic compressible gas bubble in a liquid */

#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"
#include "tension.h"
#include "compressible/tension.h"
#include "rayleigh-plesset.h"

/**
Static refinement is used with a level of refinement varying from
MINLEVEL far from the bubble to MAXLEVEL close to the bubble. */

int MINLEVEL = 4, MAXLEVEL;

/**
The liquid and gas densities as well as the initial gas pressure `pg0`
and the pressure "at infinity" `pinf`. */

double rhoL = 1., rhoG = 0.001;
double pg0 = 100., pinf = 100.;

/**
The final time and the bubble radius.*/

double tend = 1.;
double R0 = 1.;

/**
Boundary conditions: imposed pressure and inflow "at infinity". */

p[right]   = dirichlet (pinf);
q.n[right] = neumann (0.);

#if dimension > 1
p[top]     = dirichlet (pinf);
q.n[top]   = neumann (0.);
#endif

event init (i = 0)
{

  /**
  We compute the reference Rayleigh-Plesset and Keller-Miksis solutions. */
  
  static bool first = true;
  if (first) {
    struct RPdata RPd = {
      .rhol  = rhoL,
      .pliq = pinf,
      .p0 = pg0,
      .sigma = f.sigma,
      .gamma = gamma2,
      .R0 = R0,
      .visc = mu1,
      .cson = sqrt(gamma1*(pinf + PI1)/rhoL)
    };
    
    FILE * fp = fopen("RP.dat", "w");
    Integrate_RP (fp, 0., tend, &RPd);
    fclose (fp);
  
    fp = fopen("RPinviscid.dat", "w");
    RPd.visc = 0;
    Integrate_RP (fp, 0., tend, &RPd);
    fclose (fp);

    first = false;
  }

  /**
  The static mesh refinement. */

#if TREE
  for (int l = MINLEVEL ; l <= MAXLEVEL; l++)
    refine (level < l && sqrt(sq(x) + sq(y)) < (2.5*R0 + 4.*sqrt(2.)*L0/(1 << (l - 1))));
#endif
  
  /**
  The initial volume fraction, densities, energies and pressure. */

  fraction (f, sq(x) + sq(y) - sq(R0));
  
  foreach() {
    frho1[] = f[]*rhoL;
    frho2[] = (1. - f[])*rhoG;

    /** 
    The initial liquid pressure field is approximated from the
    solution in the incompressible limit. */
    
    double r = sqrt(sq(x) + sq(y));
    double pL = pinf*(1. - R0/r) + (pg0 - 2*f.sigma/R0)*R0/r;
	
    fE1[] = f[]*(pL + PI1*gamma1)/(gamma1 - 1.);
    fE2[] = (1. - f[])*pg0/(gamma2 - 1.);

    p[] = average_pressure (point);
  }
}

/**
Some statistics. */

event logfile (i++)
{
  scalar pgas[], keliq[];
  double volume = 0.;
  
  foreach (reduction(+:volume)) {
   
    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);
    Ek /= 2.*(frho1[] + frho2[]);
    
    /**
    The gas pressure is recovered from the energy. */
    
    pgas[] = average_pressure (point)*(1. - f[]);
    keliq[] = Ek*f[];
    volume += dv()*(1. - f[]);
  }

  if (i == 0)
    fprintf (stderr, "t volume statsf(keliq).sum statsf(pgas).sum/volume\n");
  fprintf (stderr, "%10.6f %10.6f %10.8f %10.6f\n",
	   t, volume, statsf(keliq).sum, statsf(pgas).sum/volume);
}

/**
The bubble shape as a function of time. */

event profiles (t += tend/10; t <= tend) {
  output_facets (f);
}
