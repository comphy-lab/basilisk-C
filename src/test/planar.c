/**
# Electrostatic in planar layers

This test reproduces partially the table of convergence of the electric
field for different electrostatic configurations shown in
[Lopez-Herrera et al,
(2011)](/src/references.bib#lopez-herrera2011). The configurations
tested are: (a) both layers are dielectric, (b) both layers are
conducting and (c) one layer conducting and the other dielectric.

We use the implicit electrodynamic solver. The explicit solver
*ehd.h*, which is fully conservative, gives similar results. We also
need the *run()* loop. */

#include "ehd/implicit.h"
#include "run.h"

/**
Variable $ic$ controls the type of electrical configuration tested.
The analytical values used for comparison are stored in E1 and E2.
Finally, *tfinal* controls the instant at which results are
written. Note that for configuration (a) to reach a stationary state
requires just one iteration while for the other configurations about
15 are required. */

int ic;
double tfinal, E1, E2;

phi[top]    = dirichlet(0);
phi[bottom] = dirichlet(1);

#define beta 3.
#define eta 2.

int main() { 
  X0 = Y0 = -0.5;
  DT = 1;
  TOLERANCE = 1e-5;
  int LEVEL;
  for (ic = 0; ic < 3; ic++) {
    switch (ic) {
    case 0: // both layers are dielectric;
      fprintf(stderr, "dielectric-dielectric\n");
      E1 = 2/(1+beta); 
      E2 = 2*beta/(1+beta);
      tfinal = 1;
      break;
    case 1: // both layers are conducting;
      fprintf(stderr, "conducting-conducting\n");
      E1 = 2/(1+eta); 
      E2 = 2*eta/(1+eta);
      tfinal = 15;
      break;
    case 2: //bottom layer conducting upper one dielectric;
      fprintf(stderr, "conducting-dielectric\n");
      E1 = 0.; 
      E2 = 2.;
      tfinal = 15;
      break;
    }
    for (LEVEL = 5; LEVEL <= 7; LEVEL++) {
      N = 1 << LEVEL;
      run(); 
    }
  }
}

event init (t = 0) {
  scalar f[];
  foreach() {
    f[] = (y > 0. ? 0. : 1.0);
    rhoe[] = 0.;
  }

  foreach_face(){
    double T = (f[]+f[-1,0])/2.;

    /**
     If in configuration (a) we use an harmonic interpolation for the
     permittivity the exact values of the electric field are
     recovered. */

    epsilon.x[] = 1./(T/beta + (1. - T)); //Exact for the diel-diel
    // epsilon.x[] = T*beta + (1. - T); // Non exact 
    K.x[] = (ic == 0 ? 0. : (ic == 1 ? T*eta+(1. - T) : eta*T));
  }
}

event set_dt (i++)
   dt = dtnext (DT);

/**
In the last step the relative error of the electric field, $E=-\nabla \phi$, 
in each medium is reported. */

event result (t = tfinal) {
  scalar E[];
  foreach()
    E[] = (phi[0,-1] - phi[0,1])/(2.*Delta);
  stats s = statsf (E);
  fprintf (stderr, "%d %g %g\n", N, 100*fabs(1.- s.max/E2),
	   E1 > 0. ? 100*fabs(1.- s.min/E1) : 0. );
}
