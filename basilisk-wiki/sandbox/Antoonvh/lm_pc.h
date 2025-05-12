/**
# A 3rd-order accurate linear multistep predictor--corrector method

With equal timesteps (`h`) the scheme consists of a predictor
(2nd-order accurate Adams-Bashforth):

$$Y_{n + 1} = y_n + \frac{h}{2}\left(3f(y_n) - f(y_{n - 1})\right), $$

and a corrector with positive coefficients only, 

$$y_{n + 1} = \frac{1}{5}\left(4y_n + y_{n - 1}\right) + \frac{h}{5}\left(2f(Y_{n+1}) + 4f(y_{n})\right). $$

Each step has two evaluations (PECE) and we must store the tendencies the the solution at the previous step.  
 */
#include "run.h"
#include "runge-kutta.h"

scalar * dfl[3] = {NULL, NULL, NULL};

void store_solution (scalar * al, scalar * bl) {
  foreach() {
    scalar a, b;
    for (a, b in al, bl)
      b[] = a[];
  }
}

void LM_PC (scalar *sl, double dt,
	    void (* LU) (scalar *ul, double t, scalar * dul),
	    int i) {
  double coefs_B[2] = {2./5., 4./5.}; // implicit coeficients for corrections
  if (dfl[0] == NULL) 
    for (int i = 0; i < 3; i++) 
      dfl[i] = list_clone(sl); // fixme boundary conditions for tendencies
  
  scalar * stmp = list_clone (sl);
  store_solution (sl, stmp);
  
  // Startup procedure
  if (i < 2)  {
    runge_kutta (sl, t, dt, LU, 4);
    LU (sl, t, dfl[i % 2]);
  }
  else {
    // Predictor
    foreach() {
      scalar s, ds, sp;
      for (s, ds in sl, dfl[(i - 1) % 2])
	s[] += 3./2.*ds[]*dt;
      for (s, ds in sl, dfl[(i) % 2])
	s[] -= 1./2.*ds[]*dt;         // P
    }
    // corrector
    LU (sl, t, dfl[i % 2]);           // E
    foreach() {
      scalar s, ds, sp, st;
      for (s, st, sp in sl, stmp, dfl[2])
	s[] = (4.*st[] + sp[])/5.;
      for (int j = 0; j < 2; j++) {
	for (s, ds in sl, dfl[(i + j)%2])
	  s[] += dt*ds[]*coefs_B[j]; // P 
      }
    }
    LU (sl, t, dfl[i % 2]);          // E
  }
  store_solution (stmp, dfl[2]);
  delete (stmp);
  free (stmp);
}

event cleanup (t = end) {
  if (dfl[0]) 
    for (int i = 0; i < 3; i++) {
      free (dfl[i]);
      dfl[i] = NULL;
    }
}
