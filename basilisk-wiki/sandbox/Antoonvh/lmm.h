/**
## 3rd-order Linear multistep methods 


Linear multistep methods could be intersting. 

*/
#include "run.h"
#include "runge-kutta.h" //Startup procedure

/**

## Adams-Bashfort

Here special care is taken to allow for a variable timestep.
 */

// Adams-Bashfort coeficients for fixed timestep
double coefs[3] = {5./12., -16./12.,  23./12.};

// The tendencies and timesteps of previous steps are stored 
double hs[5];
scalar * dfl[3] = {NULL, NULL, NULL};

void get_coef (double a[3], double dts[3], int i) {
  if (dts[0] == dts[1] && dts[1] == dts[2]) 
    for (int i = 0; i < 3; i++)
      a[i] = coefs[i]*dts[0];
  else {
    double h[3] = {dts[i % 3], (dts[i % 3] + dts[(i - 1) % 3]),
      (dts[0] + dts[1] + dts[2])};
    a[2] = (2*cube(h[0]) - 3*sq(h[0])*h[1] - 3*sq(h[0])*h[2] + 6*h[1]*h[2]*h[0])/
      (6*sq(h[0]) - 6*h[0]*h[1] - 6*h[0]*h[2] + 6*h[1]*h[2]);
    a[1] = (cube(h[0]) -3*sq(h[0])*h[2])/
      (6*h[0]*h[1] - 6*h[0]*h[2] - 6*sq(h[1]) + 6*h[1]*h[2]);
    a[0] = (-cube(h[0]) + 3*sq(h[0])*h[1])/
      (6*h[0]*h[1] - 6*h[0]*h[2] - 6*h[1]*h[2] + 6*sq(h[2]));
  }
}

void LM_AB (scalar * sl, double dt,
	 void (*LU) (scalar * ul, double t, scalar * dul), int i) {
  if (dfl[0] == NULL) {
    for (int i = 0; i < 3; i++)
      dfl[i] = list_clone(sl);
  }
  // Find f(s_n)
  hs[i % 3] = dt;
  LU (sl, t, dfl[i % 3]);
  // Startup procedure
  if (i < 2) 
    runge_kutta (sl, t, dt, LU, 4);
  // LM
  else {
    double coef[3];
    get_coef (coef, hs, i);
    foreach()
      for (int j = 0; j < 3; j++) {
	scalar s, ds;
	for (s, ds in sl, dfl[(i - j)%3])
	  s[] += ds[]*coef[(2 - j)];
      }
  }
}

event cleanup (t = end) {
  for (int i = 0; i < 3; i++) {
    free (dfl[i]);
    dfl[i] = NULL;
  }
}
/**
## Adams Moulton Predictor corrector

The function below makes an Adams bashforth prediction, and a
user-specified number of Adams Moulton corrections. 

method -- Strategy
------ -- ------
1      -- PEC
2      -- PECE
3      -- PECEC
4      -- PECECE
>4     -- etc.
<= 0   -- Converged upto `AM_TOL` 
*/

double AM_TOL = 1e-6; // This does not scale

void LM_AM (scalar *sl, double dt,
	    void (* LU) (scalar *ul, double t, scalar * dul), int i, int method) {
  double coefs_AB[3] = {5./12., -16./12.,  23./12.};
  double coefs_AM[3] = {-1./12., 2./3., 5./12.   };
  if (dfl[0] == NULL) {
    for (int i = 0; i < 3; i++)
      dfl[i] = list_clone(sl);
  }
  // Startup procedure
  if (i == 0)
    LU (sl, t, dfl[i % 3]);
  if (i < 2)  {
    runge_kutta (sl, t, dt, LU, 4);
    LU (sl, t, dfl[i + 1 % 3]);
  }
  // LM-AB predictor
  else {
    scalar * tmp = list_clone (sl);
    foreach() {
      scalar s, ds, tmps;
      double a[list_len(sl)];
      memset (a, 0, list_len(sl)*sizeof(double));
      for (int j = 0; j < 3; j++) {
	int it = 0;
	for (ds in  dfl[(i - j)%3])
	  a[it++] += dt*ds[]*coefs_AB[(2 - j)];
      }
      int it = 0;
      for (s, tmps in sl, tmp)
	tmps[] = s[] + a[it++];      // P
    }
    //LM-AM corrector
    scalar * st2 = list_clone (tmp);
    
    foreach() {
      scalar stemp, stemp2;
      for (stemp, stemp2 in tmp, st2)
	stemp2[] = stemp[];
    }
    double ch = 0;
    int EC_iter = 0;
    do {
      LU (tmp, t, dfl[(i + 1) % 3]); // E
      foreach() {
	scalar s, ds, tmps;
	double a[list_len(sl)];
	memset (a, 0, list_len(sl)*sizeof(double));
	for (int j = 0; j < 3; j++) {
	  scalar s, ds, tmps;
	  int it = 0;
	  for (s, tmps, ds in sl, tmp, dfl[(i + 1 - j)%3])
	    a[it++] += dt*ds[]*coefs_AM[(2 - j)]; 
	}
	int it = 0;
	for (s, tmps in sl, tmp)
	  tmps[] = s[] + a[it++];    // C
      }
      ch = 0;

      if (method <= 0) {
	scalar stemp, stemp2;
	for (stemp, stemp2 in tmp, st2) {
	  double ch1 = change (stemp, stemp2);
	  ch = max (ch1, ch); // ch = max (change ...) bug! 
	}
      }
      EC_iter++;
    } while (method <= 0 ? ch > AM_TOL : 2*EC_iter <= method);
    if (method % 2)
      LU (tmp, t, dfl[(i + 1) % 3]); // E

    foreach() {
      scalar s, st;
      for  (s, st in sl, tmp)
	s[] = st[];
    }
    for (int i = 0; i < 3; i++) {
      free (tmp);
      tmp = NULL;
    }
  }
}

/**
## BDF3 predictor corrector

The BDF3 method stores the solution instead of the tendency, such that
no boundary conditions for the tendency are needed on trees. the
`method' specifies the number of corrections (1: PEC, 2: PECEC, 3:
PECECEC, etc).
 */

void store_solution (scalar * al, scalar * bl) {
  foreach() {
    scalar a, b;
    for (a, b in al, bl)
      b[] = a[];
  }
}

scalar * tend = NULL;
void LM_BDF (scalar *sl, double dt,
	    void (* LU) (scalar *ul, double t, scalar * dul), int i, int method) {
  double coefs_F[4] = {-3./2., 3, -1./2., 3.}; // Explicit coeficients for predictions
  double coefs_B[4] = {18./11., -9./11., 2./11., 6./11.}; // BDF3
  if (dfl[0] == NULL) {
    for (int i = 0; i < 3; i++) //Fix me: only 2 are needed
      dfl[i] = list_clone(sl);
    tend = list_clone (sl);
  }
  store_solution (sl, dfl[i % 3]);
  // Startup procedure
  if (i < 2)  
    runge_kutta (sl, t, dt, LU, 4);
  // Forward predictor
  else {
    if (i == 2)
      LU (sl, t, tend); // Initialize tendency for first prediction
    foreach() {
      scalar s, ds, sp;
      for (s in sl)
	s[] = 0;
      for (int j = 0; j < 3; j++) {
	for (s, sp in sl, dfl[(i - j)%3])
	  s[] += sp[]*coefs_F[j]; 
      }
      for (s, ds in sl, tend)
	s[] += ds[]*dt*coefs_F[3];   // p
    }
    // BDF corrector
    int EC_iter = 0;
    do {
      LU (sl, t, tend);              // E
      foreach() {
	scalar s, ds, sp;
	for (s in sl)
	  s[] = 0;
	for (int j = 0; j < 3; j++) {
	  for (s, sp in sl, dfl[(i - j)%3])
	    s[] += sp[]*coefs_B[j]; 
	}
	for (s, ds in sl, tend)
	  s[] += ds[]*dt*coefs_B[3]; // C
      }
      EC_iter++;
    } while (EC_iter <= method);
  }
}

/**

~~~python
#!/usr/bin/env python

from sympy import *

h1, h2, h3 = symbols('h1 h2 h3')
a = Matrix([[1, -h1, h1**2, -h1**3],
            [0, 1, -2*h1, 3*h1**2],
            [0, 1, -2*h2, 3*h2**2],
            [0, 1, -2*h3, 3*h3**2]])
XtXi = a**-1
init_printing()
pprint (XtXi[0,:], use_unicode = False)
~~~
*/
