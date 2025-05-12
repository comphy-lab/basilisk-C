/**
# Gauss-Legendre Runge-Kutta schemes

Implicit time-integration with 2nd, 4th, 6th or 8th-order accurate
time advancing.

The scheme (in principle) takes any Butcher array, the coefficients
for the aforementioned Gauss-Legendre scheme are given.
*/
#include "run.h"
#ifndef RKORDER
#define RKORDER (4)
#endif
/**
see,  
*Implicit Runge-Kutta Processes* by J.C. Butcher
[(pdf.)](https://www.ams.org/journals/mcom/1964-18-085/S0025-5718-1964-0159424-9/S0025-5718-1964-0159424-9.pdf)

*/
#if (RKORDER == 2)
#define STAGES (1)
double an[STAGES][STAGES] = {{1./2.}};
double bn[STAGES] = {1.};

#elif (RKORDER == 4)
#define STAGES (2)
#define SQRT3 (1.73205080756887729352744634150587236694280525)
double an[STAGES][STAGES] = {{1./4.           , 1./4. - SQRT3/6.},
			     {1./4. + SQRT3/6., 1./4.}};
double bn[STAGES] = {1./2., 1./2.};

#elif (RKORDER == 6)
#define STAGES (3)
#define SQRT15 (3.8729833462074168851792653997823996108329217)
double an[STAGES][STAGES] =
  {{5./36.             , 2./9. - SQRT15/15., 5./36. - SQRT15/30. },
   {5./36. + SQRT15/24., 2./9.             , 5./36. - SQRT15/24. },
   {5./36. + SQRT15/30., 2./9. + SQRT15/15., 5./36.}};
double bn[STAGES] = {5./18., 4./9., 5./18.};

#elif (RKORDER == 8)
#define STAGES (4)
#define SQRT30 (5.47722557505166113456969782800802133952744694)
#define SQRTEP (0.86113631159405257522394648889280950509572537)
#define SQRTEM (0.33998104358485626480266575910324468720057587)
#define W1  (1/8. - SQRT30/144.)
#define W1p (1/8. + SQRT30/144.)
#define W2  (SQRTEP/2.)
#define W2p (SQRTEM/2.)
#define W3  (W2* (1./6. + SQRT30/24.))
#define W3p (W2p*(1./6. - SQRT30/24.))
#define W4  (W2* (1./21. + 5.*SQRT30/168.))
#define W4p (W2p*(1./21. - 5.*SQRT30/168.))
#define W5  (W2 - 2*W3)
#define W5p (W2p- 2*W3p)
double an[STAGES][STAGES] =
  {{W1           , W1p - W3 + W4p, W1p - W3 - W4p, W1 - W5},
   {W1 - W3p + W4, W1p           , W1p - W5p     , W1 - W3p - W4},
   {W1 + W3p + W4, W1p + W5p     , W1p           , W1 + W3p - W4},
   {W1 + W5      , W1p + W3 + W4p, W1p + W3 - W4p, W1}};
double bn[STAGES] = {2*W1, 2*W1p, 2*W1p, 2*W1};
#endif

scalar * kl[STAGES] = {NULL};
// for list_len(fl) = n, -> kl[stages][n]

void A_Time_Step (scalar * fl, double dt,
		  void (* Lu) (scalar * ul, scalar * dul),
		  double Tol) {
  // Allocate the global `kl` scalars once
  scalar * kj;
  if (kl[0] == NULL) 
    for (int ii = 0; ii < STAGES; ii++) 
      for (int g = 0; g < list_len (fl); g++) {
	scalar c = new_scalar ("ki");
	kl[ii] = list_append (kl[ii], c);
	foreach() {
	  for (scalar k in kl[ii]) 
	    k[] = 0;
	}
      }
  // reuse values of k
  ;
  // Allocate and initialize fields for comparisons
  scalar * templ[STAGES];
   for (int ii = 0; ii < STAGES; ii++) {
    templ[ii] = list_clone (kl[ii]);
    foreach()
      for (scalar s in templ[ii])
	s[] = -9999;
  }
  // Solution field at stages 
  scalar * ftl = list_clone (fl);
  // Set iterative details
  double em = 0;
  int it = 0, itmin = STAGES + 1, itmax = 20;
  // Start iteration
  do {
    it++;
    em = 0;
    // Foreach stage
    for (int ii = 0; ii < STAGES; ii++) {
      foreach() {
	// Compute solution estimate at stage
	scalar f, ft;
	for (f, ft in fl, ftl) 
	  ft[] = f[];
	for (int j = 0; j < STAGES; j++) {
	  scalar ko, ft;
	  scalar * kj = kl[j];
	  for (ft, ko in ftl, kj)
	    ft[] += dt*ko[]*an[ii][j];
	}
      }
      // Compute corresponding `k`
      Lu (ftl, kl[ii]);
      // Check for convergence 
      scalar k, temp;
      for (k, temp in kl[ii], templ[ii]) {
	double et = change (k, temp);
	if (et > em)
	  em = et;
      }
    }
    // Prey for convergence...
  } while ((em > Tol || it < itmin) && it < itmax);
  if (it == itmax)
    fprintf (stderr, "IRK convergence not reached. "
	     "i: %d t: %g: change: %g\n", iter, t, em);
  // Update solution
  scalar * ki = NULL;
  foreach() {
    for (int ii = 0; ii < STAGES; ii++) {
      scalar f, k;
      scalar * ki = kl[ii];
      ki;
      for (f, k in fl, ki)
	f[] += dt*bn[ii]*k[];
    }
  }
  // Clean up temporary fields
    for (int ii = 0; ii < STAGES; ii++) {
      delete (templ[ii]); free(templ[ii]); templ[ii] = NULL;
    }
    delete (ftl); free(ftl); ftl = NULL;
}

event rm_dfl (t = end) {
  // More cleaning
  for (int ii = 0; ii < STAGES; ii++) {
    delete (kl[ii]); free (kl[ii]); kl[ii] = NULL;
  }
}

/**
## Test

* [6th-order accurate scheme](tGLrk.c)

## Usage

* [8th-order accurate scheme for the Korteweg-de Vries equation solver](KdV.h)
*/
