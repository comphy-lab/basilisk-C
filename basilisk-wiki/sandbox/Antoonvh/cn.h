/**
# Generic Crank-Nicolson time integration

Here is a implemention for the implicit scheme for time andvancement
of a pde,

$$\frac{\partial f}{\partial t} = \mathcal{L}\left(f\right).$$

The Crank-Nicolson (CA) scheme is 

$$f_{n + 1} = f_n + F \mathrm{d}t,$$

with, 

$$F = \frac{1}{2}\left(\mathcal{L}(f_n) + \mathcal{L}(f_{n+1})\right),$$

Where we have (ab)used the right-hand-side (fist Eq.) operator to also
work with the discrete solution. 

The implementation is _generic_ as it uses an inefficient but flexible
iterative approach to find the implicit update down to a
`CA_TOLERANCE`. The stability of this iterative approach is not likely
to be better than an explicit method.

 */
#include "run.h"
double CA_TOLERANCE = 1e-13;
int maxstps = 200;
/**
We have a function that calls `change()` for a list of scalars and
returns the maximum change.
 */
#include "utils.h"
double max_change (scalar * al, scalar * anl) {
  double max = -1;
  scalar a, an;
  for (a, an in al, anl) {
    double maxi = change (a, an);
    if (maxi > max)
      max = maxi;
  }
  return max;
}

int CA_step (scalar * al, double dt,
	     void (* Lu) ( scalar * al, scalar * dal)) {
  scalar * dal_ex = list_clone (al);
  scalar * dal_im = list_clone (al);
  scalar * aln = list_clone (al);
  scalar * alt = list_clone (al);
    
  scalar a, an, dae, dai;
  // Forward tendency
  Lu (al, dal_ex);
  // Forward Guess
  foreach() {
    for (a, an, dae in al, aln, dal_ex)
      an[] = a[] + dt*dae[];
    for (a, an in aln, alt)
      an[] = a[];
  }
  int stps = 0;
  do {
    stps++;
    if (stps > maxstps) {
      fprintf (stderr, "CA: Max stps reached..\n");
      continue;
    }
    Lu (aln, dal_im);
    foreach() {
      for (a, an, dai, dae in al, aln, dal_ex, dal_im) 
	an[] = a[] + dt*(dai[] + dae[])/2.;
    }
  } while (max_change(aln, alt) > CA_TOLERANCE);

  foreach() {
    for (a, an, in al, aln)
      a[] = an[];
  }
  free (dal_ex);
  free (dal_im);
  free (aln);
  free (alt);
  
  return stps;
}

/**
Here is an implentation for face vector fields:
 */
double changef (face vector a, face vector an) {
  double max = -1;
  foreach_face(reduction(max:max)) {
    if (fabs(a.x[] - an.x[]) > max)
      max = fabs(a.x[] - an.x[]);
    an.x[] = a.x[];
  }
  return max;
}

int CA_step_face (face vector a, double dt,
	     void (* Lu) ( face vector a, face vector da)) {
  face vector da_ex[];
  face vector da_im[];
  face vector an[];
  face vector at[];
  // Forward tendency
  Lu (a, da_ex);
  // Forward Guess
  foreach_face() {
    an.x[] = a.x[] + dt*da_ex.x[];
    at.x[] = an.x[];
  }
  int stps = 0;
  do {
    stps++;
    if (stps > maxstps) {
      fprintf (stderr, "CA: Max stps reached..\n");
      continue;
    }
    Lu (an, da_im);
    foreach_face() {
      an.x[] = a.x[] + dt*(da_im.x[] + da_ex.x[])/2.;
    }
  } while (changef(an, at) > CA_TOLERANCE);
  
  foreach_face()
    a.x[] = an.x[];
  
  return stps;
}
