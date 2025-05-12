/** Here we try to obtain a flow driven by a pressure drop on a tri-periodic cubic domain. 

Solution: since the pressure drop cannot be imposed in a periodic domain, the solution is to replace it with an equivalent body force i.e. we have in the r.h.s. of the momentum equation
$$
-\partial_x p + \rho a_x
$$
so that adding a constant pressure gradient is equivalent to adding an acceleration
$$
a_x = - 1/\rho\partial_x p
$$
*/

#define LEVEL 6
#include "grid/octree.h"
#include "navier-stokes/centered.h"

#define rhoval 1000.
#define muval 60.

#define mydt 0.00375
#define maxiter (20)

double deltau;
scalar un[];

int main() {		
  
  L0 = 1; 
  
  DT = mydt;
  
  init_grid(1 << (LEVEL));


  /* Tri-periodic flow driven by body force */ 
  
  foreach_dimension()
    periodic(top);

  /*
  p[left]  = dirichlet(0);
  p[right] = dirichlet(-10); */
  const face vector dp[] = {10./L0/rhoval, 0, 0};
  a = dp;
  
  run();
}


event init (i = 0) {
  origin (0, 0, 0);
  
/* set dynamic viscosity */
  const face vector muc[] = {muval, muval, muval};
  mu = muc;

/* set density of the flow */ 
  const scalar rhoc[] = rhoval;
  rho = rhoc;
  
/* The flow is at rest initially. */
  foreach(){
    u.x[] = 0.;
    un[] = u.x[];
  }
}

/* Logging the variable "deltau" that measures the difference between the current velocity to the initial one. */
event logfile (i++;i<maxiter){
  deltau = change (u.x, un);
  fprintf (ferr, "%g %g\n", t, deltau);
}

/** If the line "periodic(right)" is commented, the flow will be set in motion and deltau (defined above) will be different than 0, otherwise the flow won't move and deltau will remain null 

~~~gnuplot deltau and time
plot 'log' u 1:2 w l
~~~
*/

