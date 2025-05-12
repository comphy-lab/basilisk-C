/**
#  Periodic couette flow  with imposed shear stress at the lower wall

 
This is a simple test of solution of Navier Stokes with mixed BC at the wall. No slip at upper wall were $u=U_0$, and we impose the shear stress at the lower wall
$$\mu \frac{\partial u}{\partial y}|_0=\tau_w$$
where $\tau_w$ is given. In a classical plane Couette flow :$u=U_0$ at the upper wall and 0 at lower.

Analytical steady solution is just $\tau(y)= \tau_w$ constant in the layer and
$$ u(y) = U_0 - \frac{\tau_w h}{\mu}(1-\frac{y}{h})$$


*/

#include "navier-stokes/centered.h"
#define LEVEL 4 
double tauw;
scalar mu_eq[];
face vector muv[];
 
/**
The domain is one unit long. $0<x<1$, $0<y<1$
*/

int main() {
  L0 = 1.0;
  origin (0., 0);
  tauw=.50;
/**
 Boundary conditions are periodic
*/
    periodic (right);
/**
 classical couette OK
  no slip at the top wall which moves at unit velocity
  no slip at the bottom
*/  
#if 0
  u.t[top] = dirichlet(1);
  u.n[top] = dirichlet(0); 
  u.n[bottom] = dirichlet(0); 
  u.t[bottom] =  dirichlet(0);
#endif
 
/**
 impose slope at bottom: OK
*/
#if 0
  u.t[top] = neumann(0);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0); 
  u.t[bottom] = neumann(-tauw);
#endif  
/**
 impose shear at the wall, note the minus sign (external normal),  and the face value of $\mu$ 
*/
#if 1
  u.t[top] =  dirichlet(1);
  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0); 
  u.t[bottom] = neumann(mu.y[] ? -tauw/mu.y[] : 0.);
#endif  
/**
we tune the NS solver, as $U=u(y)e_x$. `stokes` is defined as `true`, so there is  no CFL condition
*/
 // DT = 0.1;
  NITERMAX = 100;
  TOLERANCE = 1e-3;
  stokes = true; 
  run(); 
}
/**

*/

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
  pressure gradient,  acceleration `mdpdx` if we want to do a Poiseuille flow
 $$-\frac{\partial p}{\partial x} = 0 $$
 $$-\frac{\partial p}{\partial y} = 0 $$
*/
  const face vector mdpdx[] = {0,0};
  a = mdpdx;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = y;
    u.y[] = 0;
  }
}
/**
 old value of the velocity is saved
*/
scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}
/**
 so that when it does not more change we are converged
*/
event conv (t += 1; i <= 10000) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g %g %g \n",t,interpolate (u.x, L0/2, .999),interpolate (u.x, L0/2, .999),du);
    if (i > 0 && du < 5.0e-5)
        return 1; /* stop */
}
/**
## Implementation of the Bagnold viscosity
*/
event properties (i++) {
    /** Compute viscosity, here constant and equal to one! 
     */
    foreach() {  
     mu_eq[] =1;
    }
    boundary ({mu_eq});
    foreach_face() {
        muv.x[] = (mu_eq[] + mu_eq[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}
/**
  Save profiles
*/
event profiles (t += 1  )
{
  FILE * fp = fopen("xprof", "w");
    scalar shearS[];
    foreach()
    shearS[] = mu_eq[]*(u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shearS});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g  %g  %g \n", y, interpolate (u.x, L0/2, y),interpolate (shearS, L0/2, y));
    fclose (fp);
}

event profile (t = end)
{
  scalar shearS[];
  foreach()
    shearS[] = mu.y[]*(u.x[0,0] - u.x[0,-1])/(Delta);
  boundary ({shearS});
  for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
    fprintf (stderr, "%g  %g  %g \n", y, interpolate (u.x, L0/2, y),interpolate (shearS, L0/2, y));
}

/**
## Run
 
To run the program
 
~~~bash
  qcc -g -O3 -o couette_muw  couette_muw.c -lm
 ./couette_muw
 
~~~

## Results and plots

we compare 
here with the steady solution
$$ u(y) = U_0 + \frac{\tau_w h}{\mu}(\frac{y}{h}-1)$$
where the stress is constant
$$\tau(y)= \tau_w$$
 
 
 Plots of the velocity and $\tau$ 
 with $\tau_w=0.5$, $U_0=1$, $h=1$, $\mu=1$ 


~~~gnuplot Velocity,   $\tau= \mu du/dy$ profiles computed
 set xlabel "y"
 set xlabel "U, tau"
 set key left
 p'xprof' u 1:2 t'U comp.' , 1-.5*(1-x) t'U anal.',''u 1:($3) t'tau comp.' , .5 t'imposed tau'
~~~

 
*/
