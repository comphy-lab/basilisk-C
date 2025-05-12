/**
   
# Test the boussinesq equations

This is a very simple implementation of the Boussinesq equations.  (obviously
too simple because we observe that a flow that should be at rest is not at
rest, see below )

*/

#include "navier-stokes/centered.h"
#include "output.h"

/**
   We add a buoyancy tracer `b` which is advected by the flow. In the momentum
   equations, we have
   
   du/dt = -grad P + b
   
*/

#include "tracer.h"
#include "diffusion.h"

face vector av[];
scalar b[];
scalar * tracers = {b};

event init (t = 0) {
  a = av;
}

/**
   Here is the addition of the buoyancy component to the momentum equations
*/

event acceleration (i++) {

  coord gravity_dir = {0, 1};

  foreach_face()
    av.x[] = gravity_dir.x*face_value(b, 0);
}


/**
   Main driver: define the grid and geometry (periodic left-right)
*/

int main(int argc,char* argv[]) {

  N = 32;
  periodic (right);

  run();

}


/**
## Initial conditions

Initial conditions correspond zero velocity and uniform buoyancy.
Set no flux BC at top and bottom boundaries, homogeneous neumann for buoyancy
*/

event init (t=0) {
  foreach() {
    b[] = 1.;
    u.x[] = 0.;
    u.y[] = 0.;
  }

  u.n[top] = dirichlet(0);
  u.n[bottom] = dirichlet(0);

  b[top] = neumann(0);
  b[bottom] = neumann(0);


}

/**
   Do a couple of time steps
*/
   

event output (t = 0; t <= 100;  t += 1) {
  fprintf (stdout,"i = %d, dt = %g, t = %g\n",i, dt, t);
}

/**
   Write the last buoyancy
*/

event writeout(t = 100) {
  FILE * fp = fopen ("b.dat", "w");
  output_field({b, u}, fp); 
}
/**

## Results

~~~gnuplot Buoyancy after a couple of time steps
set xlabel 'x'
set ylabel 'y'
set pm3d map
splot 'b.dat' u 1:2:3

~~~ 

~~~gnuplot Vertical velocity after a couple of time steps
set xlabel 'x'
set ylabel 'y'
set pm3d map
splot 'b.dat' u 1:2:5

~~~ 


We started from a uniform buoyancy (b=1): the fluid should be at rest. After a couple
of timesteps, it is no longer uniform and the fluid is moving. I must have forgotten a boundary
condition somewhere...

## Comments

This points was actually already investigated in

- http://basilisk.fr/sandbox/Antoonvh/afluidatrest.c
- http://basilisk.fr/sandbox/Antoonvh/thermo.h

This behavior is due to a bad convergence of the poisson solver. We need to
decrease the tolerance and the timestep for the first iterations.

*/
