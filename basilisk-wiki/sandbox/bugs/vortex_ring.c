/**
# Vorticity evaluation seems not work in axisymmetric configuration

We compute here the evolution of an axisymmetric vortex ring using axi.h
The vorticity evaluation using vorticity(u,omega) fails in the vicinity of the axis. This seems to be linked to the way it is evaluated in utils.h
contained in the patch 9c1eb2eaeebca1dc9ea05006917a2508ee605e96 that does not work in 2D axi.

For axi, the following simplest expression may be used:

~~~c
omega[] = (u.y[1,0]-u.y[-1,0]+u.x[0,-1]-u.x[0,1]) / (2.*Delta + SEPS);
~~~


Here is the code. Note that MAXLEVEL 9 may not be sufficient but this is not the issue here. 
*/

#include "axi.h"
#include "navier-stokes/centered.h"

#define MAXLEVEL 9
#define L0 20.

/**
Parameters
*/
double Gamma=1.;
double Reynolds=100.;
double xo=0.;
double yo=1.;
double ao=0.1;
double TFINAL = 1.;

face vector muv[];

// This is necessary for convergence when lowering the tolerance
uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

/**
Main
*/

int main()
{
  size(L0);
  origin(-L0/2., 0.);
  init_grid (1 << MAXLEVEL);
  mu = muv;
  run();
}

/**
Setting viscosity and initial conditions
*/
event init (t = 0)
{
   foreach_face() {
    muv.x[] = fm.x[]*(1./Reynolds); 
    muv.y[] = fm.y[]*(1./Reynolds);
   }

  scalar psi[], omega[];

  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

/**
Vorticity initialisation and Poisson solving to get velocity
*/
   foreach() {
    omega[] = Gamma/(pi*ao*ao)*(exp(-(sq(x - xo) + sq(y-yo))/(ao*ao)) - exp(-(sq(x - xo) + sq(y+yo))/(ao*ao))) ;
    psi[] = 0.;
  }
  boundary ({psi,omega});
  poisson (psi, omega);
  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
  boundary ((scalar *){u});
}

/**
Logfile settings
*/
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/** 
Movie: the max value is set according to the theoretical viscous diffusion law.
*/
event movie (t += TFINAL/20.; t <= TFINAL) {
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, linear = true, max= Gamma/(pi*ao*ao+4*pi*t/Reynolds), min=-Gamma/(pi*ao*ao+4*pi*t/Reynolds), box = {{-4.,0.},{4,2.}}, file = "vort.mp4");
}

/**
Adaptation
*/
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}

