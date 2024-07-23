/**
# Lock exchange (Kelvin--Helmoltz shear instability)

This is the [centered Navier--Stokes solver](/src/navier-stokes/centered.h)
version of the [multilayer lock-exchange test](kh.c) where a more
detailed description can be found. 

![Animation of the density field](kh-ns/T.mp4)

~~~gnuplot Ten equally-spaced contours of density at $t=$ 21.
set term svg enhanced size 1000,170 font ",10"
set size ratio -1
unset key
unset ytics
plot 'log' u 1:2 w l lc rgbcolor "black"
~~~
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "navier-stokes/perfs.h"

/**
We add a tracer and define the acceleration of gravity. */

scalar T[], * tracers = {T};
double G = 9.81;

/**
The acceleration field is allocated to store the Boussinesq buoyancy
term. */

face vector av[];

int main()
{

  /**
  We define a [rectangular domain using (eight) parallel
  domains](/src/Tips#non-cubic-domains). */
  
  dimensions (ny = 1);
  L0 = npe();
  X0 = -L0/2.;
  N = 64*npe();

  /**
  Acceleration field, maximum timestep and viscosity. */
  
  a = av;
  DT = 0.03;
  const face vector muc[] = {2e-4,2e-4};
  mu = muc;

  /**
  This is required due to a problem with the automatic increase in the
  number of relaxations in [poisson.h](/src/poisson.h#mg_solve). */

  // fixme
  TOLERANCE = HUGE;
  NITERMIN = 4;
  
  run();

  /**
  The contour lines are saved in 'log' to serve as reference
  solution. */
  
  system ("gnuplot -e 'set table' kh-ns.plot | sed '/^#.*/d' > log");
}

/**
We add diffusion for the tracer, with a Peclet number unity. */

event tracer_diffusion(i++) {
  diffusion (T, dt, mu);
}

/**
This is the Bousinesq buoyancy vertical acceleration, with a relative
density ratio $\Delta\rho/\rho = 10^{-2}$. */

event acceleration (i++)
{
  double dr = 0.01;
  foreach_face(y)
    av.y[] = G*dr*(T[] + T[0,-1])/2.;
}

/**
The initial conditions with a hyperbolic tangent initial profile, and
limited tracer gradient. */

event init (i = 0)
{
  T.gradient = minmod;
  foreach()
    T[] = - 0.5*tanh((x + 0.1*cos(pi*y/2.))/0.04);
}

/**
Density field at $t=$ 21. */

event density (t = 21)
{
  output_field ({T}, box = {{-L0/2.,0},{L0/2,1}});
}

/**
We generate an animation of the density field. */

event movie (t += 1; t <= 30)
{
  output_ppm (T, file = "T.mp4", spread = -1, linear = true, n = 128*npe(),
	      box = {{-L0/2.,0},{L0/2,1}});
}
