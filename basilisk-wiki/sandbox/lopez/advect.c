/**
# Integration of the advection equation

This test illustrates how to calculate the advection transport equation,

$$
\partial_t c + \nabla \cdot (\mathbf{u} c) = 0 \,
$$

using Basilisk being $c$ any tracer. Note that the above equation
only coincides with 

$$ 
\partial_t c + \mathbf{u} \cdot \nabla c = 0 \,
$$ 

if the velocity field is solenoidal ($\nabla \cdot \mathbf{u} =
0$).  

Basilisk provides two procedures to calculate the first equation:
[advection.h]() that it is based in the scheme of
[Bell-Collela-Glaz,1989](references.bib#bell89) and [vof.h]() in which
the face flux calculations of the tracer $c$ (associated to a volume
of fraction $f$) are based in the geometrical fluxes of $f$. We will
test the following cases,

* Case 1. 2D with velocity field $u_y=1$. Equation solved: $\partial_t c +
 \partial_y (u_y c) = 0$. Exact solution: $c(y,t) = e^{y-t}$.  

* Case 2. 2D with velocity field $u_y=1/y$. Equation solved:
 $\partial_t c + \partial_y (u_y c) =\partial_t c + \partial_y (c/y) =
 0$. Exact solution: $c(y,t) = y\, e^{y^2/2-t}$. Note that the equation
 is different to $\partial_t c + u \partial_y c = 0$ because the
 velocity field is not divergence free.

* Case 3. Axisymmetric with velocity field $u_r=1/y$. Equation solved:
 $\partial_t c + y^{-1} \partial_y (y u c) = 0$. Exact solution:
 $c(y,t) = e^{y^2/2-t}$. Note that the equation is equal to $\partial_t c + u \partial_y c = 0$ 
 (this velocity field is solenoidal in cylindrical coordinates). */

#define CASE 3

#if CASE == 1
#define exact(r,t) (exp(r-t))
#define veloc(r) (1)
#elif CASE == 2
#define exact(r,t) (r*exp(0.5*sq(r)-t))
#define veloc(r) (1./r)
#else
#include "axi.h"
#define exact(r,t) (exp(0.5*sq(r)-t))
#define veloc(r) (1./r)
#endif

#include "advection.h"
#include "vof.h"
scalar c[], f[], cv[];
scalar * tracers = {c};
scalar * interfaces = {f};

/**
To avoid singularities in the velocity field the y-coordinates span
between 1 and 2. */

int main()
{
  X0 = -0.5;
  Y0 = 1.0;
  L0 = 1.0 [0];
  DT = 1.0 [0];
  N = 64;
  f.tracers = {cv};
  run();
}

c[bottom] = dirichlet(exact(1,t));
cv[bottom] = dirichlet(exact(1,t));

event init (i = 0) {

  fraction(f, 1.2-y);
  foreach() {
    c[] = exact(y,0);
    cv[] = exact(y,0)*f[];
  }

  /**
  Important! The face vector _u_ is the velocity weighted with the
  face metric factor _fm_ to take into account properly all the metric
  terms. This prevention is unnecessary when the tracers are solved
  together with the fluid motion since the velocity-metric weighting is
  performed within the [navier-stokes/centered.h](Navier-Stokes
  solver). */

  foreach_face(y)
    u.y[] = veloc(y)*fm.y[];

  /**
  vof.h does not set the time step "per se". It leaves this task to
  other libraries (for example the [navier-stokes/centered.h](Navier-Stokes solver)).
  Therefore we have to compute it. */

  dt = dtnext (timestep (u, DT));
}

/**
## Results

We just compare analytical with numerical results in some instants... */

event profile (t += 0.4) {
  for (double y = 1.; y <= 2.; y += 1e-2)
    fprintf (stderr, "%g %.4f %.4f %.4f\n",
	     y, interpolate (c, 0., y), 
	     interpolate (cv, 0., y), exact(y,t));
}

event salvado(t = 1.2) {
  dump();
}
/**
~~~gnuplot Concentration profiles for several instants.
set terminal @PNG enhanced size 640,640 font ",12"
set xlabel 't'
set ylabel 'c[], cv[]'
plot 'log' u 1:2  t 'advection.h', 'log' u 1:3  t 'vof.h', 'log' u 1:4  lt 4 t 'Analytical'
~~~
*/
