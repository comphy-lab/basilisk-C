/**
I am using Ubuntu 20.04LTS and I updated Basilisk before this.

I tried manually setting the density to 1 (which I think is the
default value) in examples/karman.c and including it in the expression
for the viscosity, but I experinced a problem.

I added to the top of the file:\
const double rho_value = 1.;\
const scalar rhoc[] = rho_value;

and to main():\
rho = rhoc;

In "event properties", when I added rhoc[] to the expression for the
viscosity:\
muv.x[] = fm.x[]\*rhoc[]\*0.125/Reynolds;\
I think I did not notice any difference. But, when I added rho[]:\
muv.x[] = fm.x[]\*rho[]\*0.125/Reynolds;\
I think that the difference is obvious.

According to this, "rho" is not equal to "rhoc" and hence is not equal
to 1? Also, the value for the density that I manually set is different
than the default one?

I also added a movie for "rho". I also saw artefacts similar to the
ones I previously reported in:\
[http://basilisk.fr/sandbox/bugs/bug_karman.c](http://basilisk.fr/sandbox/bugs/bug_karman.c)\
so maybe the artefacts in "omega" are caused by artefacts in "rho".

Did I make a mistake or did I encounter a bug? Based on this and on
the previous report about setting the viscosity:\
[http://basilisk.fr/sandbox/bugs/porous3D_mu.c](http://basilisk.fr/sandbox/bugs/porous3D_mu.c)\
results obtained with Basilisk could be questionable. It looks to me
that both the density and the viscosity cannot be set to desired
values. */

/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];
const double rho_value = 1.; //ADDED///////////////////////////////////////////////////////////////
const scalar rhoc[] = rho_value; //ADDED///////////////////////////////////////////////////////////

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  
  rho = rhoc; //ADDED//////////////////////////////////////////////////////////////////////////////

  /**
  When using bview we can interactively control the Reynolds number
  and maximum level of refinement. */
  
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
  
  run(); 
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
//    muv.x[] = fm.x[]*0.125/Reynolds; //DEFAULT///////////////////////////////////////////////////
//    muv.x[] = fm.x[]*rhoc[]*0.125/Reynolds; //OPTION_1///////////////////////////////////////////
    muv.x[] = fm.x[]*rho[]*0.125/Reynolds; //OPTION_2////////////////////////////////////////////
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{

  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */

  solid (cs, fs, intersection (intersection (0.5 - y, 0.5 + y),
			       sq(x) + sq(y) - sq(0.125/2.)));

  /**
  We set the initial velocity field. */
  
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 15.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
  output_ppm (rho, file = "rho.mp4", box = {{-0.5,-0.5},{7.5,0.5}}, //ADDED////////////////////////
	      min = 0, max = 2, linear = true, mask = m); //ADDED//////////////////////////////////
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/