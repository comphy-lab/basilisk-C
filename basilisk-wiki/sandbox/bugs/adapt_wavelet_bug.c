/**
# Rising bug

This is a slight modification done on the rising bubble example. When compiling this example with the -catch flag, the simulation crashes as soon as adapt_wavelet is called. gdb reveals that the crash happens at:  
> Program received signal SIGFPE, Arithmetic exception.  
> 0x000000000047a2f6 in _boundary1 (point=..., neighbor=..., _s=...)  
> at /home/gerris/basilisk/src/navier-stokes/centered.h:79  
> 79	p[left]  = neumann(-a.x[]*fm.x[]/alpha.x[]);  
An inspection of the fields shows that a.x is not defined near the bubble (ie in the refined region).*/

#if 0
# include "axi.h"
#endif
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

/**
Hysing et al. consider two cases (1 and 2), with the densities, dynamic
viscosities and surface tension of fluid 1 and 2 given below. */

#define rho1 1000.
#define mu1 10.

#if 1
# define rho2 100.
# define mu2 1.
# define SIGMA 24.5
#else
# define rho2 1.
# define mu2 0.1
# define SIGMA 1.96
#endif

#define LEVEL 8
#define MAXLEVEL 10

/**
The interface between the two-phases will be tracked with the volume
fraction field *f*. We allocate the fields for the variable density
and viscosity (*alphav*, *alphacv* and *muv* respectively). */

scalar f[];
scalar * interfaces = {f};
face vector alphav[];
scalar rhov[];
face vector muv[];

/**
The boundary conditions are slip lateral walls and no-slip on the right
and left walls. */

u.t[right]  = dirichlet(0);
uf.t[right] = dirichlet(0);
u.t[left]   = dirichlet(0);
uf.t[left]  = dirichlet(0);

uf.n[top]    = 0;
uf.n[bottom] = 0;
p[top]       = neumann(0);
p[bottom]    = neumann(0);

int main() {

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  $256\times 64$ grid points. */

  size (2);
  init_grid (1 << LEVEL);
  
  /**
  The density and viscosity are defined by the variable fields we
  allocated above. We also set the surface tension for interface *f*.
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  alpha = alphav;
  rho = rhov;
  mu = muv;
  f.sigma = SIGMA;
  TOLERANCE = 1e-4;

  run();
}

event init (t = 0) {

  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x - 0.5) + sq(y) - sq(0.25);
  fractions (phi, f);
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}

/**
The density and viscosity are defined using the arithmetic average. */

#define rho(f) ((f)*rho1 + (1. - (f))*rho2)
#define mu(f)  ((f)*mu1 + (1. - (f))*mu2)

event properties (i++) {
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);
}

event adapt (i = 10; i++) {
  adapt_wavelet ({f,u.x,u.y}, (double[]){1e-2,5e-3,5e-3}, MAXLEVEL);
  event ("properties");
}

/**
At $t=3$ we output the shape of the bubble. */

event interface (t = 0.2) {
  output_facets (f, stderr);
}


