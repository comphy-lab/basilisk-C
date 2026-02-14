/**
# Co-flow device 

![Horizontal velocity component and interface](coflow/ux.png){width="60%"}

![Adaptive mesh refinement](coflow/cells.png){width="60%"}
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

scalar f[];

#ifdef PECLET
#  include "tracer.h"
scalar * tracers = {f};
#  include "diffusion.h"
#else // !PECLET
#  include "vof.h"
scalar * interfaces = {f};
#endif

#include "view.h"

/**
## Parameters 

The velocity unit is the "horizontal" injection velocity. The length
unit is the width of the "horizontal" injection channel. */

const double Re = 0.02;    // horizontal Reynolds
const double Uv = 10;      // relative vertical velocity
const double muv = 1./20;  // relative viscosity of the vertical fluid
const double Hv = 1.;      // relative vertical channel diameter
#ifdef PECLET
const double Pe = PECLET;
#endif

/**
## Boundary conditions 

Left boundary: inflow with unit velocity and Poiseuille flow. */

#define poiseuille(x) (6.*((x) + 0.5)*(0.5 - (x)))

u.n[left] = dirichlet(cs[]*poiseuille(y));
p[left]   = neumann(0.);
pf[left]  = neumann(0.);
f[left]   = 1.;

/**
Right boundary: outflow */

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
Top and bottom boundaries: inflow with Uv velocity and Poiseuille
flow. */

u.n[bottom] = dirichlet(+ Uv*cs[]*poiseuille(x));
p[bottom]   = neumann(0.);
pf[bottom]  = neumann(0.);
f[bottom]   = dirichlet(0.);

u.n[top] = dirichlet(- Uv*cs[]*poiseuille(x));
p[top]   = neumann(0.);
pf[top]  = neumann(0.);
f[top]   = dirichlet(0);
  
/**
Embedded boundary is no slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

const int maxlevel = 9;

face vector mut[];

int main()
{
  L0 = 3.;
  DT = HUGE;
  origin (-1.25, -L0/2.);
  N = 128;
  mu = mut;
#ifdef PECLET  
  f.gradient = minmod2;
#endif

  /**
  The dirty trick below may be necessary for large viscosity
  ratios. */
  
#if 0 // c'est du bricolage à améliorer 
  TOLERANCE = HUGE;
  NITERMIN = 10;
  NITERMAX = 10;
#endif
  
  run(); 
}

/**
We set the viscosity (and diffusion coefficient). */

#ifdef PECLET
face vector CD[];
#endif

event properties (i++) {
  double muh = 1./Re, muV = muv*muh;
  foreach_face() {
    mut.x[] = fm.x[]*((muh - muV)*((clamp(f[], 0, 1) + clamp(f[-1], 0, 1))/2.) + muV);
#ifdef PECLET
    CD.x[] = fm.x[]/Pe;
#endif
  }
}

/**
The channel geometry is a cross. */

event init (t = 0) {

  /**
  Be careful here, if the channel width is exactly one, then no-slip
  boundary conditions are not applied properly... this is a bug. */
  
  solid (cs, fs,  union (intersection(- y + 0.999/2., + y + 0.999/2.),
			 intersection(- x + Hv*0.999/2. , + x + Hv*0.999/2.)));

  /**
  The initial tracer fills the inflow channel. */

  foreach()
    f[] = (cs[] > 0.)*(x < 0.5)*(fabs (y) - Delta/2. < 0.5);
}

/**
Diffusion, if necessary. */

#ifdef PECLET
event tracer_diffusion (i++)
{
  diffusion (f,dt,CD);
}
#endif // PECLET

/**
Mesh adaptation. */

#if 1
event adapt (t = 0.3; i++) {
  const double mmax = 1e-2, fmax = 1e-2, uemax = 0.1;
  adapt_wavelet ({cs, f, u}, (double[]){mmax, fmax, uemax, uemax}, maxlevel);
}
#endif

/**
We make images. */

event movie (t = 0.35)
{
  view (tx = -0.118, ty = -0.000, tz = -3.038,
	width = 918, height = 890);
  draw_vof (c = "cs", s = "fs", edges = true, filled = -1);
  squares (color = "u.x", spread = -1, linear = true);
  box ();
  //  vectors (u = "u", scale = 0.001);
  draw_vof (c = "f");
  save ("ux.png");

  
  draw_vof (c = "cs", s = "fs", edges = true, filled = -1);
  squares (color = "u.x", spread = -1, linear = true);
  box ();
  draw_vof (c = "f");
  cells ();
  save ("cells.png");
}
