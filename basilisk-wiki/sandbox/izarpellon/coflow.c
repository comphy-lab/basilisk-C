/**
# Co-flow device 

Note that the [viscous term](/src/viscosity-embed.h) in the embedded
boundary solver assumes constant viscosity, which is not the case
here. So the results below need to be interpreted with caution.

The left figure is the norm of the velocity and the interface, right
figure shows adaptive mesh refinement.

<table>
<tr><td>
![](coflow/u.png){width="100%"}
</td><td>
![](coflow/cells.png){width="100%"}
</td></tr>
</table>

~~~gnuplot Average width of the thread close to the exit ($x > 1.5$) as a function of time
set xlabel 'Time'
set logscale y
aw = 1.*system("tail -n1 log | cut -d' ' -f2")
title(a) = sprintf("measured final average width = %.3g", a)
set grid
plot 'log' w l t '', aw t title(aw) lt 3
~~~

~~~gnuplot Final interface and measured average width close to the exit
reset
set grid
plot 'out' w l t 'interface', aw/2. lt 3 t title(aw), -aw/2. lt 3 t ''
~~~

~~~gnuplot Velocity profile at $x=1.5$
unset grid
set xlabel 'y'
set arrow from -aw/2.,0 to -aw/2.,35 nohead lt 3
set arrow from aw/2.,0 to aw/2.,35 nohead lt 3
poiseuille(x)=(6.*((x) + 0.5)*(0.5 - (x)))
Uv = 10.
set key bottom left
plot [-0.5:0.5][0:]'profile' w p t 'measured', (1 + 2*Uv)*poiseuille(x) t '(1 + 2*Uv)*poiseuille(y)', 1 + 3*Uv lt 4
~~~

As illustrated in the figure the outflow profile is well matched by the
Poiseuille flow
$$
u(y) = 6 (1 + 2 \text{Uv}) (y + 1 / 2) (1 / 2 - y)
$$
where $1 + 2 \text{Uv}$ is the total flow rate. The maximum velocity near the
outflow is thus
$$
u_{\max} = u (0) = 3 / 2 (1 + 2 \text{Uv})
$$
If the thread diameter is small we can assume that the exit flow rate is given
by $u_{\max} w$ and since this must match the inflow rate (unity) we get
$$
w = \frac{1}{u_{\max}} = \frac{2}{3 (1 + 2 \text{Uv})}
$$
For $\text{Uv} =$ 10 this gives $w \approx$ 0.0317 which is not far from the
measured 0.0323. 

**Exercise**: show that a more accurate approximation (valid also for a finite-width very-viscous thread) is
$$
u_{\max} = 1 + 3 \text{Uv}, \;\;\; w = \frac{1}{1 +3 \text{Uv}}
$$
which gives $w \approx$ 0.0322.
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
f[left]   = (cs[] > 0.);

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
f[bottom]   = 0.;

u.n[top] = dirichlet(- Uv*cs[]*poiseuille(x));
p[top]   = neumann(0.);
pf[top]  = neumann(0.);
f[top]   = 0.;

/**
Embedded boundary is no slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

const int maxlevel = 9;

face vector mut[];

int main()
{
  L0 = 3.000001; // just to make sure that the solid boundaries do not exactly coincide with cell boundaries
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
  solid (cs, fs,  union (intersection(- y + 1./2., + y + 1./2.),
			 intersection(- x + Hv/2. , + x + Hv/2.)));

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

event adapt (t = 2./Uv; i++) {
  const double mmax = 1e-2, fmax = 1e-2, uemax = 0.1;
  adapt_wavelet ({cs, f, u}, (double[]){mmax, fmax, uemax, uemax}, maxlevel);
}

/**
We measure the average width of the thread close to the exit. */

event thread_width (i++)
{
  double s = 0., sf = 0.;
  foreach (reduction(+:s) reduction(+:sf))
    if (x > L0/2.)
      s += dv(), sf += f[]*dv();
  fprintf (stderr, "%g %g\n", t, sf/s);
}

/**
We make images. */

event movie (t = 4./Uv)
{
  view (tx = -0.118, ty = -0.000, tz = -3.038,
	width = 918, height = 890);

  draw_vof (c = "f", lc = {1,1,1}, lw = 2);
  draw_vof (c = "cs", s = "fs", edges = true, filled = -1, fc = {1,1,1});
  squares (color = "sqrt(u.x*u.x + u.y*u.y)", spread = -1, linear = true);
  box ();
  //  vectors (u = "u", scale = 0.001);
  save ("u.png");

  
  draw_vof (c = "f", lc = {1,1,1}, lw = 2);
  draw_vof (c = "cs", s = "fs", edges = true, filled = -1, fc = {1,1,1});
  squares (color = "sqrt(u.x*u.x + u.y*u.y)", spread = -1, linear = true);
  box ();
  cells ();
  save ("cells.png");

  /**
  We also save the shape of the thread. */

  output_facets (f, stdout);
  
  /**
  And a velocity profile. */
  
  FILE * fp = fopen ("profile", "w");
  for (double y = -0.5; y < 0.51; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, L0/2., y));
  fclose (fp);
}
