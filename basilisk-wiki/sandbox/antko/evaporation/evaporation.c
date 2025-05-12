/**
# Drop evaporation

We investigate the evaporation of a single droplet in a still and
relatively dry environment.

![Evaporating droplet (interface tracked in black) and
 vapour concentration field](evaporation/evap.gif)

For this investigation we use -- in an axisymmetric setting -- the VOF
description of the interface, the Bell-Collela-Glaz advection solver
and the fully implicit reaction-diffusion solver. */

#include "axi.h"
#include "advection.h"
#include "vof.h"
#include "diffusion.h"

/**
We allocate several scalars and vector fields to describe both the
interface and the concentration field. */

scalar f[];
scalar * tracers = NULL, * interfaces = {f};
scalar concentration[];

/**
We non-dimensionalise the problem with the initial radius of the drop
$R$, the diffusion coefficient of liquid vapour in air $D$ and the
saturation concentration of water in air $c_0$ (this quantity having
the dimension of a density $M\cdot L^{-3}$).

Thus the diffusion equation for the vapour in the non-dimensional
space reads:
$$
\frac{\partial c}{\partial t} = \Delta c,
$$
appropriately completed with Dirichlet conditions at the surface of
the drop and at infinity:
$$
\left\{
\begin{array}{ll}
c = 1 & \text{at the drop surface,}\\
c \to c_\infty/c_0 & \text{far from the drop.}
\end{array}
\right.
$$
As the interface recedes at a (dimensional) velocity $v_\text{e} =
\frac{D}{\rho}\nabla c \sim \frac{D}{\rho}\frac{c_0}{R} \equiv V$, a
Peclet number Pe$= \frac{VR}{D} = \frac{c_0}{\rho}$ also enters in the
problem description. Typically, for the problem of a water droplet
evaporating in dry air, this Peclet number is $O(10^{-5})$, meaning
that the problem is really dominated by diffusion. Here, we choose
Pe=$10^{-3}$, and set the relative humidity of the surroundings to
20\%: */

const double Peclet_number = 1e-3;
const double cs = 1.0;
const double cinf = 0.2;

/**
We set some numerical constants and placeholders for e.g. the
effective radius of the droplet. */

const double tend = 900.;
int LEVEL = 9, MINLEVEL = 5, MAXLEVEL = 10;
double effective_radius = 1.;
double dt, dirichlet_time_constant;

mgstats mgconcentration;

/**
Thanks to symmetry, we only solve a quarter of the domain, and
requires the concentration to drop at its asymptotic value "at
infinity" (that is, at the box boundary) */

concentration[right]  = dirichlet(cinf);
concentration[top]    = dirichlet(cinf);

/**
The main function of the program (where we set the domain geometry to
be ten times larger than the drop) */

int main()
{
  size (10.);
  origin (0., 0.);
  N = 1 << LEVEL;
  init_grid (N);
  run();
}

/**
The initial position of the interface is defined with this function */

#define circle(x,y) (sq(1.) - sq(x) -sq(y))

/**
Before the first step, we initialise the concentration field to $c_0$
within the drop, and $c_\infty$ everywhere else. The diffusion
coefficient is set to one, but beware that in Basilisk the diffusion
coefficient also incorporates information about the metrics, namely
the face length fm.x(): */

event init (t = 0) {
  fraction (f, circle(x,y));
  foreach()
    concentration[] = f[]*cs + (1 - f[])*cinf;
  boundary({concentration});
  CFL = 0.2;
}

/**
A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */

coord normal (Point point, scalar c)
{
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

/**
The core routine of the code. Several key operations are performed here:

* Computation of the concentration field,
* Derivation of the concentration gradient,
* Setting of the advection velocity field for the interface according
  to Fick's law.

An important point is how the immersed Dirichlet boundary condition is
handled here. To ensure that the concentration within the drop stays
at $c_0$ and is not washed away, we introduce a source term in the
diffusion equation. This source term is only active within the drop
and is as large as the difference between the actual concentration and
$c_0$ is. This modified equation therefore reads:
$$
\frac{\partial c}{\partial t} = \Delta c + \frac{f}{\tau} \left(c_0 - c\right).
$$
Here $f$ is the liquid phase fraction, and $\tau$ is a time constant,
which has to be shorter than the typical diffusive timescale in order
for the control to be effective. In the following, we set $\tau$ to a
tenth of the smallest diffusive time computed at the cell level. */

event velocity (i++) {

  /**
  The processes at play are slow for most of the simulation, and the
  solver is fully implicit so we can go full throttle and set huge
  timesteps like 10.  However in the final step the velocity is expected
  to diverge, so we have to be careful and bound the timestep: */

  dt = dtnext (timestep (u, 10.));
  dirichlet_time_constant = 1e-1*sq(L0/(1 << MAXLEVEL));

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[], concentration_flux[];
  foreach() {
    volumic_metric[]          =   cm[];
    dirichlet_source_term[]   =   cm[]*cs*f[]/dirichlet_time_constant;
    dirichlet_feedback_term[] = - cm[]*f[]/dirichlet_time_constant;
  }

  foreach_face()
    diffusion_coef.x[] = fm.x[];

  boundary({volumic_metric,dirichlet_source_term,dirichlet_feedback_term,diffusion_coef});

  mgconcentration = diffusion (concentration, dt, D = diffusion_coef,
			       r = dirichlet_source_term,
			       beta = dirichlet_feedback_term,
			       theta = volumic_metric);

  /**
  Having the concentration we derive the concentration gradient, which
  is a face vector field: */

  foreach_face()
    concentration_flux.x[] = (concentration[] - concentration[-1])/Delta;
  boundary((scalar*){concentration_flux});

  /**
  With the concentration gradient we can now advance the interface,
  following the lines drawn in
  [meanflow.c](/sandbox/popinet/meanflow.c). We define in particular
  the velocity of the interface as the product between the Peclet
  number, the local normal flux and the normal. Note that *u.x* is
  weighted by the metric (embodied in the diffusion coefficient). */

  foreach_face() {
    u.x[] = 0.;
    if (f[] > 0. && f[] < 1.) {
      coord n = normal (point, f);
      u.x[] = (n.x > 0. ?
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[1] :
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[-1]);
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      coord np = normal (neighborp(-1), f);
      u.x[] = (np.x > 0. ?
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[1] :
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[-1]);
    }
  }
  boundary((scalar*){u});
}

/**
We now write several post-processing events, e.g. to track the
effective radius of the drop: */

event compute_radius (i++) {
  stats s = statsf (f);
  effective_radius = pow(3.*s.sum, 1./3.);
}

event logfraction (t += 1.) {
  fprintf (stderr, "%f %.12f\n", t, effective_radius);
}

/**
During the lifetime of the drop, we also track the concentration and
fraction profiles. */

event logfile (t = 30; t += 10.; t < 580) {
  static FILE * fpinterf = fopen("interfaces", "w");
  output_facets (f, fpinterf);

  static FILE * fp = fopen("cprof", "w");
  static FILE * fpcresc = fopen("cprof_resc", "w");
  static FILE * fpf = fopen("fprof", "w");
  for (double y = 0.; y <= 5.; y += 0.01) {
    fprintf (fp, "%g %g\n", y,
      interpolate (concentration, 0., y));
    fprintf (fpf, "%g %g\n", y,
      interpolate (f, 0., y));
  }
  for (double y = effective_radius; y <= 5.; y += 0.01) {
    fprintf (fpcresc, "%g %g %g\n", effective_radius, y/effective_radius,
      interpolate (concentration, 0., y));
  }
  fprintf (fp, "\n");
  fprintf (fpcresc, "\n");
  fprintf (fpf, "\n");
}

#if 0
event gfsview (i += 10) {
  foreach()
    f[] = clamp(f[],0,1);
  boundary ({f});
  static FILE * fpview = popen ("gfsview2D evaporation.gfv", "w");
  output_gfs (fpview);
}
#endif

#if 1
event movies (t += 15; t <= tend) {
  static FILE * fpmov =
    popen ("gfsview-batch2D evaporation.gfv | ppm2gif > evap.gif", "w");

  // we clamp volume fraction because gfsview does not like over/undershoots
  foreach()
    f[] = clamp(f[],0,1);
  boundary ({f});
  output_gfs (fpmov);
  
  fprintf (fpmov, "Save stdout { format = PPM width = 300 height = 300}\n");
}
#endif

#if 1
event adapt (i++) {
  adapt_wavelet ({concentration}, (double[]){1e-3},
		 maxlevel = MAXLEVEL);
}
#endif

/**
# Results

## Radius evolution: the $\text{d}^2$ law

This problem was first investigated independently by Maxwell and
Langmuir. They found that quite counter-intuitively, the evaporation
mass rate does not scale with the surface of the drop, but with its
radius (as a consequence of evaporation being a diffusion-limited
process).

It follows that the squared-radius of the drop decreases linearly with
time (so-called "$\text{d}^2$ law"), according to:
$$
r^2(t) = r_0^2 - 2 \frac{D(c_0-c_\infty)}{\rho}t
$$

Here is represented the time evolution of the effective radius of the
drop, compared with the theoretical prediction:

~~~gnuplot Time evolution of the squared radius $r^2(t)$.
plot [0:600] [-0.1:1] "log" u 1:($2*$2), 1.-0.002*0.8*x lw 2
~~~

## Concentration profiles

Now checking the concentration profiles vs time we have

~~~gnuplot Raw concentration profiles vs time
plot "cprof" w l
~~~

The theoretical concentration profiles follow:
$$
c(r,t) = \left(\frac{R(t)}{r}\right) (c_0-c_\infty) + c_\infty
$$
This suggests to replot the concentration profiles versus the rescaled
variable $r/R(t)$:

~~~gnuplot Tentative rescaled concentration profiles
plot "cprof_resc" u 2:3, 0.8*(1./x) + 0.2 lw 2
~~~

This is not entirely satisfying.

The reason is that our infinity is 10; confinement effects, if not
dominant, are nonetheless not entirely negligible. Working out these
effects, we get:
$$
c(r,t) = c_0 + \frac{c_0-c_\infty}{1-\frac{R(t)}{R_\infty}}
\left(\frac{R(t)}{r} - 1\right)
$$

suggesting a slightly different rescaling:

~~~gnuplot Better rescaling for the concentration profiles, accounting for confinement effects:
plot "cprof_resc" u 2:(($3-1.)*(1-$1/10.)+0.8), 0.8*(1./x) lw 2
~~~

But wait. If taking confinement effects into account yield more neat
agreement, we should modify our $\text{d}^2$ law accordingly? Let's
give it a try. The confined $\text{d}^2$ law reads:
$$
R^2(t) - \frac{2}{3R_\infty} R^3(t) = R_0^2 - \frac{2}{3R_\infty}
R_0^3 - 2 \frac{D(c_0-c_\infty)}{\rho}t
$$

~~~gnuplot A "better" $\text{d}^2$ law
plot [0:600]"log" u 1:($2*$2-2./3./10.*$2*$2*$2), 1.-2./30.-0.002*0.8*x
~~~

That's it!

# References

~~~bib
@article{maxwell1878,
  title={Diffusion},
  author={Maxwell, James Clerk},
  journal={Encyclopedia britannica},
  volume={7},
  pages={214--221},
  year={1878}
}

@article{langmuir1918,
  title={The evaporation of small spheres},
  author={Langmuir, Irving},
  journal={Physical review},
  volume={12},
  number={5},
  pages={368},
  year={1918}
}
~~~
*/
