/**
# Evaporation of the mixture drop

We investigate the evaporation of a single droplet of a mixture in a still and
dry environment. In this example the solute diffuse very slowly and limits the
evaporation rate.

![Evaporating droplet with the vapor concentration
field](high_peclet_drop/video_solute.mp4)

We didn't put surface tension and we didn't make depend the equilibrium vapor
concentration on the curvature. Therefore, the interface is completely unstable
and the smallest wavelength (the grid size) appears. I didn't investigate
this instability, but we can understand it easily. The solvent, which diffuses
slowly, has to reach of the interface to evaporate. Where the interface is
closer to the core of the drop, the solvent concentration is higher, so the
evporation is quicker.

Options: */

#define VIDEO 1
#define ADAPT 1

/**
Geometry and resolution */

#define R0 1.
#define L 10. // size of the box

#define MIN_LEVEL 5
#define LEVEL 8
#define MAX_LEVEL 9
#define dR_refine (2.*L0/(1 << LEVEL))

#define MY_TOLERANCE 1e-6
#define F_ERR 1e-10
#define ADAPT_f 1e-3
#define ADAPT_solute 1e-2
#define ADAPT_vapor 1e-1

#define T_END 750.
#define DT_MAX 0.05
#define DELTA_T (T_END/200.) // for videos

/**
Includes */

#include "../advection_Q.h"
#include "../elementary_body.h"
#include "../mixtures.h"

#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

#ifndef AXI
  #define AXI 0
#endif

/**
Physical parameters

In the liquid:

* diffusion coefficient of the solute in the liquid, $D_\mathcal{L}$;
* initial concentration (molar fraction) of the solute, solute_0. */

#define D_L 1e-5
#define solute_0 0.5

/**
In the vapor: 

* the ratio of density between the vapor of solvent in the air and the density
of the liquid solvent, which corresponds to the Péclet number comparing
the evaporation velocity to the diffusion velocity of the vapor, vapor_peclet;
* diffusion coefficient of the solvent vapor in the air, $D_\mathcal{V}$;
* the saturation concentration of vapor when the solvent is pure, vcs;
* the boundary condition of the solvent vapor, c_\mathrm{inf};
* the time factor we need in [elementary_body.h](../elementary_body.h). */

#define vapor_peclet 1e-3
#define D_V 1.
#define vcs 1.0
#define cinf 0.
#define dirichlet_time_factor 10.

/**
We allocate several scalars and vector fields to describe both the
interface and the concentration field. */

scalar f[], solute[], vapor[];
scalar * tracers = {vapor}, * interfaces = {f};

/**
During the step of advection of the interface (because of the evaporation), the
solute is concentrated. If we do not take care, at equal evaporative rate, the
less liquid fraction is in the cell, the more the solute will be concentrated.
This not physical: a cell should be concentrated more quickly than its neighbor
just because its liquid fraction is smaller. To solve this problem we define a
function to compute a concentration at the interface which average the
concentration over a reconstructed cell, which have always the same volume (or
area in 2D). This reconstructed cell share its volume between the cell at the
interface and its neighbors accondingly to the normal to the interface. */

double interfacial_concentration (Point point, scalar f, scalar solute) {
  double local_solute = solute[]*f[];
  coord n = interface_normal (point, f);
  foreach_dimension() {
    int l = (n.x > 0. ? -1 : 1);
    local_solute += fabs(n.x)*solute[l]*(1. - f[]);
  }
  return local_solute;
}

/**
Thanks to symmetry, we only solve a quarter of the domain, and
requires the concentration to drop at its asymptotic value *at
infinity* (that is, at the box boundary) */

vapor[right] = dirichlet(cinf);
vapor[top] = dirichlet(cinf);

solute[right] = dirichlet(0.);
solute[top] = dirichlet(0.);

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the drop: */

int main() {

  size (L);
  origin (0., 0.);

  N = 1 << LEVEL;
  init_grid (N);
  TOLERANCE = MY_TOLERANCE;
  NITERMIN = 2;
  DT = DT_MAX;

  run();
}

/**
The initial position of the interface is defined with this function */

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

/**
Before the first step, we initialise the concentration field to $c_0$
within the drop, and $c_\infty$ everywhere else: */

double initial_solute_mass;

event init (i = 0) {
#if TREE && ADAPT
  refine (level < MAX_LEVEL && circle(x, y, (R0 - dR_refine)) < 0.
          && circle(x, y, (R0 + dR_refine)) > 0.);
#endif
  fraction (f, circle(x, y, R0));
  foreach() {
    solute[] = solute_0*f[];
    vapor[] = f[]*vcs + (1. - f[])*cinf;
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({f, solute, vapor, uf});
  initial_solute_mass = (AXI ? 4.*pi : 1.)*statsf(solute).sum;

  CFL = 0.2;
}

/**
## Evaporation velocity

The velocity due to evaporation is computed in the *stability()* event to take
into account this velocity in the CFL condition. But before, we declare two
placeholders to save both the velocity field without evaporation and the
evaporation velocity. */

face vector uf_save[], evaporation_velocity[];

event stability (i++) {
  if (i > 0) {
    vapor.D = D_V;
    vapor.peclet = vapor_peclet;
    vapor.inverse = true;
    phase_change_velocity (f, vapor, evaporation_velocity);
    foreach_face() {
      uf_save.x[] = uf.x[];
      uf.x[] += evaporation_velocity.x[];
    }
  }
  boundary((scalar *){uf, evaporation_velocity});
}

static scalar * interfaces_save = NULL;

event vof (i++) {

  /**
  ## Evaporation advection
  
  We save the VOF field describing the interface before the advection due to
  the evaporation. */
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  
  /**
  We advect just the VOF tracer alone with only the evaporation velocity. */
  
  foreach_face()
    uf.x[] = evaporation_velocity.x[];
  boundary ({f, solute, uf, vapor});
  
  f.tracers = NULL;
  
  vof_advection({f}, i);

  boundary ({f});

  /**
  We redistribute now the solute remaining in the dry cells. */
  
  distribution (f, solute, uf);
  
  /**
  Since we consider the concentration field as a fraction, it can't exceed 1.
  Therefore, we redistribute the solute of the overloaded cells. */
  
  distribution_over_2 (f, solute);
  
  /**
  We set $\mathbf{u_f}$ back to the incompressible Navier-Sokes velocity. */
  
  foreach_face()
    uf.x[] = uf_save.x[];
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, solute, uf});
  
  /**
  ## Incompressible advection
  
  To advect the solute with respect to the incompressible Navier-Stokes velocity,
  we use the scheme implemented in [vof.h](/src/vof.h). To do so, we associate
  the solute to the VOF tracer $f$. */
  
  solute.inverse = false;
  solute.gradient = minmod2;
  f.tracers = {solute};

  vof_advection({f}, i);

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, solute});

  /**
  We set the list of interfaces to NULL to prevent vof.h to do a third
  advection. */

  interfaces_save = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces_save;
}

event tracer_diffusion(i++) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  /**
  ## Diffusion of the solute with a no flux condition */
  
  foreach()
    solute[] = (f[] > F_ERR ? solute[]/f[] : 0.);
  boundary({solute});
  solute.D = D_L;
  solute.inverse = false;
  no_flux_diffusion (solute, f, dt);
  
  /**
  ## Diffusion of the vapor with Dirichlet condition
  
  Following the Raoult law, the equilibrium concentration of vapor at the
  interface depends on the concentration of solvent. We deduce the concentration
  of solvant from the concentration of solute, since we have consider that the
  concentration of solute was a molar fraction. */
  
  scalar solvant[];

  foreach() {
    if (interfacial(point, f)) {
      double local_solute = interfacial_concentration (point, f, solute);
      solvant[] = clamp(1. - local_solute, 0., 1.);
    }
    else
      solvant[] = (f[] <= F_ERR ? 0. : clamp(1. - solute[], 0., 1.));
  }
  
  /**
  
  We don't need anymore the *concentration* field of solute, so we multiply
  it back by $f$ to recover the *quantity* field we need for the next advection.
  */
  
  foreach() {
    solute[] *= f[];
    vapor[] = (f[] > F_ERR ? vcs*solvant[] : vapor[]);
  }
  boundary({solvant, solute, vapor});

#if TREE && ADAPT
  int max_level = MAX_LEVEL;
#else
  int max_level = LEVEL;
#endif

  vapor.D = D_V;
  vapor.tr_eq = vcs;
  vapor.inverse = true;

  dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor,
                       tr_op = solvant);
}

#if TREE && ADAPT
event adapt (i++) {
  adapt_wavelet ({f, solute, vapor}, (double[]){ADAPT_f, ADAPT_solute, ADAPT_vapor},
		             minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif

/**
We now juste have to write several post-processing events. */

event outputs (t = 0.; t += max(DELTA_T, DT); t <= T_END) { 

  double effective_radius, total_solute, drop_mass;

  drop_mass = statsf(f).sum;
  total_solute = statsf(solute).sum;
#if AXI
  effective_radius = pow(3.*drop_mass, 1./3.);
  total_solute *= 4.*pi;
  drop_mass *= 4.*pi;
#else
  effective_radius = pow(drop_mass/pi, 1./2.);
#endif

  /**
  We save some outputs. */

  static FILE * fp = fopen("outputs", "w");
  fprintf (fp, "%.17g %.17g %.17g %.17g %.17g\n",
           t, effective_radius, drop_mass, initial_solute_mass, total_solute);
  fflush(fp);

#if VIDEO
  scalar solute_c[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    solute_c[] = (f[] > F_ERR ? solute[]/f[] : 0.);
  }
  boundary({f, solute_c});
  
  /**
  We save a video with the *concentration* value of the solute. */

  view (fov = 6, width = 640, height = 640, samples = 1, relative = false,
        tx = 1e-3, ty = 1e-3, bg = {BG, BG, BG});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = -1,
           fc = {1., 1., 1.});
  squares ("solute_c", min = -1., max = 1., linear = false, map = cool_warm);
  mirror (n = {1., 0., 0.}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = -1,
             fc = {1., 1., 1.});
    squares ("solute_c", min = -1., max = 1., linear = false, map = cool_warm);
  }
  mirror (n = {0., 1., 0.}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = -1,
             fc = {1., 1., 1.});
    squares ("solute_c", min = -1., max = 1., linear = false, map = cool_warm);
    mirror (n = {1., 0., 0.}, alpha = 0.) {
      draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = -1,
               fc = {1., 1., 1.});
      squares ("solute_c", min = -1., max = 1., linear = false,
      map = cool_warm);
    }
  }
  save ("video_solute.mp4");
#endif 
}
