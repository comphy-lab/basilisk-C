/**
The original centered-evaporation.h, [copied from Edoardo's sandbox](http://basilisk.fr/sandbox/ecipriano/README), is modified from the previous one fluid method-based solver to a ghost fluid method-based solver. In this file, two velocities for each phase are employed. For a given phase, velocities located in the other phase are meaningless. We set them to the correct value to impose the correct jump conditions. By setting the ghost velocities accordingly, for each phase's velocity, the discontinuity is removed, thus improving numerical stability.

# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations with phase change
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = \dot{m} \left(\dfrac{1}{\rho_g}
- \dfrac{1}{\rho_l}\right)\delta_\Gamma
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
# include "viscosity-embed.h"
#else
# include "myviscosity.h"
#endif

#if USE_MY_SOLID
extern scalar is_solid;
extern face vector is_solid_face;
void boundarySolidNeummanNoauto (scalar s);
void boundarySolidVelCNoauto (vector v);
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

face vector uf2[];
vector u2[];
mgstats mgu2;

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;


/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

void setCenterGhostVel(vector u, vector u2)
{
  extern scalar is_interfacial;
  extern scalar color_cc;
  extern scalar mdot;
  extern vector n_pc;
  extern double rho1, rho2;
  double dendif = (1.0 / rho1 - 1.0 / rho2);
  foreach()
  {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif
    //mdot will be zero beyond the interfacial region
    // we set the ghost velocity
    // real liquid region, set ghost vapor velocity
    if (color_cc[] > 0.0)
    {
      foreach_dimension()
      {
        u2.x[] = u.x[] + mdot[] * n_pc.x[] * dendif;
      }
    }
    else
    {
      // real vapor region, set ghost liquid  velocity
      foreach_dimension()
      {
        u.x[] = u2.x[] - mdot[] * n_pc.x[] * dendif;
      }
    }
  }

  boundary({u, u2});
}
void setFaceGhostVel()
{
  extern scalar is_interfacial;
  extern scalar color_cc;
  extern scalar mdot;
  extern vector n_pc;
  extern face vector face_fraction;
  extern double rho1, rho2;
  double dendif = (1.0 / rho1 - 1.0 / rho2);

  foreach_face()
  {
#if USE_MY_SOLID
    if ((int)is_solid_face.x[] == 1)
      continue;
#endif
    bool is_inf_m = (int)is_interfacial[-1] == 1;
    bool is_inf_p = (int)is_interfacial[] == 1;
    int num = 0;
    double mdot_f = 0.0;
    double nf = 0.0;
    if(!is_inf_m && !is_inf_p)
    {
      if (face_fraction.x[] > 0.5)
      {
        uf2.x[] = uf.x[];
      }
      else
      {
        uf.x[] = uf2.x[];
      }
      continue;
    }
    if(is_inf_m)
    {
      mdot_f += mdot[-1];
      nf += n_pc.x[-1];
      num++;
    }
    if(is_inf_p)
    {
      mdot_f += mdot[];
      nf += n_pc.x[];
      num++;
    }

    if(num == 0)
    {
      printf("EORRO!\n");
      exit(0);
    }
    
    mdot_f /= (double)num;
    nf /= (double)num;

    if (face_fraction.x[] > 0.5)
    {
      uf2.x[] = uf.x[] + fm.x[] * mdot_f * nf * dendif;
    }
    else
    {
      uf.x[] = uf2.x[] - fm.x[] * mdot_f * nf * dendif;
    }
  }

  boundary({uf, uf2});
}

struct ProjectDouble {
  face vector uf;
  face vector uf2;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats project_sf_double (struct ProjectDouble q)
{
  face vector uf = q.uf;
  face vector uf2 = q.uf2;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  extern scalar f;
  foreach() {
    div[] = 0.;
    if (f[] > 0.5)
    {
      foreach_dimension()
          div[] += uf.x[1] - uf.x[];
    }
    else
    {
      foreach_dimension()
          div[] += uf2.x[1] - uf2.x[];
    }
    div[] /= dt*Delta;
  }

  // for 1D problem, for a better convergence, we solve the poisson equaiton on the whole domain
#if USE_MY_SOLID && IS_1D
  int Ncell = 1 << grid->maxdepth;
  double *div_top = (double *)calloc(Ncell, sizeof(double));
  double *f_top = (double *)calloc(Ncell, sizeof(double));
  extern scalar f;
  extern double rho1, rho2;
  scalar f2D[];
  for (int ii = 0; ii < Ncell; ++ii)
  {
    div_top[ii] = 0.0;
    f_top[ii] = 0.0;
  }

  foreach_boundary(top, reduction(+:div_top[:Ncell]) reduction(+:f_top[:Ncell]))
  {
    div_top[point.i - 2] = div[];
    f_top[point.i - 2] = f[];
  }

  foreach ()
  {
    div[] = div_top[point.i - 2];
    f2D[] = f_top[point.i - 2];
  }
  face vector alpha2D[];
  foreach_face()
  {
    double ff = (f2D[] + f2D[-1]) / 2.;
    double rhof = ff * (rho1 - rho2) + rho2;
    alpha2D.x[] = fm.x[] / rhof;
  }
  alpha = alpha2D;
  free(div_top);
  free(f_top);
#endif

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

#if USE_MY_SOLID && !IS_1D
  mgstats mgp = poisson_solid (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
  boundarySolidNeummanNoauto(p);
#else
  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);
#endif 

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
  {
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
    uf2.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
  }

  boundary({uf, uf2});

  return mgp;
}

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
			      a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0)
{

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE
}


/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  trash ({uf});
  foreach_face()
  {
    uf.x[] = fm.x[]*face_value (u.x, 0);
    uf2.x[] = fm.x[]*face_value (u2.x, 0);
  }
  boundary({uf, uf2});
  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : min(timestep (uf, dtmax),timestep (uf2, dtmax)));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction_double(vector u, face vector uf)
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
#if USE_MY_SOLID
    if(is_solid_face.x[] || is_solid_face.x[1])
      du.x[] = 0.;
    else
#endif
	  du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
#if USE_MY_SOLID
    if(is_solid_face.x[] || is_solid_face.x[1])
      du.x[] = 0.;
    else
#endif
	  du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar *){du});
}

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++,last)
{
  if (!stokes) {
    prediction_double(u, uf);
    prediction_double(u2, uf2);
    boundary({uf, uf2, alpha});
    setFaceGhostVel();
    mgpf = project_sf_double (uf, uf2, pf, alpha, dt/2., mgpf.nrelax);
    setFaceGhostVel();
    advection ((scalar *){u}, uf, dt, (scalar *){g});
    advection ((scalar *){u2}, uf2, dt, (scalar *){g});
#if USE_MY_SOLID
    boundarySolidVelCNoauto(u);
    boundarySolidVelCNoauto(u2);
#endif
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (vector u, double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
#if USE_MY_SOLID
  boundarySolidVelCNoauto(u);
#endif
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (u, dt);
    correction (u2, dt);
    setCenterGhostVel(u, u2);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    mgu2 = viscosity (u2, mu, rho, dt, mgu.nrelax);
    correction (u, -dt);
    correction (u2, -dt);
    setCenterGhostVel(u, u2);
  }

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  trash({uf, uf2});
  boundary({uf, uf2, a});
  foreach_face()
  {
    uf.x[] = fm.x[] * (face_value(u.x, 0) + dt * a.x[]);
    uf2.x[] = fm.x[] * (face_value(u2.x, 0) + dt * a.x[]);
  }
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;
  boundary({gf});
  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);

#if USE_MY_SOLID
  foreach()
  {
    if(is_solid[] == 1.0)
    {
      foreach_dimension()
        g.x[] = 0.0;
    }
  }
#endif
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  boundary({uf, uf2, alpha});
  setFaceGhostVel();
  mgp = project_sf_double (uf, uf2, p, alpha, dt, mgp.nrelax);
  setFaceGhostVel();
  centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (u, dt);
  correction (u2, dt);
  setCenterGhostVel(u, u2);
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
#endif
  event ("properties");
}
#endif

event outLog(i++, last)
{
  // do nothing;
}

event dumpDebug(i++, last)
{
  // do nothing;
}

event movie(i++, last)
{
  // do nothing;
}

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
