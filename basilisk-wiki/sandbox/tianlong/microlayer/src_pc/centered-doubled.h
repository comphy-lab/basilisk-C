/**
# Double Pressure-Velocity Coupling Model

This is a copy of the centered solver, adding the suffix *ext*
to the name of the quantities to be resolved. We want to find
a divergence-free velocity field for the advection of the volume
fraction field in phase change simulations.
Using this approach, we solve the pressure-velocity coupling twice,
first including the volume expansion term in the continuity equation
[navier-stokes/centered-evaporation.h](centered-evaporation.h), and
then with the classical incompressible continuity equation. This
method works well also for droplets in static conditions, with strong
density ratio or with high vaporization rate, since the field velocity
which contains the Stefan flow, and the velocity used for the advection
of the volume fraction field *ufext* are decoupled.

# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}^E+\nabla\cdot(\mathbf{u}^E\otimes\mathbf{u}^E) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u}^E = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u}^E + (\nabla\mathbf{u}^E)^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar pext[];
vector uext[], gext[];
scalar pfext[];
face vector ufext[];

/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

extern scalar stefanflowext;

trace
mgstats project_sfext (struct Project q)
{
  face vector uf = q.uf;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  /**
  We add the volume expansion contribution. */

  foreach() {
    div[] += stefanflowext[]/dt;
#ifdef VARPROP
    div[] += drhodt[]/dt;
#endif
  }

#if USE_MY_SOLID
  foreach()
  {
    div[] *= (1.0 - is_solid[]);
  }

  face vector alpha2[];
  foreach_face()
  {
    alpha2.x[] = (1.0 - is_solid_face.x[]) * alpha.x[];
  }
  alpha = alpha2;
#endif 
  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats mgp = poisson (p, div, alpha,
       tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

#if USE_MY_SOLID
  boundarySolidNeummanNoauto(p);
#endif 

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
  boundary((scalar *){uf});
  return mgp;
}

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

mgstats mgp_ext, mgpf_ext, mgu_ext;

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

pext[right] = neumann (neumann_pressure(ghost));
pext[left]  = neumann (- neumann_pressure(0));

#if AXI
ufext.n[bottom] = 0.;
ufext.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                                // is zero on the axis of symmetry
pext[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
pext[top]    = neumann (neumann_pressure(ghost));
pext[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
pext[front]  = neumann (neumann_pressure(ghost));
pext[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient_ext (Point point, scalar p, coord * g)
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

  pext.nodump = pfext.nodump = true;
  
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
  ufext.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  ufext.x.refine = refine_face;
  foreach_dimension()
    ufext.x.prolongation = refine_embed_face_x;
  for (scalar s in {pext, pfext, uext, gext}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {pext, pfext})
    s.embed_gradient = pressure_embed_gradient_ext;
#endif // EMBED
#endif // TREE
}

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  trash ({ufext});
  foreach_face()
    ufext.x[] = fm.x[]*face_value (uext.x, 0);
  boundary((scalar *){ufext});
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
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
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

void prediction_ext()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (uext.x.gradient)
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
	  du.x[] = uext.x.gradient (uext.x[-1], uext.x[], uext.x[1])/Delta;
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
	  du.x[] = (uext.x[1] - uext.x[-1])/(2.*Delta);
    }

  trash ({ufext});
  foreach_face() {
    double un = dt*(uext.x[] + uext.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    ufext.x[] = uext.x[i] + (gext.x[] + gext.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = uext.y[i] < 0. ? uext.x[i,1] - uext.x[i] : uext.x[i] - uext.x[i,-1];
      ufext.x[] -= dt*uext.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = uext.z[i] < 0. ? uext.x[i,0,1] - uext.x[i] : uext.x[i] - uext.x[i,0,-1];
      ufext.x[] -= dt*uext.z[i]*fzz/(2.*Delta);
    }
    #endif
    ufext.x[] *= fm.x[];
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
    prediction_ext();
    boundary((scalar *){ufext, alpha});
    mgpf = project_sfext (ufext, pfext, alpha, dt/2., mgpf_ext.nrelax);
    advection ((scalar *){uext}, ufext, dt, (scalar *){gext});
#if USE_MY_SOLID
    boundarySolidVelCNoauto(uext);
#endif
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction_ext (double dt)
{
  foreach()
    foreach_dimension()
      uext.x[] += dt*gext.x[];
#if USE_MY_SOLID
  boundarySolidVelCNoauto(uext);
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
    correction_ext (dt);
    mgu = viscosity (uext, mu, rho, dt, mgu_ext.nrelax);
    correction_ext (-dt);
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
  trash ({ufext});
  boundary((scalar *){ufext, a});
  foreach_face()
    ufext.x[] = fm.x[]*(face_value (uext.x, 0) + dt*a.x[]);
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient_ext (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gfext[];
  foreach_face()
    gfext.x[] = fm.x[]*a.x[] - alpha.x[]*(pext[] - pext[-1])/Delta;
  boundary((scalar *){gfext});
  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({gext});
  foreach()
    foreach_dimension()
      gext.x[] = (gfext.x[] + gfext.x[1])/(fm.x[] + fm.x[1] + SEPS);

#if USE_MY_SOLID
  foreach()
  {
    if(is_solid[] == 1.0)
    {
      foreach_dimension()
        gext.x[] = 0.0;
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
  boundary((scalar *){ufext, alpha});
  mgp_ext = project_sfext (ufext, pext, alpha, dt, mgp_ext.nrelax);
  centered_gradient (pext, gext);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction_ext (dt);
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
    if (ufext.x[] && !fs.x[])
      ufext.x[] = 0.;
#endif
  event ("properties");
}
#endif

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
