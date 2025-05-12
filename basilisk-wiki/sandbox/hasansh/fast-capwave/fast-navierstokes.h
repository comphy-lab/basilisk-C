/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations using the "constant density"
method of [Dodd and Ferrante,
2014](http://dx.doi.org/10.1016/j.jcp.2014.05.024)
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(\mu\nabla\mathbf{u})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$

The scheme implemented here is close to the method implemented in
[Basilisk](http://basilisk.fr/src/navier-stokes/centered.h) except for
the projection of velocity field.

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "viscosity.h"
#include "fast-poisson.h"

/**
The primary variables are the centered pressure field $p$, the
auxiliary centered pressure field $pn$ for storing pressure at
prevoius timestep, and the centered velocity field $\mathbf{u}$. The
centered vector field $\mathbf{g}$ will contain pressure gradients and
acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
scalar pn[];
vector u[], g[];
scalar pf[];
face vector uf[];

/**
We define an enumerator for determining projection method. For FASTPn and
FASTP* methods we need $\rho_0$ which is the menimum of $\rho_1$ and
$\rho_2$. If not specified by the user, this variable is set to one. */

enum { unsplit, fastpn, fastpstar};
double rho0 = 1.0;
int poisson_solver = unsplit;

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

(const) face vector mu = zerof, a = zerof;
(const) face vector alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

p[right] = neumann(a.x[ghost]*fm.x[ghost]/alpha.x[ghost]);
p[left]  = neumann(-a.x[]*fm.x[]/alpha.x[]);
#if AXI
uf.n[bottom] = 0.;
#else // !AXI
#  if dimension > 1
p[top]    = neumann(a.y[ghost]*fm.y[ghost]/alpha.y[ghost]);
p[bottom] = neumann(-a.y[]*fm.y[]/alpha.y[]);
#  endif
#  if dimension > 2
p[front]  = neumann(a.z[ghost]*fm.z[ghost]/alpha.z[ghost]);
p[back]   = neumann(-a.z[]*fm.z[]/alpha.z[]);
#  endif
#endif

/**
## Initial conditions */

event defaults (i = 0)
{

  /**
  The default velocity and pressure are zero. */
  
  CFL = 0.8;
  foreach() {
    foreach_dimension()
      u.x[] = g.x[] = 0.;
    p[] = pn[] = pf[] = 0.;
  }
  foreach_face()
    uf.x[] = 0.;

  /**
  The default density is one. */
  
  if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = 1;
    boundary ((scalar *){alpha});
  }

  /**
  The default acceleration is zero. */
  
  if (!is_constant(a.x)) {
    face vector av = a;
    foreach_face()
      av.x[] = 0.;
    boundary ((scalar *){av});
  }
  
  boundary ({p,pf,u,g,uf});

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;
#endif
}

/**
We initialise the face velocity field and apply boundary conditions
after user initialisation. We also define fluid properties. */

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
  boundary ((scalar *){uf});

  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
    
  event ("properties");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

double dtmax;

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last) {
  boundary ({alpha, mu, rho});
}

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension()
        du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
  else
    foreach()
      foreach_dimension()
        du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    #endif
    #if dimension > 2
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    #endif
    uf.x[] *= fm.x[];
  }
  boundary ((scalar *){uf});

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
    prediction();
    mgpf = project (uf, pf, alpha, dt/2.);
    advection ((scalar *){u}, uf, dt, (scalar *){g});
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
  boundary ((scalar *){u});  
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
    correction (dt);
    mgu = viscosity (u, mu, rho, dt);
    correction (-dt);
  }

  /**
  The (provisionary) face velocity field at time $t+\Delta t$ is
  obtained by simple interpolation. We also reset the acceleration
  field (if it is not a constant). */

  face vector af = a;
  trash ({uf,af});
  foreach_face() {
    uf.x[] = fm.x[]*(u.x[] + u.x[-1])/2.;
    if (!is_constant(af.x))
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
averaging. */

event acceleration (i++,last)
{
  boundary ((scalar *){a});
  foreach_face()
    uf.x[] += dt*fm.x[]*a.x[];
  boundary ((scalar *){uf});
}

/**
## Approximate projection

To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). */

event projection (i++,last)
{

  /**
  We use switch case for projecting with on of the usplit, fastPn,
  fastpstar methods which is specified by the user. */

  boundary ({p}); 
  switch(poisson_solver){
	case 0:
	mgp = project (uf, p,alpha, dt);
	break;
	case 1:
	mgp = fast_project (fastpn, uf, p, pn, alpha, rho0, dt);
	break;
	case 2:
	mgp = fast_project (fastpstar, uf, p, pn, alpha, rho0, dt);
	break;
  }   

  /**
  We then compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */
  
  face vector gf[];
  foreach_face()
  gf.x[] = a.x[] - alpha.x[]/fm.x[]*(p[] - p[-1])/Delta;
  boundary_flux ({gf});

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/2.;
  boundary ((scalar *){g});

  /**
  And finally add this term to the centered velocity field. */

  correction (dt);
}
