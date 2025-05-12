/**
# Incompressible Navier--Stokes solver (MAC formulation)

We wish to approximate numerically the incompressible Navier--Stokes
equations with phase change
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
-\nabla p + \nabla\cdot(\nu\nabla\mathbf{u})
$$
$$
\nabla\cdot\mathbf{u} = \dot{m} \left(\dfrac{1}{\rho_g}
- \dfrac{1}{\rho_l}\right)\delta_\Gamma
$$

We will use the generic time loop, a CFL-limited timestep and we will
need to solve a Poisson problem. */

#include "run.h"
#include "timestep.h"
#include "poisson.h"

/**
The Markers-And-Cells (MAC) formulation was first described in the
pioneering paper of [Harlow and Welch,
1965](/src/references.bib#harlow1965). It relies on a *face*
discretisation of the velocity components `u.x` and `u.y`, relative to
the (centered) pressure `p`. This guarantees the consistency of the
discrete gradient, divergence and Laplacian operators and leads to a
stable (mode-free) integration. */

scalar p[];
vector u[];
face vector uf[];

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;

/**
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;

/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

trace
mgstats project_sf (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  
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

  foreach()
    div[] += stefanflow[]/dt;

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

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}

/**
The only parameter is the viscosity coefficient $\nu.y[-1,0]*nu$.

The statistics for the (multigrid) solution of the Poisson problem are
stored in `mgp`. */

face vector nu[];
mgstats mgp;

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] : \
                             a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

double dtmax;

event init (i = 0)
{
  trash ({uf});

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

### Advection--Diffusion 

In a first step, we compute $\mathbf{u}_*$
$$
\frac{\mathbf{u}_* - \mathbf{u}_n}{dt} = \nabla\cdot\mathbf{S}
$$
with $\mathbf{S}$ the symmetric tensor
$$
\mathbf{S} = - \mathbf{u}\otimes\mathbf{u} + \nu\nabla\mathbf{u} =
\left(\begin{array}{cc}
- u_x^2 + 2\nu\partial_xu_x & - u_xu_y + \nu(\partial_yu_x + \partial_xu_y)\\
\ldots & - u_y^2 + 2\nu\partial_yu_y
\end{array}\right)
$$ 

The timestep for this iteration is controlled by the CFL condition
(and the timing of upcoming events). */

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
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

event advance (i++,last)
{
  /**
  We update the dynamic viscosity. */

  foreach_face()
    nu.x[] = mu.x[]*alpha.x[];

  /**
  We allocate a local symmetric tensor field. To be able to compute the
  divergence of the tensor at the face locations, we need to
  compute the diagonal components at the center of cells and the
  off-diagonal component at the vertices. 
  
  ![Staggering of $\mathbf{u}$ and $\mathbf{S}$](/src/figures/Sxx.svg) */
  
  scalar Sxx[], Syy[];
  vertex scalar Sxy[];

  /**
  We average the velocity components at the center to compute the
  diagonal components. */

  foreach() {
    Sxx[] = - sq(uf.x[] + uf.x[1,0])/4. + 2.*(nu.x[1,0]*uf.x[1,0] - nu.x[]*uf.x[])/Delta;
    Syy[] = - sq(uf.y[] + uf.y[0,1])/4. + 2.*(nu.y[0,1]*uf.y[0,1] - nu.y[]*uf.y[])/Delta;
  }

  /**
  We average horizontally and vertically to compute the off-diagonal
  component at the vertices. */

  foreach_vertex()
    Sxy[] = 
      - (uf.x[] + uf.x[0,-1])*(uf.y[] + uf.y[-1,0])/4. +
      (nu.x[]*uf.x[] - nu.x[0,-1]*uf.x[0,-1] + nu.y[]*uf.y[] - nu.y[-1,0]*uf.y[-1,0])/Delta;

  /**
  Finally we compute
  $$
  \mathbf{u}_* = \mathbf{u}_n + dt\nabla\cdot\mathbf{S}
  $$ */

  foreach_face(x)
    uf.x[] += dt*(Sxx[] - Sxx[-1,0] + Sxy[0,1] - Sxy[])/Delta;

  foreach_face(y)
    uf.y[] += dt*(Syy[] - Syy[0,-1] + Sxy[1,0] - Sxy[])/Delta;

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

event acceleration (i++,last)
{
  foreach_face()
    uf.x[] += dt*a.x[]*fm.x[];
}

/**
### Projection 

In a second step we compute
$$
\mathbf{u}_{n+1} = \mathbf{u}_* - \Delta t\nabla p
$$
with the condition
$$
\nabla\cdot\mathbf{u}_{n+1} = \dot{m} \left(\dfrac{1}{\rho_g}
- \dfrac{1}{\rho_l}\right)\delta_\Gamma
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\nabla p) = \frac{\nabla\cdot\mathbf{u}_*}{\Delta t}
$$ */

event projection (i++,last)
{
  mgp = project_sf (uf, p, alpha, dt, mgp.nrelax);
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

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last) {
  /**
  We reconstruct a colocated velocity using a linear
  interpolation just for visualization. */

  foreach()
    foreach_dimension()
      u.x[] = 0.5*(uf.x[1] + uf.x[]);
}

