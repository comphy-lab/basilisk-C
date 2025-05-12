/**
# Incompressible Navier--Stokes solver with Phase Change Jump Condition

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
for viscosity. This scheme is extended considering the phase change
between two phases described by a VOF or CLSVOF approach. Using the
jump condition is convenient because it does not require an explicit
source term, localized at the gas-liquid interface, in the projection
step, which causes oscillations in the velocity field.

The method that we use here to enforce the velocity jump condition is
the same proposed by [Tanguy et al. 2007](#tanguy2007level) for droplet
evaporation problems, and by [Tanguy et al. 2014](#tanguy2014benchmarks)
for boiling simulations. It is based on the combination of the ghost fluid
velocity method poposed by [Nguyen et al. 2001](#nguyen2001boundary),
which defines ghost velocities as:

$$
  \mathbf{u}^{ghost}_l = \mathbf{u}_g -
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}
$$

$$
  \mathbf{u}^{ghost}_g = \mathbf{u}_l +
  \dot{m}\left(\dfrac{1}{\rho_g} - \dfrac{1}{\rho_l}\right)\mathbf{n}
$$

The ghost velocities impose the continuity of the velocity fields
across the interface. Two different momentum equations for the gas
and liquid phase velocities are solved, a single projection step is
used to update the velocities at the new time step. Additional
velocity extensions are used to obtain a divergence-free velocity for
the advection of the volume fraction field.
*/

#define VELOCITY_JUMP 1

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
# include "viscosity-embed.h"
#else
# include "myviscosity.h"
#endif
#include "redistance.h"
#include "aslam.h"
#include "vof2front.h"

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

/**
Other variables specific to this algorithm. */

vector u1[], u2[];
face vector uf1[], uf2[], ufext[];
scalar mEvapTotE[];
scalar is_interfacial[];
vector n[];
scalar ls[];

extern scalar f;
extern scalar mEvapTot;
extern double rho1, rho2, mu1, mu2;

#if USE_MY_SOLID
extern scalar is_solid;
extern face vector is_solid_face;
void boundarySolidNeummanNoauto (scalar s);
void boundarySolidVelCNoauto (vector v);
#endif

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
The volume expansion term is declared in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

extern scalar stefanflow;
#ifdef VARPROP
scalar drhodt[];
#endif

/**
## Helper functions

We define the function that performs the projection
step with the volume expansion term due to the phase
change. */

struct ProjectTwo {
  face vector uf1;
  face vector uf2;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats project_sf_twofield (struct ProjectTwo q)
{
  face vector uf1 = q.uf1;
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
  foreach() {
    div[] = 0.;
    double div1 = 0., div2 = 0.;
    foreach_dimension() {
      div1 += uf1.x[1] - uf1.x[];
      div2 += uf2.x[1] - uf2.x[];
    }
    div[] = (f[] > 0.5) ? div1 : div2;
    div[] /= dt*Delta;
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
  We add the density lagrangian derivative. */

//#ifdef VARPROP
//  foreach() {
//    div[] += drhodt[]/dt;
//  }
//#endif

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

  foreach_face() {
    uf1.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
    uf2.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);
  }

  boundary((scalar *){uf1, uf2});

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
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :  \
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
uf1.n[bottom] = 0.;
uf1.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
uf2.n[bottom] = 0.;
uf2.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
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
  uf1.x.refine = refine_face_solenoidal;
  uf2.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  uf1.x.refine = refine_face;
  foreach_dimension() {
    uf.x.prolongation = refine_embed_face_x;
    uf1.x.prolongation = refine_embed_face_x;
    uf2.x.prolongation = refine_embed_face_x;
  }
  for (scalar s in {p, pf, u, g, u1, u2}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf, ps, pg})
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
  trash ({uf1, uf2, uf});
  foreach_face() {
    uf.x[] = fm.x[]*face_value (u.x, 0);
    uf1.x[] = fm.x[]*face_value (u1.x, 0);
    uf2.x[] = fm.x[]*face_value (u2.x, 0);
  }

  boundary((scalar *){uf, uf1, uf2});

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
  dt = dtnext (stokes ? dtmax : min (timestep (uf1, dtmax), timestep
        (uf2, dtmax)));
}

/**
## Extrapolations

We use PDE-based Aslam's extrapolations, to extrapolate the
vaporization rate $\dot{m}$, from the interface to the liquid and gas
phases. */
trace
void extrapolations (void)
{

  scalar color_cc[];
  foreach()
  {
    color_cc[] = f[] > 0.5 ? 1.0 : -1.0;
  }

  /**
  We store the extrapolated vaporization rate $\hat{m}$ on the fields
  `mEvapTot1` and `mEvapTot2`. Using `constant_extrapolations` the
  liquid and gas phase cells are populated with the vaporization rate
  without changing its values. */

  foreach() {
    mEvapTotE[] = mEvapTot[];
  }


  scalar is_cut[];
  //instead of using alsma, we employ a simple average strategy
  //we should label the interfacial region
  foreach()
  {
    is_interfacial[] = 0.;
    is_cut[] = 0.;
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif
    if (f[] > F_ERR && f[] < 1.0 - F_ERR)
      is_cut[] = 1.0;

    bool isn = false;

    foreach_neighbor(1)
    {
      if (f[] > F_ERR && f[] < 1.0 - F_ERR)
        isn = true;
    }
    
    is_interfacial[] = isn ? 1.0 : 0.0;
  }

  //paranoid coding
  //set the normal outside to zero
  foreach ()
  {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif

    foreach_dimension()
    {
      n.x[] *= is_interfacial[];
    }

  }

  getSignedDistance(color_cc, n, f, ls);
  extendVariables(1.0, color_cc, n, is_interfacial, is_cut, true, mEvapTotE);
  extendVariables(-1.0, color_cc, n, is_interfacial, is_cut, true, mEvapTotE);

//   //here we extend the mEvapTot to the neighbors via a simple average stragegy
//   foreach ()
//   {
// #if USE_MY_SOLID
//     if ((int)is_solid[] == 1)
//       continue;
// #endif
//     //the values inside the cut cells will not be changed
//     if (f[] > F_ERR && f[] < 1.0 - F_ERR)
//       continue;

//     int num_ave = 0;
//     double total_ave = 0;
//     if((int)is_interfacial[] == 1)
//     {
//       foreach_neighbor(1)
//       {
//         if((int)is_interfacial[] == 1 && fabs(mEvapTot[]) > 0.0)
//         {
//           total_ave += mEvapTot[];
//           num_ave ++;
//         }
//       }
//     }
//     mEvapTotE[] = num_ave == 0 ? 0.0 : total_ave / (double) num_ave;
//   }
}

/**
We perform the extrapolations after the `phasechange` event in
[evaporation.h](/sandbox/ecipriano/src/evaporation.h). */

event phasechange (i++,last);

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last)
{
  extrapolations();
}
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);


void setGhostVelCenter()
{
  double den_diff = 1./rho2 - 1./rho1;
  foreach()
  {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif
    //mdot will be zero beyond the interfacial region
    // we set the ghost velocity
    // real liquid region, set ghost vapor velocity
    if (f[] > 0.5)
    {
      foreach_dimension()
      {
        u2.x[] = u1.x[] - den_diff * mEvapTotE[]*n.x[];
      }
    }
    else
    {
      // real vapor region, set ghost liquid  velocity
      foreach_dimension()
      {
        u1.x[] = u2.x[] + den_diff * mEvapTotE[]*n.x[];
      }
    }
  }

  boundary((scalar *){u1, u2});
}

void setGhostVelFace()
{
  double den_diff = 1. / rho2 - 1. / rho1;
  foreach_face()
  {
#if USE_MY_SOLID
    bool is_solidf = ((int)is_solid[] == 1 && (int)is_solid[-1] == 1);
    if (is_solidf)
      continue;
#endif

    double vof_f = 0.5 * (f[] + f[-1]);
    double ls_f = 0.5 * (ls[] + ls[-1]);
    bool is_inf_m = (int)is_interfacial[-1] == 1;
    bool is_inf_p = (int)is_interfacial[] == 1;
    int num = 0;
    double mdot_f = 0.0;
    double nf = 0.0;
    if(!is_inf_m && !is_inf_p)
    {
      if (ls_f > 0.0)
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
      mdot_f += mEvapTotE[-1];
      nf += n.x[-1];
      num++;
    }
    if(is_inf_p)
    {
      mdot_f += mEvapTotE[];
      nf += n.x[];
      num++;
    }

    if(num == 0)
    {
      printf("EORRO!\n");
      exit(0);
    }
    
    mdot_f /= (double)num;
    nf /= (double)num;

    if (ls_f > 0.0)
    {
      uf2.x[] = uf.x[] - fm.x[] * mdot_f * nf * den_diff;
    }
    else
    {
      uf1.x[] = uf2.x[] + fm.x[] * mdot_f * nf * den_diff;
    }
  }

  boundary((scalar *){uf1, uf2});
}

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction(vector u, face vector uf)
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
We solve the advection equations for $\mathbf{u}_l$ and
$\mathbf{u}_g$. */

event advection_term (i++, last)
{
  if (!stokes) {
    prediction (u1, uf1);
    prediction (u2, uf2);

    setGhostVelFace();
    boundary((scalar *){uf1, uf2, alpha});
    mgpf = project_sf_twofield (uf1, uf2, pf, alpha, dt/2., mgpf.nrelax);
    setGhostVelFace();

    setGhostVelCenter();
    advection ((scalar *){u1}, uf1, dt, (scalar *){g});
    advection ((scalar *){u2}, uf2, dt, (scalar *){g});

#if USE_MY_SOLID
    boundarySolidVelCNoauto(u1);
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
Solving the viscous term we obtain the temporary velocities
$\mathbf{u}_l^*$ and $\mathbf{u}_g^*$. */

event viscous_term (i++, last)
{
  if (constant(mu.x) != 0.) {
    correction (u1, dt);
    correction (u2, dt);
    setGhostVelCenter();
    mgu = viscosity (u1, mu, rho, dt, mgu.nrelax);
    mgu = viscosity (u2, mu, rho, dt, mgu.nrelax);
    correction (u1, -dt);
    correction (u2, -dt);
    setGhostVelCenter();
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
  trash ({uf1,uf2});
  boundary((scalar *){uf1, uf2, a});
  foreach_face() {
    uf1.x[] = fm.x[]*(face_value (u1.x, 0) + dt*a.x[]);
    uf2.x[] = fm.x[]*(face_value (u2.x, 0) + dt*a.x[]);
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
  boundary((scalar *){gf});
  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);

#if USE_MY_SOLID
  foreach ()
  {
    if (is_solid[] == 1.0)
    {
      foreach_dimension()
          g.x[] = 0.0;
    }
  }
#endif
}

/**
## Projection

To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  boundary((scalar *){uf1, uf2, alpha});
  setGhostVelFace();
  mgp = project_sf_twofield (uf1, uf2, p, alpha, dt, mgp.nrelax);
  setGhostVelFace();

  centered_gradient (p, g);

  correction (u1, dt);
  correction (u2, dt);
  setGhostVelCenter();
}

event reconstructions (i++, last) {
  //u.x will not be used
  foreach()
    foreach_dimension()
      u.x[] = u1.x[]*f[] + u2.x[]*(1. - f[]);

  //uf.x shouldn't be set in this way as 
  foreach_face()
    uf.x[] = uf1.x[];

#if USE_MY_SOLID
    boundarySolidVelCNoauto(u);
#endif

    boundary((scalar *){ufext, uf});
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

/**
## References

~~~bib
@article{tanguy2007level,
  title={A level set method for vaporizing two-phase flows},
  author={Tanguy, S{\'e}bastien and M{\'e}nard, Thibaut and Berlemont, Alain},
  journal={Journal of Computational Physics},
  volume={221},
  number={2},
  pages={837--853},
  year={2007},
  publisher={Elsevier}
}

@article{tanguy2014benchmarks,
  title={Benchmarks and numerical methods for the simulation of boiling flows},
  author={Tanguy, S{\'e}bastien and Sagan, Micha{\"e}l and Lalanne, Benjamin and Couderc, Fr{\'e}d{\'e}ric and Colin, Catherine},
  journal={Journal of Computational Physics},
  volume={264},
  pages={1--22},
  year={2014},
  publisher={Elsevier}
}

@article{nguyen2001boundary,
  title={A boundary condition capturing method for incompressible flame discontinuities},
  author={Nguyen, Duc Q and Fedkiw, Ronald P and Kang, Myungjoo},
  journal={Journal of Computational Physics},
  volume={172},
  number={1},
  pages={71--98},
  year={2001},
  publisher={Elsevier}
}
~~~
*/

