/**
# Incompressible Navier--Stokes solver (centered formulation)
The developed module is based on the Immersed Boundary Method (IBM), namely **Brinkman Penalization Method** (BPM) [Liu & Vasilyev 2007](http://www.skoltech.ru/app/data/uploads/sites/19/2017/02/JCP_2007.pdf).
Briefly, the idea is next: let us consider a viscous incompressible flow around a bunch of porous obstacles $O_i$. Usually the boundary condition on the surface is no-slip, therefore, the velocity of the flow on the surface is equal to the velocity of the solid $\mathbf{U}_t$.

$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) =
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] +
\overbrace{\mathbf{a}}^{\text{volume force}} - \overbrace{\frac{f_s}{\eta_s}\left(\mathbf{u}-\mathbf{U}_t\right)}^{\text{penalization term}}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$, $\eta_s$ is a penalization coefficient, characterizes ability of solids to permeate. The less $\eta_s$, the more precise.

In order to distinguish solids from fluids we will use mask array $f_s$:
$$
f_s (\mathbf{x})=
\begin{cases}
1, \text{in solids}\\
0, \text{out of solids}\\
\end{cases}
$$

Let's consider the last term in detail. If a cell is not in an obstacle, then this term is equal to 0, but if a cell belongs to solids, then the term  $-\frac{f_s}{\eta_s}\left(\mathbf{u}-\mathbf{U}_t\right)$ is dominant (we need to guarantee it choosing appropriate penalization coefficient $\eta_s$)
$$
|RHS|\ll\frac{|\mathbf{u}-\mathbf{U}_t|}{\eta_s}
$$

In this case the system of Navier--Stokes equations is reduced to 
$$
\partial_t\mathbf{u} = -\frac{\mathbf{u}-\mathbf{U}_t}{\eta_s}
$$
which has next analytical solution of exponential convergence of velocity $\mathbf{u}$ to preset value $\mathbf{U}_t$ with characteristic time $\eta_s$.
$$
\mathbf{u} = \mathbf{U}_t + \left(\mathbf{u}_0-\mathbf{U}_t \right)\exp\left( -t/\eta_s\right)
$$

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver.
There are 4 ways of solid treatment: 

(i) mask (no MPI, multiphase flow)

(ii) embedded boundaries method (EBM) (with MPI, one phase flow, very accurate)

(iii) Popinet's trick (with MPI, multiphase flow, inaccurate)

(iv) Brinkman Penalization method (BPM) (with MPI, multiphase flow, accurate, moving solids)

More detailed information in [compatability table](http://basilisk.fr/src/COMPATIBILITY). As you can see Brinkman Penalisation method is appropriate for multiphase flows.

In the Brinkman Penalization method the penalization term is handled implicitly together with the viscous term. Fortunately, this linear term lead to diagonally dominant "matrix", which will converges quickly.
Detailed information about verification you can find [here](http://basilisk.fr/sandbox/weugene/cylinder_penalization.c).

To activate Brinkman Penalization Method it is necessary to%

1) define BRINKMAN_PENALIZATION

2) need to set solid mask field: scalar fs[] and velocity of solid vector Us[]

3) Note that penalization coefficient $\eta_s$ default value is 1e-15
 */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#ifndef BRINKMAN_PENALIZATION
#include "viscosity.h"
#else
#include "viscosity-weugene.h"
#endif
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

(const) face vector mu = zerof, a = zerof, alpha = unityf, kappa = zerof;
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
    boundary ((scalar *){alpha});
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
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE
}

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  boundary ((scalar *){uf});

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

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last) {
  boundary ({alpha, mu, rho, kappa}); // Weugene: kappa added
}

event tracer_diffusion (i++,last);

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
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
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
	  du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }
  boundary ((scalar *){du});

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
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
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
time $t+\Delta t$. Note that if BRINKMAN_PENALIZATION variable is defined then the penalization term is implemented into the **new viscosity function**.*/

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    correction (-dt);
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
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);
  boundary ((scalar *){uf, a});
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
  boundary_flux ({gf});

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
  boundary ((scalar *){g});
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (dt);
}
/**
  The projection operation can spoil the boundary conditions on velocity, therefore it is necessary to correct it each time step.*/
#if BRINKMAN_PENALIZATION
event brinkman_penalization(i++, last){
    brinkman_correction(u, uf, rho, dt);
}
#endif
/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);


/**
Output vtk files*/
event vtk_file (i++, last);// Added by Weugene
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
  boundary ((scalar *){uf});
#endif
  event ("properties");
}
#endif

/**
## Useful utilities
These methods can be moved to utility.h...

void MinMaxValues calculates threshold $\epsilon$ based on min and max value. This utility can be helpful with  adapt_wavelet function. If EPS_MAXA=1, then $\epsilon_A=\epsilon_{A0} \max A$, otherwise $\epsilon_A=\epsilon_{A0} \frac{\min A +\max A}{2}$.
*/

void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
  double arr[10][2];
  int ilist = 0;
  for (scalar s in list) {
    double mina= HUGE, maxa= -HUGE;
    foreach( reduction(min:mina) reduction(max:maxa) ){
      if (fabs(s[]) < mina) mina = fabs(s[]);
      if (fabs(s[]) > maxa) maxa = fabs(s[]);
    }
    arr[ilist][0] = mina;
    arr[ilist][1] = maxa;
    ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
  }

  for (int i = 0; i < ilist; i++){
#if EPS_MAXA == 1
    arr_eps[i] *=arr[i][1];
#else
    arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
#endif
#ifdef DEBUG_MINMAXVALUES
    fprintf(stderr, "MinMaxValues: i=%d, min=%g, max=%g, eps=%g\n", i, arr[i][0], arr[i][1], arr_eps[i]);
#endif
  }
}

/**
## Useful utilities for Brinkman penalization method
Redefined functions which neglect values in solids (if fs[] =1).
*/

stats statsf_weugene (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] < 1. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

stats statsf_weugene2 (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] == 0. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

norm normf_weugene (scalar f, scalar fs)
{
  double dvr, avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg)
  reduction(+:rms) reduction(+:volume))
  if (fs[] < 1. && f[] != nodata) {
    dvr = dv()*(1. - fs[]);
    double v = fabs(f[]);
    if (v > max) max = v;
    volume += dvr;
    avg    += dvr*v;
    rms    += dvr*sq(v);
  }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

double change_weugene (scalar s, scalar sn, scalar fs)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    if (fs[] < 1) {
      double ds = fabs (s[] - sn[]);
      if (ds > max)
        max = ds;
    }
    sn[] = s[];
  }
  return max;
}

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/