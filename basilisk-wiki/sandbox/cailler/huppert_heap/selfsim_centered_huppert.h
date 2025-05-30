/**
# Self-similar centered Navier--Stokes solver ([1st Huppert problem](http://basilisk.fr/sandbox/M1EMN/Exemples/column_viscous.c))

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

To clearly understand this file, please refer to 
[this dedicated documentation](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README). 
In particular, all the notations employed are defined there.
</div>

We wish to approximate numerically the incompressible 
Navier--Stokes equations into the related self-similar space
of the *Huppert's heap* problem:

$$
\partial_\tau \mathbf{\overline{u}} 
+ \textcolor{Orchid}{
    \boldsymbol{\nabla} \cdot 
    \left(
      \mathbf{\overline{\Lambda}} 
      \otimes \overline{\mathbf{u}} 
    \right)
  }
= 
- \textcolor{purple}{ \overline{\textbf{\textsf{A}}} } \,
  \left( \boldsymbol{\nabla} \overline{p} \right)
+ \dfrac{1}{\mathrm{Re}_1} 
  \left[
    \textcolor{red}{ \overline{\textbf{\textsf{M}}} }
    .
    \boldsymbol{\nabla} 
  \right] \cdot \left(
  2 \, \overline{\textbf{\textsf{S}}}
\right) 
+ \mathbf{\overline{a}}
+ \textcolor{green}{
  \overline{\textbf{\textsf{B}}}
  \,         
  \mathbf{\overline{u}}
  }
- \mathbf{e}_\eta
$$

in the conservative continuous form, or in its discretized formulation:

$$
\dfrac{\mathbf{\overline{u}}^{**} 
- \mathbf{\overline{u}}^n}{\Delta \tau} 
+ \textcolor{Orchid}{\left[
  \boldsymbol{\nabla} \cdot 
  \left(\mathbf{\overline{\Lambda}} 
  \otimes \overline{\mathbf{u}} \right) 
\right]^{n+1/2}
}
= 
\dfrac{1}{\mathrm{Re}_1} 
  \left[
    \textcolor{red}{ \overline{\textbf{\textsf{M}}} }
    .
    \boldsymbol{\nabla} 
  \right] \cdot \left(
  2 \, \overline{\textbf{\textsf{S}}}^{**}
\right) 
+ \textcolor{green}{
  \left[
    \overline{\textbf{\textsf{B}}}  \,         
    \mathbf{\overline{u}}
  \right]^{n+1/2}
}
+ [ \mathbf{\overline{a}} - \mathbf{e}_\eta ]^{n+1/2}
$$
$$
\left.\begin{array}{rcl}
  \boldsymbol{\nabla} \cdot \left( 
    \textcolor{purple}{ \overline{\textbf{\textsf{A}}} } \,
    \boldsymbol{\nabla} \overline{p}^{n+1}
  \right)
  &=& \dfrac{\boldsymbol{\nabla} \cdot \mathbf{\overline{u}}^{**}}{\Delta \tau} \\
  \mathbf{\overline{u}}^{n+1} 
  &=& \mathbf{\overline{u}}^{**} 
  - \Delta \tau 
    \textcolor{purple}{ \overline{\textbf{\textsf{A}}} } \, 
    \boldsymbol{\nabla} \overline{p}^{n+1}
\end{array}
\right\} 
$$

where:

$$
\left\{ 
  \begin{array}{rcl}
  \tau &=& \ln t \\
  \textcolor{Orchid}{\mathbf{\overline{\Lambda}}} 
  &\textcolor{Orchid}{=}& 
  \textcolor{Orchid}{\overline{\mathbf{u}} + 
  \, \boldsymbol{\xi}}  \\
  \boldsymbol{\xi} &=& {}^{\mathrm{t}}(- \xi, \, \eta) / 5
    = {}^{\mathrm{t}}(-x.t^{-1/5}, \, y.t^{1/5}) / 5 \\ 
  \textcolor{red}{ \overline{\textbf{\textsf{M}}} }
  &\textcolor{red}{=}&
  \textcolor{red}{\begin{pmatrix}
    \mathrm{e}^{3\tau/5} & 0 \\
    0 & \mathrm{e}^{7\tau/5}
  \end{pmatrix}}
  \\
  \overline{\textbf{\textsf{S}}} 
  &=& \left( 
    \boldsymbol{\nabla} \overline{\mathbf{u}} 
    + {}^{\mathrm{t}} \boldsymbol{\nabla} \overline{\mathbf{u}} 
  \right)/2 \\
  \textcolor{green}{ \overline{\textbf{\textsf{B}}} } 
  &\textcolor{green}{=}& 
  \textcolor{green}{\begin{pmatrix}
    4/5 & 0 \\
    0 & 6/5
  \end{pmatrix}} 
  \\
  \textcolor{purple}{ \overline{\textbf{\textsf{A}}} }
  &\textcolor{purple}{=}&
  \textcolor{purple}{\begin{pmatrix}
    \mathrm{e}^{7\tau/5} & 0 \\
    0 & \mathrm{e}^{11\tau/5}
  \end{pmatrix}}
  \end{array}
\right.
$$


The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, a *slightly modified*
Bell-Colella-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
#include "timestep.h"
#include "selfsim_bcg_huppert.h"
#include "viscosity.h"


/**
The primary variables are the centered self-similar 
pressure field $\overline{p}$, the
centered self-similar advected velocity field $\mathbf{\overline{u}}$, 
the centered *position vector* $\boldsymbol{\xi}$. 
The centered self-similar vector field
$\mathbf{\overline{g}}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face *advection velocity* field 
$\mathbf{\overline{\Lambda}}_f$, an auxilliary face *advected velocity* field 
$\mathbf{\overline{u}}_f$, as well as an associated 
centered pressure field $\overline{p}_f$ and a *temporary* 
auxilliary face advected velocity field $\mathbf{\overline{u}}_f^{temp}$. */

scalar p[];
vector u[], g[], xi[];
scalar pf[];
face vector uf[], uf_temp[], lambdaf[];

/**
This problem is a special case where 
$\boldsymbol{\nabla} \cdot \overline{ \mathbf{\Lambda} } = 0$, therefore 
we can simply take into account this property by setting to zero 
the parameter for the number of dimensions, in order to let 
unchanged the files for self-similar solver, 
compared to previous problems 
(*Keller \& Miskis*, *Lamb--Oseen*): */
double Nd = 0.;



/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{\overline{a}}$ defines additional acceleration terms; 
default is zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\boldsymbol{\nabla}\cdot(\mathbf{\overline{\Lambda}}\otimes\mathbf{\overline{u}})$ 
is omitted. This is a reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

#define NEW_VERSION_CENTERED 1

#if NEW_VERSION_CENTERED
  /* New version of centered.h */
  (const) face vector mu = zerof, a = zerof, alpha = unityf;
  (const) scalar rho = unity;
  mgstats mgp = {0}, mgpf = {0}, mgu = {0};
  bool stokes = false;
#else
  /* Old version of centered.h */
  (const) face vector mu[] = zerof; 
  (const) face vector alpha[] = unityf;

  (const) face vector a = zerof;
  (const) scalar rho = unity;
  mgstats mgp, mgpf, mgu;
  bool stokes = false;
#endif


/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{\overline{a}}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{\overline{a}}$, this can be written */

#define neumann_pressure(i) (a.n[i]*fm.n[i]/(alpha.n[i]))

p[right] = neumann (neumann_pressure(ghost) 
  + (fm.n[ghost]/alpha.n[ghost])*((-.2)*x*(u.n[1]-u.n[])/Delta)) ;
  
p[left]  = neumann (- neumann_pressure(0)
  - (fm.n[0]/alpha.n[0])*((-.2)*x*(u.n[]-u.n[-1])/Delta));


xi.n[right] = -.2*x ; 
xi.t[right] = .2*y ; 
xi.n[left] = -.2*x ;
xi.t[left] = .2*y ; 

lambdaf.n[right] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[right] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ; 
lambdaf.n[left] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[left] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ;

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost)
  + (fm.n[ghost]/alpha.n[ghost])*((.2)*y*(u.n[0,1]-u.n[])/Delta));

xi.n[top] = .2*y ;
xi.t[top] = -.2*x ; 
xi.n[bottom] = .2*y ;
xi.t[bottom] = -.2*x ; 

lambdaf.n[top] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[top] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ;
lambdaf.n[bottom] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[bottom] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ; 

#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost)
  + (.2)*y*(u.n[0,1]-u.n[])/(Delta*alpha.n[ghost]));
p[bottom] = neumann (- neumann_pressure(0) 
  + (.2)*y*(u.n[]-u.n[0,-1])/(Delta*alpha.n[0])) ;

xi.n[top] = .2*y ;
xi.t[top] = -.2*x ; 
xi.n[bottom] = .2*y ;
xi.t[bottom] = -.2*x ; 

lambdaf.n[top] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[top] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ;
lambdaf.n[bottom] = uf.n[ghost] + fm.n[ghost]*xi.n[ghost] ;
lambdaf.t[bottom] = uf.t[ghost] + fm.t[ghost]*xi.t[ghost] ; 

#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI


/**
## Initial conditions */

event defaults (i = 0)
{

  /**
  We reset the multigrid parameters to their default values. */
  
  #if NEW_VERSION_CENTERED
    /* New version of centered.h */
    mgp = (mgstats){0};
    mgpf = (mgstats){0};
    mgu = (mgstats){0};  
  #endif

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric and the
  solid factors). */

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
#endif // TREE
}


/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocities and fluid
properties. */

double dtmax;

event init (i = 0)
{
  trash ({uf});
  trash ({lambdaf});

  // defining position vector ξ
  foreach() {
    xi.x[] = -.2*x;
    xi.y[] = .2*y;
  }
  boundary({xi});

  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);

  boundary((scalar *){uf});
  foreach_face()
    lambdaf.x[] = uf.x[] + fm.x[]*face_value (xi.x, 0);
  boundary((scalar *){lambdaf});

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
applied to the advecting face centered velocity field 
$\mathbf{\overline{\Lambda}}_f$; and the timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : timestep (lambdaf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $\tau+\Delta \tau/2$) here. Note that this assumes that tracer fields
are defined at time $\tau-\Delta \tau/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$\tau+\Delta \tau/2$ -- can be defined by overloading this event. */

event properties (i++, last);


/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\boldsymbol{\nabla}\cdot(\mathbf{\overline{\Lambda}}\otimes\mathbf{\overline{u}})$, 
we need to define the face velocity field $\mathbf{\overline{u}}_f$ 
at time $\tau+\Delta \tau/2$. 
We use a **modified** version of the Bell-Colella-Glaz 
[advection scheme](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_bcg_keller.h) 
that is adapted for the advection velocity \mathbf{\overline{\Lambda}}:

$$
\mathbf{\overline{u}}_{f,p}^{n+1/2} 
    = 
    \mathbf{BCG} \left( 
      \mathbf{\overline{u}}^n, \,
      \textcolor{Orchid}{\left.\overline{\lambda}_d^n\right|_{c \to f}}, \,
      \textcolor{Orchid}{\overline{\lambda}_{\perp d}^n}, \, 
      \left.\overline{\mathbf{g}}^n\right|_{c \to f} 
    \right)
$$

We also use the pressure gradient and acceleration terms at time $\tau$ 
(stored in vector $\mathbf{\overline{g}}$). */

void prediction_uf()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
	      du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
	      du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }

  trash ({uf});

  foreach_face() {
    double un = (dt/Delta)*((u.x[] + u.x[-1]) + (xi.x[] + xi.x[-1]))/2.;
    double s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + ((g.x[] + g.x[-1]))*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.; 
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = (u.y[i] + xi.y[i]) < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= (dt/2.)*(u.y[i] + xi.y[i])*fyy/Delta; 
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

We predict the face velocity field $\mathbf{\overline{u}}_f$ at time $\tau+\Delta
\tau/2$ then project it to make it divergence-free. This allows us to 
build the advection face velocity field $\mathbf{\overline{\Lambda}}_f$.  */

event advection_term (i++,last)
{
  if (!stokes) {
    prediction_uf();
    boundary({xi});
    boundary((scalar *){uf});

    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    boundary((scalar *){uf});

    foreach_face(){
      lambdaf.x[] = uf.x[] + fm.x[]*face_value (xi.x, 0);
    }

/**
We store the value of the face velocity field $\mathbf{\overline{u}}_f$ 
at time $\tau+\Delta \tau/2$ in a temporary face vector field, as it will be 
needed later for the treatment of the source term along with 
the acceleration terms:*/

    trash({uf_temp});
    foreach_face(){
      uf_temp.x[] = uf.x[];
    }

/** 
We can now use the face velocity field $\mathbf{\overline{u}}_f$ to
compute the velocity advection term, using the **modified** standard
Bell-Colella-Glaz advection scheme for each component of the velocity
field, with the advecting velocity $\mathbf{\overline{\Lambda}}_f$:

$$
\textcolor{blue}{\overline{\mathbf{U}}_{f}^{n+1/2}}
=
\mathbf{BCG} \left( 
    \mathbf{\overline{u}}^n, \,
    \textcolor{Orchid}{\overline{\lambda}_{f,d}^{n+1/2}}, \,
    \textcolor{Orchid}{\left.\overline{\lambda}_{f, \perp d}^{n+1/2}\right|_{f \to c}}, \, 
    \left.\overline{\mathbf{g}}^n\right|_{c \to f} 
\right)
$$
where $\left.{}\right|_{f \to c}$ describes a linear interpolation from *faces* 
towards *cell centers* (and vice-versa).
*/

    boundary({xi});
    boundary((scalar *){lambdaf});
    boundary((scalar *){uf_temp});

    selfsim_advection ((scalar *){u}, lambdaf, dt, Nd, (scalar *){g});
    boundary((scalar *){lambdaf});
    boundary((scalar *){u});
    boundary({xi});
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
      u.x[] += dt*(g.x[]);
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $\tau$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $\tau+\Delta \tau$. */

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

The acceleration term $\mathbf{\overline{a}}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{\overline{g}}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $\tau+\Delta \tau$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added, **along with the source term defined with 
the temporary face velocity field at time $\tau + \Delta \tau/2$**. */

event acceleration (i++,last)
{
  trash ({uf});
  foreach_face(x){
      uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]) 
        + dt*(.8)*uf_temp.x[];
  }
  foreach_face(y){
      uf.y[] = fm.y[]*(face_value (u.y, 0) + dt*a.y[]) 
        + dt*((1.2)*uf_temp.y[] - fm.y[]*exp(2.2*t));
  }
  boundary((scalar *){uf});
}

/**
## Approximate projection

This function constructs the self-similar centered pressure gradient and
acceleration field $\overline{g}$ using the self-similar face-centered 
acceleration field $\overline{a}$
and the self-similar cell-centered pressure field $\overline{p}$. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{\overline{g}}_f$ combining both
  acceleration and pressure gradient, **and including the self-similar 
  source term defined at $\tau + \Delta \tau/2$** (*cf.* green term in 
  the Navier--Stokes equation at the beginning of this page). */

  face vector gf[];
  foreach_face(x){
      gf.x[] = fm.x[]*(a.x[]) + (.8)*uf_temp.x[]
        - alpha.x[]*(p[] - p[-1])/Delta ;
  }
  foreach_face(y){
      gf.y[] = fm.y[]*(a.y[]) + (1.2)*uf_temp.y[] - fm.y[]*exp(2.2*t)
        - alpha.y[]*(p[] - p[-1])/Delta ;
  }

  /**
  We average these face values to obtain the centered, combined
  acceleration, source term and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

/**
To get the pressure field at time $\tau + \Delta \tau$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field $\overline{g}$. 
We ensure that $\mathbf{\overline{\Lambda}}_f$ is also up-to-date.*/

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);

  trash ({lambdaf});
  boundary ((scalar *){uf});
  boundary({xi});
  foreach_face(){
    lambdaf.x[] = uf.x[] + fm.x[]*face_value (xi.x, 0);
  }
  boundary ((scalar *){lambdaf});

  centered_gradient (p, g);

  /**
  We add the gradient field $\overline{g}$ to the centered velocity field. */

  correction (dt);
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
  event ("properties");
}
#endif

