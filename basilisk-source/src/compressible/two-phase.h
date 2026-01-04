/**
# A compressible two-phase flow solver

This solves the two-phase compressible Navier-Stokes equations
including the total energy equation.
$$
\frac{\partial (f \rho_i)}{\partial t } + 
\nabla \cdot (f \rho_i \mathbf{u}) = 0 
$$
$$
\frac{\partial (\rho_i \mathbf{u})}{\partial t } + 
\nabla \cdot ( \rho_i  \mathbf{u} \mathbf{u}) = 
-\nabla p + \nabla \cdot \tau_i' 
$$
$$
\frac{\partial [\rho_i (e_i + \mathbf{u}^2/2)]}{\partial t } 
+ \nabla \cdot [ \rho_i \mathbf{u} (e_i  + \mathbf{u}^2/2)] =
-\nabla \cdot (\mathbf{u} p_i) + 
\nabla \cdot \left( \tau'_i \mathbf{u} \right)
$$
an advection equation for the volume fraction $f$
$$
\frac{\partial f}{\partial t} + \mathbf{u} \cdot \nabla f = 0
$$
and an Equation Of State (EOS) that needs to be defined. 

The method is described in detail in [Fuster & Popinet,
2018](#fuster2018) and relies on a time splitting technique were we
first solve the advection step using the VOF method for $f$, $f_i
\rho_i$, $f_i \rho_i \mathbf{u}$ $f_i \rho_i e_{tot,i}$ and then
perform the projection step using the all-Mach solver. */

#include "all-mach.h"

/**
We transport VOF tracers without the one-dimensional compressive term. */

#define NO_1D_COMPRESSION 1
#include "vof.h"

/**
The primary fields are:
$$
\begin{aligned}
\text{frho1} & = f \rho_1, \\
\text{frho2} & = (1-f) \rho_2, \\
\text{q} & = f \rho_1 \mathbf{u} + (1-f) \rho_2 \mathbf{u}, \\
\text{fE1} & = f \rho_1 (e_1 + \mathbf{u}^2/2), \\
\text{fE2} & = (1-f) \rho_2 (e_2 + \mathbf{u}^2/2)
\end{aligned}
$$
*/

scalar f[], * interfaces = {f};
scalar frho1[], frho2[], fE1[], fE2[];

/**
The timestep can be limited by a CFL based on the speed of sound. */

double CFLac = HUGE;

/**
The dynamic viscosities for each phase, as well as the volumetric
viscosity coefficients. */

double mu1 = 0., mu2 = 0.;
double lambdav1 = 0., lambdav2 = 0. ;

/**
These functions are provided by the Equation Of State. See for example
[the Mie--Gruneisen Equation of State](Mie-Gruneisen.h). */

extern double sound_speed          (Point point);
extern double average_pressure     (Point point);
extern double bulk_compressibility (Point point);
extern double internal_energy      (Point point, double fc);

/**
By default the Harmonic mean is used to compute the phase-averaged
dynamic viscosity. */

#ifndef mu
# define mu(f) (mu1*mu2/(clamp(f,0,1)*(mu2 - mu1) + mu1))
#endif

/**
The volumetric viscosity uses arithmetic average by default. */

// fixme: Include the term depending on $\nabla \cdot u$. It is tricky because this
// quantity can be very different upon the phase
#ifndef lambdav
# define lambdav(f)  (clamp(f,0,1)*(lambdav1 - lambdav2) + lambdav2)
// # define lambdav(f)  (clamp(f,0.,1.)*(lambdav1 - lambdav2 - 2./3*(mu1 - mu2)) + lambdav2 - 2./3*mu2)
#endif

/**
## Auxilliary fields

Auxilliary fields need to be allocated. The quantity $\rho c^2$, the
average density `rhov` and its inverse `alphav`. */

scalar rhoc2v[];
face vector alphav[];
scalar rhov[];
const face vector lambdav0[] = {0,0,0};
(const) face vector lambdav = lambdav0;

/**
## Time step restriction based on the speed of sound */

event stability (i++)
{
  if (CFLac < HUGE)
    foreach (reduction (min:dtmax)) {
      double c = sound_speed (point);
      if (CFLac*Delta < c*dtmax)
	dtmax = CFLac*Delta/c;
    }
}

#if TREE
/**
## Energy refinement function

The energy is refined from the refined pressures, momentum and
densities, using the equation of state. */
void fE_refine (Point point, scalar fE)
{
  foreach_child() {
    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);
    Ek /= 2.*(frho1[] + frho2[]);
    fE[] = fE.inverse ?
      (internal_energy (point, 0.) + Ek)*(1. - f[]) :
      (internal_energy (point, 1.) + Ek)*f[];
  }
}
#endif // TREE

/**
## Initialisation

We set the default values. */

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2) {
    mu = new face vector;
    lambdav = new face vector;
  }

  rhoc2 = rhoc2v;  
  foreach () {
    frho1[] = 1., fE1[] = 1.;
    frho2[] = 0., fE2[] = 0.;
    p[] = 1.;
    f[] = 1.;
  }
  
  /** 
  We also initialize the list of tracers to be advected with the VOF
  function $f$ (or its complementary function). */
  
  f.tracers = list_copy ({frho1, frho2, fE1, fE2});
  for (scalar s in {frho2, fE2})
    s.inverse = true;

  /**
  We set limiting. */
  
  for (scalar s in {frho1, frho2, fE1, fE2, q}) {
    s.gradient = minmod2;
#if TREE
    /**
    On trees, we ensure that limiting is also applied to prolongation
    and refinement. */

    s.prolongation = s.refine = refine_linear;
#endif
  }
  
  /**
  We add the interface and the density to the default display. */

  display ("draw_vof (c = 'f');");
  display ("squares (color = 'rhov', spread = -1);");
}

event init (i = 0)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(q.x[]/(frho1[] + frho2[]) + q.x[-1]/(frho1[-1] + frho2[-1]))/2.;

  /**
  We update the fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");

  /** 
  For the associated tracers we use the gradient defined by
  f.gradient. */

  if (f.gradient)
    for (scalar s in {frho1, frho2, fE1, fE2, q})
      s.gradient = f.gradient;  
}

/**
## VOF advection of momentum 

We overload the *vof()* event to transport consistently the volume
fraction and the momentum of each phase. */

static scalar * interfaces1 = NULL;

event vof (i++)
{

  /**
  We split the total momentum $q$ into its two components $q1$ and
  $q2$ associated with $f$ and $1 - f$ respectively. */

  vector q1 = q, q2[];
  foreach()
    foreach_dimension() {
      double u = q.x[]/(frho1[] + frho2[]);
      q1.x[] = frho1[]*u;
      q2.x[] = frho2[]*u;
    }

  /**
  Momentum $q2$ is associated with $1 - f$, so we set the *inverse*
  attribute to *true*. We use the same limiting for q1 and q2. */

  foreach_dimension() {
    q2.x.inverse = true;
    q2.x.gradient = q1.x.gradient;
  }

#if TREE
  /**
  The refinement function is modified by *vof_advection()*. To be able
  to restore it, we store its value. */
  
  void (* refine) (Point, scalar) = q1.x.refine;
#endif
  
  /**
  We associate the transport of $q1$ and $q2$ with $f$ and transport
  all fields consistently using the VOF scheme. */

  scalar * tracers = f.tracers;
  f.tracers = list_concat (tracers, (scalar *){q1, q2});
  vof_advection ({f}, i);
  free (f.tracers);
  f.tracers = tracers;
  
  /**
  We recover the total momentum. */
  
  foreach() {
    foreach_dimension()
      q.x[] = q1.x[] + q2.x[];
    
    /**
    We avoid negative densities and energies which may have been
    caused by round-off during VOF advection. */
    
    if (f[] <= 0.){
      frho1[] = 0.;
      fE1[] = 0.;
    }
    if (f[] >= 1.){
      frho2[] = 0.;
      fE2[] = 0.;
    }
  }
  
#if TREE
  /**
  We restore the refinement function for the total momentum. */

  for (scalar s in {q}) {
    s.refine = s.prolongation = refine;
    s.dirty = true;
  }

#if 0
  /**
  This is switched off by default for now as the standard refinement
  in [/src/vof.h#vof_concentration_refine]() seems to work fine. */
  
  for (scalar s in {fE1,fE2}) {
    s.refine = s.prolongation = fE_refine;
    s.dirty = true;
  }
#endif
#endif // TREE  

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */
  
  interfaces1 = interfaces, interfaces = NULL;  
}

/**
We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
}

/**
## Pressure and density

During the projection step we compute the provisional pressure from
the EOS using the values after advection. */

event properties (i++)
{
  foreach() {
    rhov[] = frho1[] + frho2[];
    ps[] = average_pressure (point);

    /** 
    We also compute $\rho c^2$. */
    
    rhoc2v[] = bulk_compressibility (point);
  }
  
  foreach_face() {
    
    /** 
    If viscosity is present we obtain the averaged viscosity at the cell faces. */
  
    if (mu1 || mu2) {
      face vector muv = mu, lambdavv = lambdav;
      double ff = (f[] + f[-1])/2.;
      muv.x[] = fm.x[]*mu(ff);
      lambdavv.x[] = lambdav(ff);
    }
    
    /** 
    We also compute the averaged density at the cell faces. */
    
    alphav.x[] = 2.*fm.x[]/(rhov[] + rhov[-1]);

  }

  /**
  The all-Mach solver needs rho*cm. */
  
  foreach()
    rhov[] *= cm[];
  
  /* I still miss a source term in the divergence related to viscous
   * dissipation. Usually it is small and a little bit complicated to
   * calculate. In the current form of allmach it cannot be added
   * because it does not allow user-defined divergence sources. */
}

/**
## Update of the total energy

After projection we update the values of the total energy adding the
term missing from the projection step.

For a Newtonian fluid, the stress tensor comprises the pressure term,
the viscous uncompressible term and the viscous compressible term,

$$ \tau = -p \mathbf{I} + \mu (\nabla \mathbf{u} + \nabla \mathbf{u}^T)
+ \lambda_v  (\nabla \cdot \mathbf{u}) \mathbf{I}
$$

where $\mathbf{I}$ is the unity tensor and $\lambda_v = \mu_v - 2/3
\mu$. The volumetric viscosity $\mu_v$ is negligible in most cases. */

event end_timestep (i++)
{
#if 0
  /**
  We first compute the divergence of the velocity at cell centers
  using the cell face velocity. Note that the face vector *uf*
  incorporates the metric factor. */

  scalar divU[];
  foreach () {
    divU[] = 0.;
    foreach_dimension ()
      divU[] += (uf.x[1] - uf.x[])/(Delta*cm[]);
  }

  /**
  We compute explicitly the contribution of the compressible viscous
  term to the momentum equation and the total energy equation,
  $$
  \partial_t \mathbf{q} = \nabla \cdot [\lambda_v(\nabla \cdot \mathbf{u}) \mathbf{I}] 
  $$
  $$
  \partial_t [\rho (e + \mathbf{u}^2/2)] = 
  \nabla \cdot [\lambda_v (\nabla \cdot \mathbf{u}) \mathbf{u}] 
  $$

  The axysimmetric case allows some simplifications. Since
  $$
  \mathbf{S} = \lambda_v(\nabla \cdot \mathbf{u}) \mathbf{I} =
  \left( 
  \begin{array}{ccc}
  S_{xx} & 0 &  0 \				\
  0 &  S_{yy} &  0 \				\
  0 & 0 & S_{\theta \theta}  
  \end{array}
  \right)
  $$
  with $S = S_{xx} = S_{yy} = S_{\theta \theta}= \lambda_v \nabla
  \cdot \mathbf{u}$. Hence, the radial component of $\nabla \cdot
  \mathbf{S}$ is simply
  $$
  (\nabla \cdot \mathbf{S}) \cdot \mathbf{e}_y  = 
  \frac{1}{y} \frac{\partial(y S_{yy}) }{\partial y} 
  - \frac{S_{\theta \theta}}{y} 
  = \frac{\partial S }{\partial y}
  $$
  Also $u_\theta = 0$. */ 

  foreach() {
    double momentum = 0., energy = 0.;
    foreach_dimension() {
      double right = lambdav.x[1]*(divU[1]  + divU[])/2.;
      double left = lambdav.x[]*(divU[-1] + divU[])/2.;
      momentum += right - left;
      energy += uf.x[1]*right - uf.x[]*left;
    }
    foreach_dimension ()
      q.x[] += dt*momentum/Delta;
    energy *= dt/(Delta*cm[]);
    double fc = clamp(f[],0,1);
    fE1[] += energy*fc;
    fE2[] += energy*(1. - fc);
  }
#endif
  
  /**
  The contribution of the pressure to the energy of each phase is
  lacking. */

  {
    face vector upf[];
    foreach_face()
      upf.x[] = - uf.x[]*(p[] + p[-1])/2.;
 
    foreach () {
      double energy = 0.; 
      foreach_dimension()
	energy += upf.x[1] - upf.x[];
      energy *= dt/(Delta*cm[]);
      double fc = clamp(f[],0,1);
      fE1[] += energy*fc;
      fE2[] += energy*(1. - fc);
    }
  }
  
  /**
  This is the contribution of the incompressible viscous term. */

  if (mu1 || mu2) {
    
    /**
    The velocity $\mathbf{u}$ is first obtained from the updated
    momentum $\mathbf{q}$. */
    
    vector u = q;
    foreach()
      foreach_dimension() 
        u.x[] = q.x[]/(frho1[] + frho2[]);

    /**
    We add the incompressible contribution of the $\nabla \cdot (u_i
    \tau'_i)$ term. Note that the compressible contribution, which is
    typically small, is neglected. */

    // fixme: add the compressible contribution (careful, $\nabla \cdot
    // u$ can be very different upon the phase and also very different
    // from the value defined for the mixture with an averaged velocity
    face vector ueijk[];
    foreach_dimension() {
      foreach_face(x)
	ueijk.x[] = 2.*uf.x[]*(u.x[] - u.x[-1])/Delta;
#if dimension > 1
      foreach_face(y)
	ueijk.y[] = uf.y[]*(u.x[] - u.x[0,-1] + 
			    (u.y[1,-1] + u.y[1,0])/4. -
			    (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
#endif
#if dimension > 2
      foreach_face(z)
	ueijk.z[] = uf.z[]*(u.x[] - u.x[0,0,-1] + 
			    (u.z[1,0,-1] + u.z[1,0,0])/4. -
			    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
#endif
      foreach () {
	double energy = 0.; 
	foreach_dimension()
	  energy += ueijk.x[1] - ueijk.x[];
	energy *= dt/(Delta*cm[]);
	double fc = clamp(f[],0,1);
	fE1[] += mu1*energy*fc;
	fE2[] += mu2*energy*(1. - fc);
      }   

      // fixme: Formally, we need to include a correction term due to
      // the normal viscous stress jump (to be done)
    }

    /**
    Finally we recover momentum.*/

    foreach()
      foreach_dimension() 
        q.x[] *= frho1[] + frho2[];
  }
}

/**
At the end of the simulation we clean the tracer list. */

event cleanup (i = end) {
  free (f.tracers);
  f.tracers = NULL;
}

/**
## References

~~~bib
@hal{fuster2018, hal-01845218}
~~~

## See also

* [The Mie--Gruneisen Equation of State](Mie-Gruneisen.h)
*/
