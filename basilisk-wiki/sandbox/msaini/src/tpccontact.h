/**
# A compressible two-phase flow solver

This solves the compressible Navier-Stokes equations for a
compressible multiphase flow including the total energy equation

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
and an Equation of State (EOS). Here we use the general EOS 
written in the Mie--Gruneisen form

$$
\rho_i e_i = \frac{p_i + \Gamma_i \Pi_i}{\Gamma_i - 1}
$$

The method, described in detail in [Fuster et al., 2018](/src/references.bib#fuster2018), resorts to a
time splitting technique were we first solve the advection step using
the VOF method for $f$, $f_i \rho_i$, $f_i \rho_i \mathbf{u}$ $f_i
\rho_i e_{tot,i}$ and then we perform the projection step using the
all-mach solver.
*/
#include "./all_mach.h"
#include "contact.h"
#include "./old-vof.h"

/** Nomenclature used: frho1 = $f \rho_1$, frho2=$(1-f) \rho_2$, fE1 = $f \rho_1 e_1 + f \rho_1 \mathbf{u}^2$, fE2 = $(1-f) \rho_2 e_2 + (1-f) \rho_2 \mathbf{u}^2$
   */

scalar f[], * interfaces = {f};
scalar frho1[], frho2[], fE1[], fE2[];

face vector alphav[];
scalar rhov[];
(const) face vector lambdav = zerof;

/** Variables of the EOS */
double gamma1 = 1.4, gamma2 = 1.4;
double PI1 = 0., PI2 = 0.;
double mu1 = 0., mu2 = 0.;
double lambdav1 = 0., lambdav2 = 0. ;

#ifndef mu
//#define mu(f) ((f) > (0.5) ? (mu1) : (mu2)) //This is added to better calculate the corrections in energies
#define mu(f) (mu1*mu2/(clamp(f,0.,1.)*(mu2 - mu1) + mu1)) //HM can be used for the averaging of mu at interface
//# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
# define lambdav(f)  (clamp(f,0.,1.)*(lambdav1 - lambdav2 - 0*2./3*(mu1 - mu2)) + lambdav2 - 0*2./3*mu2)
#endif

/**
Auxilliary fields need to be allocated. In particular the quantity $\rho c^2$ */

vector q2[];
scalar rhoc2v[];

/** Naive refinement method */
#if TREE
void conservative_refine (Point point, scalar s)
{  
  double cc = f[];
  double scc = s[];
  if (cc <= 0. || cc >= 1.)
    /* foreach_child() */
    /*   s[] = scc; */
    refine_bilinear(point,s);
    
  else {

    /**
    Otherwise, we reconstruct the interface in the parent cell. */

    coord n = mycs (point, f);
    double alpha = plane_alpha (cc, n);

    /**
    And compute the volume fraction in the quadrant of the coarse cell
    matching the fine cells. We use symmetries to simplify the
    combinations. */

    coord a, b;
    foreach_dimension() {
      a.x = 0.; b.x = 0.5;
    }   
        
    foreach_child() {
      coord nc; 
      foreach_dimension()
        nc.x = child.x*n.x;
      double crefine = rectangle_fraction (nc, alpha, a, b); 
      if (s.inverse)
        s[] = scc/(1. - cc)*(1. - crefine);
      else
        s[] = scc/cc*crefine;
    }   
  }   
}
//Different refine and prolongation functions are defined to respect equation of state
//when refining the cells, The energy is calculated from interpolated pressure to preserve smoothness.
// Similar to f[] fE1[] and fE2[] are pushed to end of list all in event defaults
void fE_refine(Point point, scalar s)
{
 foreach_child(){
   double Ek = 0.;
   foreach_dimension()
     Ek += sq(q.x[]);
   if (s.inverse)
      s[] = (p[] + gamma2*PI2)/(gamma2 - 1.)*(1. - f[])+(0.5*Ek/(frho1[]+frho2[]))*(1.-f[]);
   else
      s[] = (p[] + gamma1*PI1)/(gamma1 - 1.)*f[]+(0.5*Ek/(frho1[]+frho2[]))*f[];

  }
}

#endif

/** We set the default values */



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
    frho1[] = fE1[] = p[] = f[] = 1.;
    frho2[] = fE2[] = 0.;
  }
  boundary ({frho1, frho2, p, f, fE1, fE2});


  /** 
  We also initialize the list of tracers to be advected with the VOF
  function $f$ (or its complementary function) */
  #if TREE

  p.refine = p.prolongation = refine_bilinear;
  for (scalar s in {frho1, frho2, q})
    s.refine = s.prolongation = conservative_refine;
    fE1.refine = fE1.prolongation = fE_refine;
    fE2.refine = fE2.prolongation = fE_refine; // The refine and prolongation functions for fE1 and fE2 are added as above
#endif

    
  f.tracers = list_copy ({frho1, frho2, fE1, fE2, q, q2});
  
  for (scalar s in {frho2, fE2, q2})
    s.inverse = true;


  
}

/** 
For the associated tracers we use the advection scheme provided in
f.gradient */

event init (i = 0) {

  boundary ((scalar *){q});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(q.x[]/(frho1[] + frho2[]) + q.x[-1]/(frho1[-1] + frho2[-1]))/2.;
  boundary ((scalar *){uf});

  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");


  for (scalar s in {frho1, frho2, fE1, fE2, q, q2})
    s.gradient = f.gradient;
  
}

/** 
Before the VOF advection step we obtain the momentum of each phase
separately and then we perform the advection step (vof.h)*/

event vof (i++)
{
  foreach()
    foreach_dimension() {
      double u = q.x[]/(frho1[] + frho2[]);
      q.x[]  = frho1[]*u;
      q2.x[] = frho2[]*u;
    }
  boundary ((scalar *) {q, q2});

}

/** 
After advection we compute the overall averaged momentum */

event tracer_advection (i++)
{
  foreach(){
    foreach_dimension()
      q.x[] += q2.x[];
     if (f[] <= 0.){
      frho1[] = 0.;
      fE1[] = 0.;
    }
    if (f[] >= 1.){
      frho2[] = 0.;
      fE2[] = 0.;
    } 

  }
  boundary ((scalar *){q,frho1,frho2,fE1,fE2});
  //boundary ((scalar *){q});
}

/** 
During the projection step we compute the provisional pressure from
the EOS using the values after advection */

event properties (i++)
{
  foreach() {
    rhov[] = (frho1[] + frho2[]);
   
    double Ek = 0.;
    foreach_dimension()
      Ek += sq(q.x[]);

    double fc = clamp (f[],0.,1.);
    double invgammaavg = fc/(gamma1 - 1.) + (1. - fc)/(gamma2 - 1.);
    double PIGAMMAavg = (fc*PI1*gamma1/(gamma1 - 1.) +
			 (1. - fc)*PI2*gamma2/(gamma2 - 1.));
    ps[] = (fE1[] + fE2[] - Ek/(frho1[] + frho2[])/2. - PIGAMMAavg)/invgammaavg;
    /** We also compute $\rho c^2$ */
    rhoc2v[] = (p[]*(invgammaavg + 1.) + PIGAMMAavg)/invgammaavg;

  }
  boundary ({rhov,frho1,frho2,ps,rhoc2v});
  
  /** If viscosity is present we obtain the averaged viscosity at the cell faces */
  foreach_face() {
    if (mu1 || mu2) {
      face vector muv = mu, lambdavv = lambdav;
      double ff = (f[] + f[-1])/2.;
      muv.x[] = fm.x[]*mu(ff);
      lambdavv.x[] = lambdav(ff);
    }
    /** We also compute the averaged density at the cell faces */
    alphav.x[] = 2.*fm.x[]/(rhov[] + rhov[-1]);
  }

  /* The all-mach solver needs rho*cm */
  foreach() 
    rhov[] *= cm[];
  boundary ((scalar *){alphav,lambdav,mu});
  
  /* I still miss a source term in the divergence related to viscous
   * dissipation. Usually it is small and a little bit complicated to
   * calculate. In the current form of allmach it cannot be added
   * because it does not allow user-defined divergence sources. */
}

/** 
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
/**
First it is computed the divergence of the velocity at the cell centers
using the cell face velocity. Note that the face vector *uf*
incorporates the metric factor. */ 

  scalar divU[];
  foreach () {
    divU[] = 0.;
    foreach_dimension ()
      divU[] += (uf.x[1] - uf.x[])/(Delta*cm[]);
  }
  boundary({divU});

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
S_{xx} & 0 &  0 \\
0 &  S_{yy} &  0 \\
0 & 0 & S_{\theta \theta}  
\end{array}
\right)
$$

with $S = S_{xx} = S_{yy} = S_{\theta \theta}= \lambda_v \nabla \cdot \mathbf{u}$. Hence, the radial
component of $\nabla \cdot \mathbf{S}$ is simply,

$$
(\nabla \cdot \mathbf{S}) \cdot \mathbf{e}_y  = 
\frac{1}{y} \frac{\partial(y S_{yy}) }{\partial y} 
- \frac{S_{\theta \theta}}{y} 
= \frac{\partial S }{\partial y}
$$ 

Also, $u_\theta = 0$.  */ 

  foreach () {
    double momentum = 0.;
    double energy = 0.;
    foreach_dimension () {
      double right = lambdav.x[1]*(divU[1]  + divU[])/2.;
      double left = lambdav.x[] *(divU[-1] + divU[])/2.;
      momentum  += right - left;
      energy += uf.x[1]*right - uf.x[]*left;
    }
    foreach_dimension ()
      q.x[] += dt*momentum/Delta;

    double fc = clamp(f[],0,1);
    fE1[] += fc    *dt*energy/Delta/cm[];
    fE2[] += (1-fc)*dt*energy/Delta/cm[];

  }

  /**
  The velocity $\mathbf{u}$ is obtained from the updated momentum,
  $\mathbf{q}$.*/

  vector u = q;

  foreach() 
    foreach_dimension() 
      u.x[] = q.x[]/(frho1[] + frho2[]);
  boundary((scalar *){u});
  /**
  It is lacking the contribution of the pressure and incompressible
  viscous term to the total energy. Since the face velocity has
  already the metric factor those in $\mu$ has been removed.*/

  face vector taux[];
  // I use alphav to store pressure at faces 
  foreach_face()
      taux.x[] = -(p[]+p[-1])/2.;
  boundary_flux ({taux});
 
  foreach () {
    double energy = 0.; 
    foreach_dimension()
      energy += uf.x[1]*taux.x[1] - uf.x[]*taux.x[];

    double fc = clamp(f[],0,1);
    fE1[] += fc*dt*energy/Delta/cm[];
    fE2[] += (1-fc)*dt*energy/Delta/cm[];
  }   

  //alphav is now the deformation tensor 2e_ijk
  foreach_dimension() {
    foreach_face(x)
      taux.x[] = 2.*(u.x[] - u.x[-1])/Delta;
    #if dimension > 1
      foreach_face(y)
        taux.y[] = (u.x[] - u.x[0,-1] + 
		      (u.y[1,-1] + u.y[1,0])/4. -
		      (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    #endif
    #if dimension > 2
      foreach_face(z)
        taux.z[] = (u.x[] - u.x[0,0,-1] + 
		      (u.z[1,0,-1] + u.z[1,0,0])/4. -
		      (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
    #endif
    boundary_flux ({taux});

    foreach () {
      double energy = 0.; 
      foreach_dimension()
        energy += uf.x[1]*taux.x[1] - uf.x[]*taux.x[];

      double fc = clamp(f[],0,1);
      fE1[] += fc*dt*mu1*energy/Delta/cm[];
      fE2[] += (1-fc)*dt*mu2*energy/Delta/cm[];
    }   
  }


  /**
  Resume momentum and apply boundary conditions.*/

  foreach() 
    foreach_dimension() 
      q.x[] *= (frho1[] + frho2[]);

  boundary ({q, fE1, fE2});


}

/** At the end of the simulation we clean the tracer list */
event cleanup (i = end) {
   free (f.tracers);
   f.tracers = NULL;
}
