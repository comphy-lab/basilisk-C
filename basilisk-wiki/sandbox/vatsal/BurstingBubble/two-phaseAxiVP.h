/**
# Two-phase interfacial flows
This is a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). It contains the implementation of
Viscoplastic Fluid (Bingham Fluid).<br/>
This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "vof.h"

scalar f[], * interfaces = {f};
scalar D2[];
face vector D2f[];
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double mumax = 0., tauy = 0.;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */

  if (mu1 || mu2)
    mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(muTemp, mu2, f)  (clamp(f,0.,1.)*(muTemp - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event properties (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif

#if TREE
  sf.prolongation = refine_bilinear;
  boundary ({sf});
#endif

  /**
  This is part where we have made changes.
  $$\mathcal{D}_{11} = \frac{\partial u_r}{\partial r}$$

  $$\mathcal{D}_{22} = \frac{u_r}{r}$$

  $$\mathcal{D}_{13} = \frac{1}{2}\left( \frac{\partial u_r}{\partial z}+ \frac{\partial u_z}{\partial r}\right)$$

  $$\mathcal{D}_{31} = \frac{1}{2}\left( \frac{\partial u_z}{\partial r}+ \frac{\partial u_r}{\partial z}\right)$$

  $$\mathcal{D}_{33} = \frac{\partial u_z}{\partial z}$$

  $$\mathcal{D}_{12} = \mathcal{D}_{23} = 0.$$

  The second invariant is $\mathcal{D}_2=\sqrt{\mathcal{D}_{ij}\mathcal{D}_{ij}}$ (this is the Frobenius norm)
  $$\mathcal{D}_2^2= \mathcal{D}_{ij}\mathcal{D}_{ij}= \mathcal{D}_{11}\mathcal{D}_{11} + \mathcal{D}_{22}\mathcal{D}_{22} + \mathcal{D}_{13}\mathcal{D}_{31} + \mathcal{D}_{31}\mathcal{D}_{13} + \mathcal{D}_{33}\mathcal{D}_{33}$$
  
  **Note:** $\|\mathcal{D}\| = D_2/\sqrt{2}$.<br/>

  We use the formulation as given in [Balmforth et al. (2013)](https://www.annualreviews.org/doi/pdf/10.1146/annurev-fluid-010313-141424), they use $\dot \gamma$ 
  which is by their definition $\sqrt{\frac{1}{2}\dot \gamma_{ij} \dot \gamma_{ij}}$
  and as $\dot \gamma_{ij}=2 D_{ij}$


  Therefore, $\dot \gamma$ = $\sqrt{2} \mathcal{D}_2$, that is why we have a $\sqrt{2}$ in the equations.

  Factorising with $2  \mathcal{D}_{ij}$ to obtain a equivalent viscosity
    $$\tau_{ij} = 2( \mu_0 + \frac{\tau_y}{2 \|\mathcal{D}\| } ) D_{ij}=2( \mu_0 + \frac{\tau_y}{\sqrt{2} \mathcal{D}_2 } ) \mathcal{D}_{ij} $$
  As  defined by [Balmforth et al. (2013)](https://www.annualreviews.org/doi/pdf/10.1146/annurev-fluid-010313-141424)
  $$\tau_{ij} = 2 \mu_{eq}  \mathcal{D}_{ij} $$
  with
  $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} \mathcal{D}_2 }$$

  Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$.

  The fluid flows always, it is not a solid, but a very viscous fluid.
  $$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$

  Reproduced from: [P.-Y. Lagrée's Sandbox](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_simple.c). Here, we use a face implementation of the regularisation method, described [here](http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/Couette_NonNewtonian.c).
  */

  foreach_face(x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;
    double D11 = 0.5*( (u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(2.*Delta) );
    double D22 = (u.y[] + u.y[-1, 0])/(2*max(y, 1e-20));
    double D33 = (u.x[] - u.x[-1,0])/Delta;
    double D13 = 0.5*( (u.y[] - u.y[-1, 0])/Delta + 0.5*( (u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(2.*Delta) ) );

    double D2temp = sqrt( sq(D11) + sq(D22) + sq(D33) + 2*sq(D13) );
    if (D2temp > 0. && tauy > 0.){
      double temp = tauy/(sqrt(2.)*D2temp) + mu1;
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.){
        muTemp = mumax;
      } else {
        muTemp = mu1;
      }
    }
    muv.x[] = fm.x[]*mu(muTemp, mu2, ff);
    D2f.x[] = D2temp;
  }

  foreach_face(y) {
    double ff = (sf[0,0] + sf[0,-1])/2.;
    alphav.y[] = fm.y[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;
    double D11 = (u.y[0,0] - u.y[0,-1])/Delta;
    double D22 = (u.y[0,0] + u.y[0,-1])/(2*max(y, 1e-20));
    double D33 = 0.5*( (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(2.*Delta) );
    double D13 = 0.5*( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*( (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(2.*Delta) ) );

    double D2temp = sqrt( sq(D11) + sq(D22) + sq(D33) + 2*sq(D13) );
    if (D2temp > 0. && tauy > 0.){
      double temp = tauy/(sqrt(2.)*D2temp) + mu1;
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.){
        muTemp = mumax;
      } else {
        muTemp = mu1;
      }
    }
    muv.y[] = fm.y[]*mu(muTemp, mu2, ff);
    D2f.y[] = D2temp;
  }
  /**
  I also calculate a cell-centered scalar D2, where I store $\|\mathbf{\mathcal{D}}\|$. This can also be used for refimnement to accurately refine the fake-yield surfaces.
  */
  foreach(){
    rhov[] = cm[]*rho(sf[]);
    D2[] = (D2f.x[]+D2f.y[]+D2f.x[1,0]+D2f.y[0,1])/4.;
    if (D2[] > 0.){
      D2[] = log(D2[])/log(10);
    } else {
      D2[] = -10;
    }
  }
  boundary(all);
#if TREE
  sf.prolongation = fraction_refine;
  boundary ({sf});
#endif
}
