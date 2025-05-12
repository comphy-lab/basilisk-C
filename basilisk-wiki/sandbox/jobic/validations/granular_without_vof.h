/**
 
# Granular rheology, $\mu(I)$

## Purpose
This is the implementation of $\mu(I)$ rheology (GDR MiDi, Jop Pouliquen Forterre).
 
 
 
## Implementation
 We use
*/
#include "navier-stokes/centered.h"

#define mug  1e-5

double Co=0.,dg=0.04,mus=0.383,dmu=0.26,I0=0.279,lambdaR=1e-6;
scalar eta[],foo[];

/**
  the regularisation value of viscosity 
*/
double mumax=10000;

face vector muv[];

/**
 
## parameters of the rheology
 
 by definition  (By GDR MiDi) from a scalar analysis of a unidirectional shear flow $\overrightarrow{u}=u(y)\overrightarrow{e}_x$,
 the tangential and normal stress are linked through a Coulomb relation :
 $\tau = \mu(I ) p$, the friction law is:
 $$\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}\text{ and where } I=D \frac{\partial u/\partial y}{\sqrt(p)}$$
 with in the definition of $I$,
the diameter of grains  $D=1./30$
(30 diameters of grains in the unit), and $\partial u/\partial y$ the shear of velocity and $p$
 pressure supposed to be zero at the free surface.
 The friction coefficients are
 $\mu_s=0.32$, $\Delta \mu=.28$, and $I_0=.4$, the cohesion is $C_o=0$, it can be changed.
 

In this file, we code the tensorial version $\tau_{ij}$, the shear will be replaced by $\sqrt{2} \sqrt{D_{ij}D_{ji}}$ where   $D_{ij}$  is the shear strain rate tensor
 (see for analytical solutions [the basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 and  [Bagnold with cohesion](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c)

All those numerical values are to be changed in the main.
 
*/

/**
 beware that  $\alpha = 1/\rho$ and the $\mu$ is the dynamic viscosity or absolute viscosity.
 `mu` is a reserved variable of `basilisk`.
 To avoid confusion we use symbol $\eta$ in papers and in the redefinition
 of dynamic viscosity.
*/
event init_granul (t = 0) {
    mu = muv;
}
/**
## computing the viscosity
 
 
 The total stress is written in a pressure part and a dissipative part
 $$\sigma_{ij}= - p \delta_{ij} + \tau_{ij}$$
 
 To write this in tensorial form, $\tau_{ij}$ is linked with
 $D_{ij}$  the shear strain rate tensor
 (*tenseur des taux de déformation*)  $D_{ij}=(u_{i,j}+u_{j,i})/2$,
 
 $$ \tau_{ij}  =2 \eta D_{ij}$$
 
 where $\eta$ is the equivalent dynamic viscosity
 
 
The components of $D_{ij}$ in 2D are:
 $D_{11}=\frac{\partial u}{\partial x}$,
 $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{22}=\frac{\partial v}{\partial y}$
 
 And  the second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
 hence
 $$D_2^2= D_{ij}D_{ij}= ( \frac{\partial u}{\partial x})^2 + (\frac{\partial v}{\partial y})^2 +  \frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})^2$$
 
 
 The $\mu(I)$ is a non newtonain viscosity, it looks like the Bingham one, but with a threshold depending on the pressure.
 
 The viscosity is computed as a fonction of  $D_2=\sqrt{D_{ij}D_{ji}}$ and $p$;
 
  The inertial number $I=D \sqrt{2} D_2/\sqrt(p)$
 and $\mu(I) = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$, the equivalent dynamic viscosity is $$\eta = \frac{\mu(I)p}{\sqrt{2} D_2}$$
 by Jop et al hypothesis,
 tangential stress is linked to the shear rate tensor by
 $$\tau_{ij} = 2 \mu(I) p \frac{D_{ij}}{ \sqrt{2} D_2}$$
 

 This equivalent viscosity will be coded next.
 If the shear is too small, we have to do a regularisation to avoid infinite viscosity at rest.
 so $\eta=min(\eta,\eta_{max})$.
 
Note that if $\eta$ is too small an artificial small viscosity $\rho D \sqrt{gD}$
 is taken see Lagrée et al. 11 § 2.3.1
 
 
 
 */
 
#define regul 1
event properties (i++) {
  foreach() {
    #if regul  
      double valRegul;
    #endif
    eta[] =  mumax;
    if (p[] > 0.) {
      double D2 = 0.,In,muI;
      foreach_dimension() {
        double dxx = u.x[1,0] - u.x[-1,0];
        double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
        D2 += sq(dxx) + sq(dxy);
      }
      #if regul
        D2 = sqrt(2.)*sqrt(D2)/(2.*Delta) ;
	foo[] = D2;
	valRegul=sqrt(D2*D2+lambdaR*lambdaR);
	In = D2/sqrt(fabs(p[]));
	muI = mus + (dmu)*In/(I0 + In);
	eta[] =  (Co+muI*p[])/valRegul;
      #else
        if (D2 > 0.) {
          D2 = sqrt(2.)*sqrt(D2)/(2.*Delta) ;
	  foo[] = D2;
          In = D2/sqrt(fabs(p[]));
          muI = mus + dmu*In/(I0 + In);
          eta[] =  min( (Co + muI*fabs(p[]))/D2 , mumax );
        } else{
          eta[] =  mumax;
        }
      #endif
    }
  }
  boundary ({eta});
  foreach_face() {
      muv.x[] = (eta[] + eta[-1,0])/2.;
  }
  boundary ({muv});
}

/**
 
# Biblio
 
 * Jop
 
 * GDR MiDi

*/
