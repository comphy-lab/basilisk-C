/**
 
# Granular rheology, $\mu(I)$

## Purpose
This is the implementation of $\mu(I)$ rheology (GDR MiDi, Jop Pouliquen Forterre).
 
 
 
## Implementation
 We use
*/
#include "navier-stokes/centered.h"
/**
 with two phases,
*/
#include "vof.h"
/**
 One phase is a
 passive fluid of small density  to preserve 0 pressure at the surface of the other phase, the granular one. `f` is the `VOF`  field.
 
 possible values
 `define RHOF 1e-4`
 `define mug  1e-5`, here we take:
*/
#define RHOF 1e-3
#define mug  1e-4

face vector alphav[];
face vector muv[];
face vector eta[];
scalar rhov[];
// note the v
scalar f[];
scalar * interfaces = {f};
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
double Co=0,D=1./30,mus=0.32,dmu=.28,I0=.4;
double etamax=10000;
/**
Do not forget
 Boundary conditions for granular flow, pressure must be zero at the surface:
`p[top] = dirichlet(-RHOF*LDOMAIN);`
 
 
Total density
 */
#define rho(f) ((f) + RHOF*(1. - (f)))
/**
 beware that  $\alpha = 1/\rho$ and the $\mu$ is the dynamic viscosity or absolute viscosity.
 `mu` is a reserved variable of `basilisk`.
 To avoid confusion we use symbol $\eta$ in papers and in the redefinition
 of dynamic viscosity.
*/
event init_granul (t = 0) {
    alpha = alphav;
    mu = muv;
    rho = rhov;
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
 
 

 
*Note*
 
 
 In the pure shear flow
 $D_{11}=D_{22}=0$ et $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
 so that
 $D_2=\sqrt{D_{ij}D_{ij}} =\sqrt{ 2 D_{12}^2} = \frac{\partial u}{ \sqrt{2}\partial y}$.
 In a pure shear flow, $\partial u/\partial y= \sqrt{2} D_2$.
 The inertial number $I$ is indeed $D \partial u/\partial y/\sqrt(p)$

 */

event properties (i++) {
    
    trash ({alphav});
    scalar eta[];
    
    scalar fa[];
    foreach()
    fa[] = (4.*f[] +
            2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
            f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary ({fa});
    
    foreach() {
        eta[] = mug;
        if (p[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            if (D2 > 0.) {
                D2 = sqrt(2.*D2)/(2.*Delta); // this D2 is sqrt(2) D2
                double In = D2*D/sqrt(p[]);
                double muI = mus + dmu*In/(I0 + In);
                double etamin = sqrt(D*D*D);
                eta[] = max((Co + muI*p[])/D2, etamin);// this D2 is sqrt(2) D2
                eta[] = min(eta[],etamax);      }
        }
    }
    boundary ({eta});

    foreach_face() {
        double fm = (fa[] + fa[-1,0])/2.;
        muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
        // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug); was in Gerris
        alphav.x[] = 1./rho(fm);
    }
    foreach()
    rhov[] = rho(fa[]);
    boundary ({muv,alphav,rhov});
}

/**

## alternate possibility
 
  [http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/LidDrivenBingham.c]()
 
 
 event properties (i++) {
 
 scalar fa[];
 foreach()
 fa[] = (4.*f[] +
 2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
 f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
 foreach()
 fa[]= (2*f[] + (f[-1,0] + f[1,0] + f[0,-1] + f[0,1]))/6;
 boundary ({fa});
 
 trash ({alphav});
 foreach_face() {
 eta.x[] = mug;
 if (p[] > 0.) {
 double D2 = 0.;
 foreach_dimension() {
 double dxx = u.x[1,0] - u.x[-1,0];
 double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
 D2 += sq(dxx) + sq(dxy);
 }
 if (D2 > 0.) {
 D2 = sqrt(2.*D2)/(2.*Delta); // this D2 is sqrt(2) D2
 double In = D2*D/sqrt(p[]);
 double muI = mus + dmu*In/(I0 + In);
 double etamin = sqrt(D*D*D);
 eta.x[] = max((Co + muI*p[])/D2, etamin);// this D2 is sqrt(2) D2
 eta.x[] = min(eta.x[],100000);      }
 }
 
 double fm = (fa[] + fa[-1,0])/2.;
 muv.x[] = (fm*eta.x[] + (1. - fm)*mug);
 alphav.x[] = 1./rho(fm);
 }
 boundary ({eta});
 
 foreach()
 rhov[] = rho(fa[]);
 boundary ({muv,alphav,rhov});
 }
 
 
# Biblio
 
 * Jop
 
 * GDR MiDi

 
# Exemples of use
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column.c]()  granular collapse
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_cohesif.c]()  cohesive granular collapse
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass.c]() granular sandglass
 

 
# Exemples of implemetations not yet in `granular.h`, some en axi
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_segregation.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnoldtube.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/forchheimer.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/siloaxi.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/siloaxisl_2.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/siloaxiForchheimer.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_column_muw.c]()
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/granular_front.c]()
 
 * [http://basilisk.fr/sandbox/zhenhai/LateralSandglassAngle40.c]()
 
 * [http://basilisk.fr/sandbox/yixian/README]()
 
 Confiné 04/2020
 
 
 
*/
