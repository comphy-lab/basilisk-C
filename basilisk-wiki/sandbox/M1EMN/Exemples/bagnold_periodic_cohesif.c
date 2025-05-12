/**
#  Periodic free surface flow of a Plastic fluid

##  Problem
 
An example of 2D cohesive granular flow  along an inclined plate (angle $-\alpha$) with a free surface is presented here.
The configuration is periodic what is injected to the left comes from the right. 
On the bottom solid wall there is a no slip boundary condition, on the top fixed free surface a slip condition is imposed (and zero pressure).
We use the centered Navier-Stokes solver with regularization for viscosity.
 This flow is a simple model for rock/gravels/sand wet avalanches along a slope.

 It is a kind of Nusselt film (or half Poiseuille) non newtonian counterpart.
 
##  Equations for a granular cohesive fluid: $\mu(I)$ rheology
 
 
 The $\mu(I)$ cohesive rheology supposes that the tangential and normal stress are
 proportional (GDR MiDi):
$$\tau = \tau_c + \mu(I) P  $$
with a coefficient of proportionality that is a function of a single dimensionless number, called the inertial number $I$ (square root of Savage Number...), and $\tau_c$ the yield stress:
 $$I = \frac{d_g\dot \gamma}{\sqrt{p/\rho}}$$
 where $d_g$ is the grain diameter and $\dot\gamma$ reduces to $\partial u/\partial y$ for an  unidirectional flow (note that the final analytical result and numercial computations we present are a for a linearized case $\mu=\mu_0+ a   \frac{d \frac{\partial u}{\partial y}}{\sqrt{p/\rho}}$ with $a=\Delta \mu/I_0$).
 
 
 
 To write this in tensorial form, $\tau_{ij}$ is linked with
 $D_{ij}$  the shear strain rate tensor
 (tenseur de taux de déformation)  $D_{ij}=(u_{i,j}+u_{j,i})/2$,
whose components in 2D are:
 
 $$D_{11}=\frac{\partial u}{\partial x}, \text{ }
 D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}), \text{ }
 D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}),\text{ and }D_{22}=\frac{\partial v}{\partial y}$$
 
 And where the second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
 hence
 $$D_2^2= D_{ij}D_{ij}= ( \frac{\partial u}{\partial x})^2 + (\frac{\partial v}{\partial y})^2 +  \frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})^2$$

for a unidirectional flow  $\partial u/\partial y$ is $\sqrt{2} D_2$
 
The total stress is written in a pressure part and a dissipative part
  $$\sigma_{ij}= - p \delta_{ij} + \tau_{ij}$$
 
 tangential stress is linked to the shear rate tensor by
  $$\tau_{ij} =  \sqrt{2}(\tau_c +\mu(I) p)\frac{D_{ij}}{D_2}$$
 where  the friction law   compatible with the experiments and linearised for smalll $I$  is:
 $$\mu(I)= \mu_0 + \frac{\Delta \mu}{ I_0} I$$
 the coefficients depend on the nature of the granular media
 $\mu_0=0.38$ $\Delta \mu = 0.26$ and $I_0=0.3$
 
 
##  Equivalent rheology
 From this, one can exhibit an equivalent viscosity if we write:
 $$\tau_{ij} = 2 \mu_{eq}  D_{ij} $$
we thus define
 $$\mu_{eq}=  \frac{\tau_c  +\mu(I) p}{\sqrt{2} D_2 }$$
 This equivalent viscosity will be coded next (with a regularisation to avoid infinite viscosity at rest).

 
##  Exact solution in the proposed case
 
 This problem admits an analytical solution called "Bagnold" profile in a steady case with no cohesion.  We introduce the cohesion yield stress and obtain a kind of mofified Bingham fluid.
 This is a steady flow invariant in $x$.
 We look at an unidirectional flow, a pure shear flow  $u(y)$, $v=0$, so
 $D_{11}=D_{22}=0$ and $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
this gives  (mind square root of 2):
 $D_2=\sqrt{D_{ij}D_{ij}} = \frac{1}{\sqrt{2}}\frac{\partial u}{\partial y}$.
 
 so that we have, as we wish:
$$ \tau_{12} =   2 (\tau_c  +\mu(I) p) \frac{D_{12}}{ \sqrt{2}D_2} =
  \tau_c + \mu(I)p$$
 The component  $\tau_{12}$ is as we wanted $\tau_c+\mu(I)p$.
 
 
Equilibrium between pressure gradient and viscosity (writing $\tau$ for a shorthand of $\tau_{12}$), the projection of the equations
$$0=\rho g \sin(\alpha) + \frac{\partial \tau}{\partial y}$$
$$0=-\rho g \sin(\alpha) -\frac{\partial p}{\partial y} $$
as there is no stress and no pressure at the free surface $y=h$, the stress is
$$ \tau = \rho g \sin(\alpha) (h-y)$$
 and the pressure is
$$p = \rho g \cos(\alpha) (h-y)$$
 the stress $\tau$  and the pressure increase from the free surface where they are both 0.
From those expresssions, the ratio of the $\tau$ and $p$ remains constant: $\tau / p =\tan(\alpha)$, but  by the Coulomb fricton $\tau$ is $\mu(I) p + \tau_c$, hence:
$$\mu(I) + \tau_c/p =  \tan(\alpha)$$
For sake of simplification we use here a linearised $\mu(I)$ with $a=\Delta \mu/I_0$:
 $$\mu(I)=\mu_0+ a   \frac{d \frac{\partial u}{\partial y}}{\sqrt{p/\rho}},
 \text{  with equilibrium  }\mu(I) =  \tan(\alpha)-\tau_c/p$$
 hence the shear is (as $I = d_g \frac{\partial u}{\partial y } /\sqrt{p/\rho}$ and $p=\rho g \cos(\alpha) (h-y)$) :
 $$\frac{\partial u}{\partial y }= (\frac{(\tan \alpha    - \mu_0)\sqrt{g \cos \alpha}}{a d })\sqrt{ (h-y)}
 - \frac{\tau_c \sqrt{g \cos \alpha}}{a \rho g d \cos \alpha}  \frac{1}{\sqrt{ (h-y)}}$$
 it is zero
 $\frac{\partial u}{\partial y }=0$  for $y>Y$, the yield height defined by
 $$Y= h -\frac{\tau_c}{\rho g \cos \alpha (\tan \alpha    - \mu_0)}.$$
 We decompose the shear in two contributions, the Bagnold $\frac{\partial u}{\partial y }\bigg\rvert_B$ one and the cohesive one:
 $\frac{\partial u}{\partial y }\bigg\rvert_c$
 as follows
 $$\frac{\partial u}{\partial y }=     \frac{\partial u}{\partial y }\bigg\rvert_B +  \frac{\partial u}{\partial y }\bigg\rvert_c$$
 with $K_B =(\frac{(\tan \alpha    - \mu_0)\sqrt{g \cos \alpha}}{a d })$ and $t_c=\frac{\tau_c }{\rho g h}  ( \frac{\sqrt{g }}{d})\frac{1}{a \sqrt{\cos \alpha}}$
 
 The sum :  $\frac{\partial u}{\partial y }= K_B \sqrt{ (h-y)}
 - t_c \frac{1}{\sqrt{ (h-y)}}$, with
 $\frac{\partial u}{\partial y }=0$  for $y>Y$ with $Y= h -\frac{t_c}{K_B}$
 is integrated, for $y<Y$ in:
 $$u=\frac{1}{3} \left(2 h^{3/2} K+2 \sqrt{h-y} (K_B y+3 t_c)-2 h K_B \sqrt{h-y}-6 \sqrt{h} t_c\right)$$
 or with $t_c=K_B(h-Y)$
 $$u=\frac{2}{3} \left(h^{3/2} K_B-h K_B \sqrt{h-y}+K_B y \sqrt{h-y}+3 t_c \sqrt{h-y}-3 \sqrt{h}
 t_c\right)$$
 with the plug flow velocity for $y>Y$:
 $$U=\frac{2}{3} \left(h^{3/2} K_B-3 \sqrt{h} t_c+2 t_c \sqrt{\frac{t}{K_B}}\right)$$

# Code
 
 just before, let us include the NS code, define the precision $2^N$, define the angle 28.64 degrees
*/
#include "navier-stokes/centered.h"
#define LEVEL 5
double mumax;
double dg;
# define alpha 0.6
# define tauc 0.125
    scalar mu_eq[],foo[];
/**
now, come back to the velocity shear which  is
  $$\frac{\partial u}{\partial y }= K_B \sqrt{ (h-y)}
 - t_c \frac{1}{\sqrt{ (h-y)}}$$
 with  $\frac{\partial u}{\partial y }=0$  for $y>Y$ with $Y= h -\frac{t_c}{K}$
 $K_B =(\frac{(\tan \alpha    - \mu_0)\sqrt{g \cos \alpha}}{a d })$ and $t_c=
 - \frac{\tau_c }{\rho g h}  ( \frac{\sqrt{g }}{d})\frac{1}{a \sqrt{\cos \alpha}}$
 */
double dUbc( double y){
    double dU;
    double a = 0.26/.3;
    double h = 1;
    double KB = ((tan(alpha) - 0.38)* sqrt(cos(alpha))/(a*dg));
    double tc = tauc/(sqrt(cos(alpha))*(a*dg));
    double Y = h -tc/KB;
    dU=(y<Y? KB*sqrt(h - y) - tc/sqrt(h - y):0);
    return dU;
}
/**
This gives the Bagnold's solution if $\tau_c=0$
$$ u = \frac{2}{3} \frac{I_\alpha}{d_g} \sqrt{g h^3 \cos(\alpha)}   (1 - (1 - \frac{y}{h})^{3/2})$$
 (note the  numerical value  for alpha=0.43, dg=0.04 U0 = 2.06631)
*/
double Ub( double y){
    double U0=(sqrt(cos(alpha))*(-0.114 + 0.3*tan(alpha))/(dg*(0.64 - tan(alpha))));
    return ((1. - pow(1. - y, 1.5))*2./3.*U0) ;
}
/**
 This gives the Cohesive Bagnold's solution if $\tau_c \ne0$
 $$ u =  \frac{(2(h^{3/2})K_B-3\sqrt{h} t_c-h K_B \sqrt{(h-y)}+3 t_c\sqrt{h-y}+K_B \sqrt{h-y}y)))}{3}$$
 
 $$U =\frac{2}{3} ((h^{3/2} K_B - 3 t_c \sqrt(h)+ 2 t_c \sqrt{(t_c/K_B)}))$$
 */
double Ubc( double y){
    double a = 0.26/.3;
    double h = 1;
    double KB = ((tan(alpha) - 0.38)* sqrt(cos(alpha))/(a*dg));
    double tc = tauc/(sqrt(cos(alpha))*(a*dg));
    double Y = h -tc/KB;
    double U ,Uy;
    Uy=(2*(pow(h,1.5)*KB-3*sqrt(h)*tc-h*KB*sqrt(h-y)+3*tc*sqrt(h-y)+KB*sqrt(h-y)*y))/3.;
    U = 2*((pow(h,1.5)*KB - 3*tc *sqrt(h)+ 2*tc*sqrt(tc/KB)))/3.;
    U = (y<Y? Uy : U);
    return (U) ;
}
/**
The domain is one unit long. $0<x<1$ $0<y<1$
*/

int main() {
  L0 = 1.;
  origin (0., 0);
/**
  Values of grain size
*/
  dg = 0.04;
/**
  the regularisation value of viscosity 
*/
  mumax=1000;
/**
 Boundary conditions are periodic
*/
    periodic (right);
/**
  slip at the top
*/
    u.t[top] = neumann(0);
  //  uf.t[top] = neumann(0);
    u.n[top] = dirichlet(0);
 //   uf.n[top] = neumann(0);
    /* u.t[top] = neumann(0);
    u.n[top] = neumann(0);
    uf.n[top] = neumann(0);*/
/**
 no slip at the bottom
*/
    u.n[bottom] = dirichlet(0);
   // uf.n[bottom] = dirichlet(0);
    u.t[bottom] = dirichlet(0);
/**
 note the pressure
 */
    p[top] = dirichlet(0);
 //   pf[top] = dirichlet(0);
 // with NO BC it works!
   // p[bottom] = neumann(0);
   // pf[bottom] = neumann(0);
   // p[bottom] = neumann(cos(alpha));
   // pf[bottom] = neumann(cos(alpha));
/**
 the $\Delta t_{max}$ should be enough small
*/
  DT = 0.75;
/** because $U=u(y) e_x$. `stokes true` => no CFL condition
*/
  stokes = true;
  run(); 
}
/**

*/
face vector muv[];

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
  pressure gradient, gravity acceleartion `mdpdx`
 $$-\frac{\partial p}{\partial x} = \sin(alpha) $$
 $$-\frac{\partial p}{\partial y} = -\cos(alpha) $$
*/
    const face vector mdpdx[] = {sin(alpha),-cos(alpha)};

  a = mdpdx;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = Ub(y);
    u.y[] = 0;
    p[]=cos(alpha)*(1-y);
  }
}

/**
We check the number of iterations of the Poisson and viscous
problems. */
//event logfile (i++)
// fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
/**
 old value of the velocity is saved
*/
scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}
/**
 so that when it does not more change we are converged
*/
event conv (t += 1; i <= 100000) {
    double du = change (u.x, un);
    fprintf(stdout,"#time Ucal/Utheo Ucalc err \n %g %g %g %g \n",t,interpolate (u.x, L0/2, .999)/Ubc(1),interpolate (u.x, L0/2, .999),du);
    if (i > 0 && du < .50e-4) //was -6
        return 1; /* stop */
}
/**
## Implementation of the Bagnold cohesive viscosity
*/
event bagnold(i++) {
    /** Compute
     $D_{11}=\frac{\partial u}{\partial x}$,
     $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
     $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
     $D_{22}=\frac{\partial v}{\partial y}$
     
     And where the second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
     hence
     $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
     
     with $I = d_g \sqrt{2} D_2 /\sqrt{p/\rho}$, and as
     $$\mu(I)= \mu_0 + (\Delta \mu) I/(I_0 + I)$$
     
     
     the equivalent viscosity is
     $$\mu_{eq}=   \frac{\tau_c+\mu(I) p}{\sqrt{2} D_2 }$$
     the final one is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
     */
    foreach() {
        mu_eq[] =    mumax;
        if (p[] > 0.) {
        double D2 = 0.,In,muI;
        foreach_dimension() {
            double dxx = u.x[1,0] - u.x[-1,0];
            double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
            D2 += sq(dxx) + sq(dxy);
        }
        if (D2 > 0.) {
            D2 = sqrt(D2)/(2.*Delta) ;
            In = sqrt(2.)*dg*D2/sqrt(fabs(p[]));
            muI = .38 + (.26)*In/(.3 + In);
            muI = .38 + (.26)*In/(.3 );
            foo[] = D2;
            mu_eq[] =  min( (tauc + muI*fabs(p[]))/(sqrt(2.)*D2) , mumax );}
        else{
         mu_eq[] =  mumax;
          }
    }
    }
    boundary ({mu_eq});
    foreach_face() {
        muv.x[] = (mu_eq[] + mu_eq[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}
/**
  Save profiles
*/
event profiles (t += 1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g %g %g %g %g %g \n", y, interpolate (u.x, L0/2, y), interpolate (shear, L0/2, y),
                  Ubc(y), cos(alpha)*(1-y),interpolate (p, L0/2, y),
                 interpolate (mu_eq, L0/2, y), interpolate (foo, L0/2, y),dUbc(y) );
    fclose (fp);
}
/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
 // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}
/**
## Results and plots
 
To run the program
 
~~~bash
 qcc -g -O3 -o bagnold_periodic_cohesif bagnold_periodic_cohesif.c -lm
 ./bagnold_periodic_cohesif
 
 lldb bagnold_periodic
 
 make bagnold_periodic_cohesif.tst;make bagnold_periodic_cohesif/plots;
 make bagnold_periodic_cohesif.c.html;open bagnold_periodic_cohesif.c.html;
 
~~~

 bug from MacOS
~~~bash
 cp bagnold_periodic_cohesif.c.html XXX.c.html
 
 sed -i - 's/\\)//g' XXX.c.html;sed -i - 's/\\(//g' XXX.c.html
 sed -i - 's/\\\[//g' XXX.c.html
 sed -i - 's/\\\]//g' XXX.c.html;
 mv XXX.c.html bagnold_periodic_cohesif.c.html
 open bagnold_periodic_cohesif.c.html
 
~~~

Plots of the velocity,  $\tau$ and $p$, the two last are linear as expected:

 
~~~gnuplot Velocity, pressure and $\tau= \mu_{eq} du/dy$ profiles computed
 set xlabel "y"
 set xlabel "U, tau, p"
 p'xprof' u 1:2 t'U computed' ,''u 1:($7*$3) t'tau comp.','' u 1:6 t'p'
~~~

We check the pressure $p$ is $cos(alpha)*(1-y)$
 
 
~~~gnuplot pressure profiles compared to the lithostatic for Cohesive Bagnold flow
 set xlabel "y"
 set xlabel "p"
 p[][0:]'xprof' u 1:6 t'Pression',''u 1:($5) t'cos(alpha)*(1-y)' w l
~~~
 
 we verify that the computed velocity is the Bagnold one.
 
 
 
~~~gnuplot Velocity and stress profiles computed and exact for Bagnold flow
 set xlabel "y"
 p[][0:2.5]'xprof' u 1:($2*1.) t'U computed',''u 1:4 t'U Cohesive Bagnold exact' w l linec 1,'' u 1:($6) t'pressure' w l
~~~

We check that $(\partial u/ \partial y)$ is $\sqrt{2}D_2$
 and that the computed gradient of velocity is the Bagnold one.
 
 
~~~gnuplot shear velocity profiles for Cohesive Bagnold flow
 set xlabel "y"
  p[][0:]'xprof' u 1:3 t'Computed dUbc/dy',''u 1:($8*sqrt(2))t'Sqrt(2)*D2' ,''u 1:($9) t' d Ubc/dy'w l linec 1
~~~
 
 
## Related examples, links
 

 * [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic_cohesif.c]()
 
 * [related example in Basilisk for pure Bagnold](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
 * [related example in Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/poiseuille.html#bagnold)
 
 
 
## Bibliography
 
 * Lagrée, Staron, Popinet
 ["The granular column collapse as a
 continuum: validity of a two–dimensional
 Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011
 
 * Pierre Jop, Yoël Forterre & Olivier Pouliquen
 "A constitutive law for dense granular flows", Vol 441 8 June 2006 doi:10.1038/nature04801


 

 
 Paris  dec 2019
 
*/

