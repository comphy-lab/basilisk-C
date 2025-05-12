/**
# Bagnold free surface flow with segregation: Brazil Nuts effect
 
 
 An example of 2D granular flow  of a binary mixture of small and large particules along an inclined plate
 (angle $-\alpha$) with a free surface is presented here.
 The configuration is periodic, what is injected to the left comes from the right.
 On the bottom solid wall there is a no slip boundary condition, on the top fixed free surface a slip condition is imposed (and zero pressure).
 We use the centered Navier-Stokes solver with regularization for viscosity.
 This flow is a simple model for of mixture of pebbles and gravels  avalanche along a slope.

## SeGrayGation
 
 Following Gray and Thornton we incorporate advection due to mean flow 
 ($\frac{\partial c_i}{\partial t}+ u \frac{\partial c_i}{\partial x}+ v\frac{\partial c_i}{\partial y}$),
 percolation-driven segregation $\frac{\partial }{\partial y}(c_i v_{pi})$
 and diffusion due to random particle collisions
 $\frac{\partial }{\partial x}(D \frac{\partial }{\partial x}c_i)+
 \frac{\partial }{\partial y}(D \frac{\partial }{\partial y}c_i)$. This gives a continuum transport equation for the volume concentration
 for a binary mixture  of species i=l,s (i = l and i = s represent large particles and small particles, respectively).
 The percolation velocity may be approximated by
 
 
 $v_{pl}=S_r(\frac{\partial u}{\partial y})(1-c_l)$
 and
 $v_{ps}=-S_r(\frac{\partial u}{\partial y})(1-c_s)$
 where $S_r$ is parameter depending on the granular media.
 $D$ may be approximated by $D_2 \sqrt{2}d_g^2$, here it is just a constant.
 The final equation:
 
 $$\frac{\partial c_i}{\partial t}+ u \frac{\partial c_i}{\partial x}+ v\frac{\partial c_i}{\partial y}+
 S_r \frac{\partial }{\partial y}\left(
 (\frac{\partial u}{\partial y})c_i(1-c_i) \right)=
 \frac{\partial }{\partial x}(D \frac{\partial }{\partial x}c_i)+
 \frac{\partial }{\partial y}(D \frac{\partial }{\partial y}c_i)$$
 
 Boundary condition, equilibrium at the top
 $$  S_r   (\frac{\partial u}{\partial y})|_h c_i(h)(1-c_i(h))  = - D \frac{\partial }{\partial y}c_i|_h$$
 and equilibrium at the wall
 $$  S_r   (\frac{\partial u}{\partial y})|_0 c_i(0)(1-c_i(0))  = - D \frac{\partial }{\partial y}c_i\_0$$
 
 
##  Equations for a granular fluid: $\mu(I)$ rheology
 
 We have to solve the velocity field which carries the mixture.
 By GDR MiDi, in a scalar analysis $\tau = \mu(I ) p$, by Jop et al hypothesis,
 this is generalized :
 tangential stress is linked to the shear rate tensor by
 $$\tau_{ij} =  \sqrt{2}\mu(I) p \frac{D_{ij}}{D_2}$$
 where  the friction law   compatible with the experiments is:
 $$\mu(I)= \mu_0 + \frac{\Delta \mu}{\frac{I_0}{I}+1}$$
 the coefficients depend on the nature of the granular media
 $\mu_0=0.38$ $\Delta \mu = 0.26$ and $I_0=0.3$
 
 
##  Equivalent rheology
 From this, one can exhibit an equivalent viscosity if we write:
 $\tau_{ij} = 2 \mu_{eq}  D_{ij} $
 we the define
 $\mu_{eq}=  \frac{\mu(I) p}{\sqrt{2} D_2 }$
 This equivalent viscosity will be coded next (with a regularisation to avoid infinite viscosity at rest).
 
 
##  Exact solution in the proposed case
 
 This problem admits an analytical solution called "Bagnold" profile in a steady case.
 $$ u = \frac{2}{3} \frac{I_\alpha}{d_g} \sqrt{g h^3 cos(\alpha)}   (1 - (1 - \frac{y}{h})^{3/2})$$
 
 
 
 
 */
/**
 just before, let us include the NS code, define the precision $2^N$, define the angle 24.64 degrees
 */
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
scalar c[],q0[],q[];
scalar * tracers = {c};
mgstats mgf;

#define LEVEL 6
double mumax;
double dg;
double Sr=.5;
# define alpha 0.43
scalar mu_eq[],foo[];
/**
 the exact Bagnold shear velocity is:
 */
double dUb( double y){
    double U0=(sqrt(cos(alpha))*(-0.114 + 0.3*tan(alpha))/(dg*(0.64 - tan(alpha))));
    return ((pow(1. - y, .5))*U0) ;
}
/**
 This gives the Bagnold's solution
 (note the  numerical value  for alpha=0.43, dg=0.04 U0 = 2.06631)
 */
double Ub( double y){
    double U0=(sqrt(cos(alpha))*(-0.114 + 0.3*tan(alpha))/(dg*(0.64 - tan(alpha))));
    return ((1. - pow(1. - y, 1.5))*2./3.*U0) ;
}

double zz( double t){
    return 10 ;
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
    //  alpha = 0.43;
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
    
    /**
     no slip at the bottom
     */
    u.n[bottom] = dirichlet(0);
    uf.n[bottom] = dirichlet(0);
    u.t[bottom] = dirichlet(0);
    /**
     note the pressure
     */
    p[top] = dirichlet(0);
    pf[top] = dirichlet(0);
    /**
     flux note the signs (check!)
     */
    q0[bottom] = neumann(0);
    q0[top] = neumann(0);
    
    c[bottom] = neumann(q0[]);
    c[top] = neumann(-q0[]);
    
    /**
     the $\Delta t_{max}$ should be enough small
     */
    DT = 0.01;
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
     $$-\frac{\partial p}{\partial x} = sin(alpha) $$
     $$-\frac{\partial p}{\partial y} = -cos(alpha) $$
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
        c[] = ( (y<.75) && ( y >.25) ? 1 : 0);
        c[] = ( (y<.75 && x<2./3) ? 1 : 0);
    }
}

/**
 We check the number of iterations of the Poisson and viscous
 problems. */
//event logfile (i++)
// fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
/**
 old value of the velocity and concentration is saved
*/
scalar un[],cn[];
event init_un (i = 0) {
    foreach(){
        un[] = u.x[];
        cn[] = c[];}
}
/**
 so that when it does not more change we are converged, monitor mass conservation as well
 */
event conv (t += 1; i <= 100000) {
    double du = change (u.x, un);
    double dc = change (c, cn);
    double sc = statsf(c).sum;
    
    fprintf(stdout,"t= %g u/Ub=%g mass = %g \n",t,interpolate (u.x, L0/2, .999)/Ub(1),sc);
    if (i > 0 && du < 1e-5 && dc < 1e-5)
        return 1; /* stop */
}
/**
## Implementation of the Bagnold viscosity
*/
event bagnold(i++) {
    /** Compute
     $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
     with $I = d_g \sqrt{2} D_2 /\sqrt{p/\rho}$, and as
     $$\mu(I)= \mu_0 + (\Delta \mu) I/(I_0 + I)$$
     the equivalent viscosity is
     $$\mu_{eq}=   \frac{\mu(I) p}{\sqrt{2} D_2 }$$
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
                foo[] = D2;
                mu_eq[] =  min(muI*fabs(p[])/(sqrt(2.)*D2) , mumax );}
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

## Mixing

The Gray equation is solved here
 
*/
event mixing(i++){
    const face vector D[] = {.01 , .01 };
    scalar r[],logis[],shear[];
/**
  Compute here the "logistic map source term" from percolation in computing a field
    $$  q_0   = - \frac{S_r}{D} (\frac{\partial u}{\partial y})| c_i(1-c_i)$$
 hence at the top and bottom
  $$   \frac{\partial }{\partial y}c_i = -q_0(x,y)$$
 so that at the top
 $$ \frac{\partial }{\partial y}c_i|_h = -q_0(x,h)$$
 so that at the top, as $y=h+n$
 $$ \frac{\partial }{\partial n}c_i|_h = -q_0(x,h)$$
 which was coded before as `c[top] = neumann(-q0[]);`.
 So that at the bottom, as $y=-n$
 $$ \frac{\partial }{\partial n}c_i|_0 = q_0(x,0)$$
      mind the minus sign do to the normal in the BC at the wall !

*/
    foreach() {
        shear[] = (u.x[0,1]-u.x[0,-1])/(2*Delta);
        logis[] = c[]*(1-c[])*(1-.0*c[])*shear[];
        q0[] =  -c[]*(1-c[])*Sr/.01*shear[];
        q[] = (c[0,1]-c[0,-1])/(2*Delta);}
    boundary ({logis,shear,q0,q});
    
/**
 
Compute the source term $r$ for the diffusion equation
$$r=-S_r \frac{\partial }{\partial y}\left( (\frac{\partial u}{\partial y})c_i(1-c_i) \right)$$
     
*/
    foreach() {
        r[] = -Sr*(logis[0,1]-logis[0,-1])/(2*Delta);}
    boundary ({r});
/**

Diffusion equation
     
$$\frac{\partial c_i}{\partial t}+ u \frac{\partial c_i}{\partial x}+ v\frac{\partial c_i}{\partial y}=
     \frac{\partial }{\partial x}(D \frac{\partial }{\partial x}c_i)+
     \frac{\partial }{\partial y}(D \frac{\partial }{\partial y}c_i)+r $$
 
*/
	mgf= diffusion (c, dt, D, r);
}
/**
 Save profiles + film
*/
event movies (i += 10 ) {
    static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
    scalar vorticity[];
    foreach()
    vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
    boundary ({vorticity});
    output_ppm (vorticity, fp, box = {{0,0},{1,1}},
                min = 0, max = 2, linear = true);
    static FILE * fp1 = popen ("ppm2mpeg > c.mpg", "w");
    output_ppm (c, fp1, box = {{0,0},{1,1}},
                linear = true, min = 0, max = 1);
}

event profiles (t += 1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g\n", y, interpolate (u.x, L0/2, y), interpolate (shear, L0/2, y),
                 Ub(y), cos(alpha)*(1-y),interpolate (p, L0/2, y),
                 interpolate (mu_eq, L0/2, y), interpolate (foo, L0/2, y),dUb(y), interpolate (c, L0/2, y),interpolate (q, L0/2, y));
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
 qcc -g -O3 -o bagnold_periodic_segregation bagnold_periodic_segregation.c -lm
 ./bagnold_periodic_segregation
 
 
 lldb bagnold_periodic
~~~
 
 Plots of the velocity,  $\tau$ and $p$ and concentration and flux in the middle of the domain.
 
 
~~~gnuplot  Bagnold flow and concentration of large particules
 set xlabel "y"
 set xlabel "U, tau, p, c"
 set key left
 p'xprof' u 1:2 t'U computed' ,''u 1:($7*$3) t'tau comp.','' u 1:6 t'p',''u 1:10 t'c(L0/2,y) large part.' w lp,''u 1:($11/15) t'flux dc/dy' w lp
~~~
 
 See an animation of the concentration, red corrsponds to the large particules, blue the small
 
![[Animation](bagnold_periodic_segregation/c.mpg) of concentration field.](bagnold_periodic_segregation/c.mpg)
 
 
Gentilly Avril 2015
 
## Related examples
 
* [The basic Bagnold flow](http://basilisk.fr/sandbox/M1EMN/Exemples/bagnold_periodic.c)
 
## Bibliography
  
* Gray Thornton
 ["A theory for particle size segregation in shallow granular free-surface flows"](http://rspa.royalsocietypublishing.org/content/royprsa/461/2057/1447.full.pdf)
 Proc. R. Soc. A (2005) 461, 1447–1473
 
* Lagrée, Staron, Popinet
 ["The granular column collapse as a
 continuum: validity of a two–dimensional
 Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011

* Pierre Jop, Yoël Forterre & Olivier Pouliquen
 "A constitutive law for dense granular flows", Vol 441 8 June 2006 doi:10.1038/nature04801
 
* Fan,  Schlick,   Umbanhowar,  Ottino,  Lueptow
 "Modeling size segregation of granular materials: the roles of segregation, advection, and diffusion"
 Journal of Fluid Mechanics, vol 741, pp 252-279, 2014
 
*/