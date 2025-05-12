/**
#   Free surface flow of a Bingham fluid with Navier Stokes solver

An example of 2D complex flow  over a plate with a free surface is presented here.
The configuration is periodic what is injected to the left comes from the right. 
On the bottom solid wall there is a no slip boundary condition, on the top fixed free surface a slip condition is imposed.
We use the centered Navier-Stokes solver with regularization for viscosity.

 
 
 
##  Equations for a Bingham fluid 
 
 The Bingham fluid is a non Newtonian fluid.
 The Bingham viscosity that we will put in the Navier-Stokes solver is such that if
$||\tau|| \le  \tau_y$ then there is no motion $D=0$
 if the stress is high enough
 $||\tau|| >  \tau_y$
 note that $||\tau||$ is the modulus defined as the Euclidian norm  $\sqrt{\frac{1}{2}{\tau_{ij} \tau_{ij}}}$.

It is not $\sqrt{\tau_{11}^2 + \tau_{12}^2}$ as in Balmorth 06, which is the Frobenius norm.  

 
 
 then stress tensor is linked to the shear rate tensor by
 $$\tau_{ij} = 2 \mu_0  D_{ij}  + \tau_y \frac{D_{ij}}{||D||}$$
 
 where $D_{ij}$ is the shear strain rate tensor
 (tenseur de taux de déformation)  $D_{ij}=(u_{i,j}+u_{j,i})/2$, 
 the components in 2D:
 
 $D_{11}=\frac{\partial u}{\partial x}$,
 $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{22}=\frac{\partial v}{\partial y}$
 
 
 in the Euclidian norm we have:
  $$||D||=\sqrt{\frac{D_{ij}D_{ij}}{2}}$$ 
 The second invariant defined by $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
 hence
  $$D_2^2= D_{ij}D_{ij}= ( \frac{\partial u}{\partial x})^2 + (\frac{\partial v}{\partial y})^2 +  \frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})^2$$
and we have obviously
$||D|| = D_2/\sqrt{2}$ 
 
 
 We present here the formulation in Balmforth, he uses $\dot \gamma$ 
 which is by his definition $\sqrt{\frac{1}{2}\dot \gamma_{ij} \dot \gamma_{ij}}$
 and as $\dot \gamma_{ij}=2 D_{ij}$
 then $\dot \gamma$ 
 is $\sqrt{2} D_2$, that is why we have a $\sqrt{2}$ in the equations.
 Factorising with $2  D_{ij}$ to obtain a equivalent viscosity
  $$\tau_{ij} = 2( \mu_0 + \frac{\tau_y}{2 ||D|| } ) D_{ij}=2( \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 } ) D_{ij} $$
as  defined by Balmforth
 $$\tau_{ij} = 2 \mu_{eq}  D_{ij} $$
with
 $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
 
 
 
 
 
 
##  Exact solution in the proposed case
 
 
 We look at an unidirectional flow, a pure shear flow  $u(y)$, $v=0$, so
 $D_{11}=D_{22}=0$ and $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
this gives  (mind square root of 2):
 $D_2=\sqrt{D_{ij}D_{ij}} = \frac{1}{\sqrt{2}}\frac{\partial u}{\partial y}$.
 
 
$$ \tau_{12} = 2 \mu  D_{12}  + 2 \tau_y \frac{D_{12}}{ \sqrt{2}D_2} =
   \mu    \frac{\partial u}{\partial y} + \tau_y $$
 
 Equilibrium between pressure gradient and viscosity (writting $\tau$ for a shorthand of $\tau_{12}$)
 $$0=-\frac{\partial p}{\partial x} + \frac{\partial \tau}{\partial y}$$
 as there is no stress at the free surface $y=h$, the stress is
 $$ \tau = (-\frac{\partial p}{\partial x})(h-y)$$
 the stress $\tau$ increases from the free surface, as long as $\tau<\tau_y$,
we are under the threshold,
 so shear is zero: $\frac{\partial u}{\partial y} =0$,  hence velocity is constant, say it is $U$.
 Let us define
  $Y=h-\tau_y/(-\frac{\partial p}{\partial x})$, where  $\tau=\tau_y$.
 
 
 So :
 $Y<y<h$, $\tau<\tau_y$, $\frac{\partial u}{\partial y} =0$, and $u=U$
 
 
 Then going down:
 $0<y<Y$ we have $\tau =\tau_y +\mu \frac{\partial u}{\partial y}$. 
 
 
 This gives:
 $$\tau_y +\mu \frac{\partial u}{\partial y} =  (-\frac{\partial p}{\partial x})(h-y)
 $$
and this allows to solve for the velocity profile
 $$u = (-\frac{\partial p}{\partial x})(hy-\frac{1}{2}y^2) - \tau_y y$$
which is indeed zero in $y=0$, and for  $y=Y$, we have the plug flow $u=U$ of value:
 $$U= \frac{(\tau_y - h(-\frac{\partial p}{\partial x}) )^2}{2 \mu  (-\frac{\partial p}{\partial x})}$$
 
 as in $y=0$ $u=0$ and in  $y=Y$, $\tau=\tau_y$, so  $\mu  \frac{\partial u}{\partial y} =0$ and is zero after, therefore
 $$u =U  \left(1- \left( \frac{y-Y}{Y} \right)^2\right)\;   with  \;U= \frac{(\tau_y - h(-\frac{\partial p}{\partial x}) )^2}{2 \mu  (-\frac{\partial p}{\partial x})} \;    with  \; Y = h-\tau_y/(-\frac{\partial p}{\partial x})$$
 
 flux $\int_0^h u dy=\int_0^Y u dy + U(h-Y)$  is
 $$\int_0^h u dy= \frac{hU}{3}(1-Y/h)$$
 
 
 
 */
 

#include "navier-stokes/centered.h"
#define LEVEL 6
double tau_y,mu_0,mumax;
/**
The domain is one unit long. $0<x<1$ $0<y<1$
*/

int main() {
  L0 = 1.;
  origin (0., 0);
/**
  Values of yeld stress and viscosity
*/
  tau_y=0.25;
  mu_0 = 1;
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
    u.n[top] = neumann(0);
    uf.n[top] = neumann(0);
/**
 no slip at the bottom
*/
    u.n[bottom] = dirichlet(0);
    uf.n[bottom] = dirichlet(0);
    u.t[bottom] = dirichlet(0);
/**
 note the pressure
 */
    p[top] = neumann(0);
    pf[top] = neumann(0);
    p[bottom] = neumann(0);
    pf[bottom] = neumann(0);
/** 
 the $\Delta t_{max}$ should be enough small
*/
  DT = 0.05;
  run(); 
}

/**

*/

face vector muv[];
scalar D2c[];

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
  presure gradient `mdpdx`
 $$-\frac{dp}{dx} = 1 $$
*/
  const face vector mdpdx[] = {1.0,0.};
  a = mdpdx;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
}

/**
We check the number of iterations of the Poisson and viscous
problems. */
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
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
event conv (t += 0.1; i <= 10000) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g \n",t);
    if (i > 0 && du < 1e-5)
        return 1; /* stop */
}
/**
## Implementation of the Bingham viscosity
*/
event bingham(i++) {
    scalar muB[];
/** Compute
 $D_{11}=\frac{\partial u}{\partial x}$,
 $D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{21} =D_{12} =\frac{1}{2}( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x})$,
 $D_{22}=\frac{\partial v}{\partial y}$
 
 And where the second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
 hence
 $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
the equivalent viscosity is
  $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
 the final one is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
*/
    foreach() {
        double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            if (D2 > 0.) {
                D2 = sqrt(D2)/(2.*Delta);
                D2c[] = D2;
                muB[] =  min(tau_y/(sqrt(2.)*D2) + mu_0 , mumax );    }
            }
    boundary ({muB});
    foreach_face() {
        muv.x[] = (muB[] + muB[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}
/**
  Save profiles
*/

event profiles (t += .1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g\n", y, interpolate (u.x, L0/2, y), interpolate (shear, L0/2, y),
                  interpolate (D2c, L0/2, y));
    fclose (fp);
}
/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
//  adapt_wavelet ({u}, (double[]){3e-2,3e-2}, 7, 4);
}
/**
## Results and plots
 
To run the program
 
~~~bash
 qcc -g -O3 -o bingham_simple bingham_simple.c -lm
 ./bingham_simple
 
 
 lldb bingham_simple
~~~

Plot of the velocity and shear
 
 
~~~gnuplot Velocity and stress profiles for Bingham flow
 set xlabel "y"
 set ylabel "u, shear"
 Y = 1-.25
 U = Y*Y/2/1
 p[:][:]'xprof' u 1:2 t'U computed' w lp ,''u 1:3 t'tau comp.',(x<Y ?U*(1-((x-Y)/Y)**2):U) t 'U exact',(x<Y?Y-x:0) t'shear exact'
~~~



Paris Avril 2015



## Bibliography

* [related example in Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/couette.html)

* [related example with augmeted Lagrangian](http://basilisk.fr/sandbox/popinet/poiseuille-periodic.c)

* [Application to the 1D Collapse](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)
 
*  K. F. Liu and C. C. Mei
 "Approximate equations for the slow spreading of a thin sheet of Bingham plastic fluid"
  Phys. Fluids A 2, 30 (1990); doi: 10.1063/1.857821
 
*  K. F. Liu and C. C. Mei
 Slow spreading of a sheet of Bingham fluid on an inclined plane
Fluid Mech. (1989), vol. 207. p p . 505-529
 
 
* Balmforth  Craster Rust  and Sassi
 Viscoplastic flow over an inclined surface
 J. Non-Newtonian Fluid Mech. 139 (2006) 103–127
 
 
* [see also for face implementation](http://basilisk.fr/sandbox/vatsal/GenaralizedNewtonian/Couette_NonNewtonian.c) 
 
 
*/
