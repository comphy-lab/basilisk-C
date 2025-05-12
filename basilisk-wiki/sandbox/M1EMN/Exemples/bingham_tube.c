/** 

#  Non Newtonian Bingham flow in a pipe 

Flow in a pipe of radius $R$,
in axi:
$$\tau = \tau_0 + 2 \mu_0   \frac{\partial u}{\partial r}$$
$\mu_0$ viscosity, $\tau_0$ yield stress,  
written 
$$\tau = 2 \eta_{eq} D$$
The momentum equation reduces to 
$$0 =-\partial_x p  + \frac{\partial }{r \partial r}r \tau$$
Pressure gradient $-\partial_x p =-2 \frac{\tau_w}{R}$

definition of Bingham number $Bi = \frac{\tau_0}{\mu_0 U_0/R}$


*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define LEVEL 5
double mumax,Bi;
scalar mu_eq[];

/** 
Bingham solution for comparison, 

for $r>R \tau_0/\tau_w$ 
$$ 
u = \frac{\tau_wR}{2 \mu_0}[1- (\frac{r}{R})^2]-  \frac{\tau_0R}{ \mu_0}[1- (\frac{r}{R})]
$$
for $r<R \tau_0/\tau_w$ 
$$ 
u = \frac{\tau_wR}{2 \mu_0}[1- \frac{\tau_0}{\tau_w}]^2
$$

Written without dimension, with $U_0= \frac{2 \tau_wR}{ \mu_0}$ and $r=R y$, so
 $Bi= \frac{\tau_0 R}{\mu_0 U_0}=  \frac{\tau_0}{2 \tau_w}$.

 So if $y> 2 Bi$ the velocity without dimension is:
$$ 
u = \frac{1}{4}(1- y^2)-  Bi(1- y)
$$

*/
double Ubi( double y, double bi){
     return (y < 2*bi? (1-4*bi*bi)/4-bi*(1-2*bi) :(1-y*y)/4-bi*(1-y)); }
/**
Main with parameters
*/
 int main() {
  L0 = 1.;  
  Bi = 0.2;
/**
  the regularisation value of viscosity 
*/
  mumax=1000;
/**
 Boundary conditions are periodic
*/
    periodic (right);
/**
  no slip at the top
*/
    u.t[top] = dirichlet(0);
    u.n[top] = dirichlet(0);
/**
  symmetry at the bottom
*/  
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = neumann(0); 
 // stokes = true; // because U=u(y)e_x. stokes true > no CFL condition
  run(); 
}


face vector muv[];

event init (t = 0) {
/** 
 prepare viscosity
*/
  mu = muv;
/**
  minus pressure gradient, or equivalent gravity acceleration `mdpdx` 
 $$-\frac{\partial p}{\partial x} = 1 $$
 $$-\frac{\partial p}{\partial y} = 0 $$
*/
    const face vector mdpdx[] = {1,0};
/**
 note that in "accceleration" in "navier-stokes/centered.h" there is the `fm`metric term in front.

 `event acceleration (i++,last)`

  `uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);`

  will be the same for `alphav.x[] = fm.x[]/rho(ff); ` next...
*/
  a = mdpdx;
/**
 Initialy at rest
*/
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    p[]=1;
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
event conv (t += .01; i <= 1000000) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g \n",t,interpolate (u.x, L0/2, 0));
    if (i > 0 && du < 1.0e-6)
        return 1; /* stop */
}  

/**
## Implementation of the Bingham viscosity
*/

event nonnewviscosity(i++) {
    scalar eta_eq[]; 

/** computation of $-II_2 = 2 D:D$  and $D_2=\sqrt{D:D}$
$$ 2 D:D = (2 [(\frac{\partial v}{\partial y})^2  + (\frac{ v}{ y})^2) +(\frac{\partial u}{\partial x})^2] +
   [\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}]^2) $$
   Note that $y$ is $r$

so viscosity is 
$$
\eta_{eq} = 1 + \frac{Bi}{\sqrt{2 D2}}
$$
with regularisation
*/
    
 foreach() { 
      double mI2 = 0.,D2 = 0;
      double duxx = (u.x[1,0] - u.x[-1,0])/(2 * Delta);
      double duxy = (u.x[0,1] - u.x[0,-1])/(2 * Delta);
      double duyx = (u.y[1,0] - u.y[-1,0])/(2 * Delta);
      double duyy = (u.y[0,1] - u.y[0,-1])/(2 * Delta); 
	    mI2 =  sq(duyx+duxy) + 2*(sq(duyy) + sq(duxx) + sq(u.y[]/ max(y, 1e-20)));
      D2 = sqrt(mI2/2.);
      if(D2>0){
        eta_eq[] = min( 1 +  Bi/sqrt(mI2),mumax);}
       else {
        eta_eq[]=mumax;	
       }
    }   
      boundary ({eta_eq});

 
    boundary ({mu_eq});
    foreach_face() {
        muv.x[] = fm.x[]*(eta_eq[] + eta_eq[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}

/**
  Save profiles computed, shear and exact
*/
event profiles (t += 1)
{
    FILE * fp = fopen("xprof", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g %g  \n", y, interpolate (u.x, L0/2, y), interpolate (u.y, L0/2, y),
        	   interpolate (shear, L0/2, y),
                  Ubi(y,Bi));
    fclose (fp);
}
/**
We adapt according to the error on the velocity field. 
*/
event adapt (i++) {
	 fprintf (stderr," %g \n",t);
	  // adapt_wavelet ({u}, (double[]){3e-3,3e-3}, 8, 6);
}

event profile (t = end) {
  foreach()
    printf ("%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
}

/**

## Compilation 

~~~bash
make nonnewtube.tst;
make nonnewtube/plots;
make nonnewtube.c.html;
~~~



## Results and plots

Plot of the exact and computed velocities

~~~gnuplot profiles
set ylabel "u(y)";set xlabel "y" 
p'xprof' t'u comp',''u 1:5 t' exact' w l
~~~




## Bibliography 

* R. Darby Viscoelastic fluids, Dekker ed. (1976) p 223-225, p 194 

* [Application to the 1D Collapse](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)
 
* [Application to the 2D Collapse](http://basilisk.fr/sandbox/M1EMN/Exemples/column_SCC.c)
 
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


