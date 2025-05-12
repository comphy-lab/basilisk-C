/**
<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/jsct9cv17k22l0y/SchematicCouette.png?dl=0" width="30%">
  <figcaption><p align="center">Schematic of the Planar Couette Flow.</figcaption>
</figure>

#  Planar Couette flow of Generalized Newtonian Fluid

This code extends the method used in [/sandbox/M1EMN/Exemples/bingham_simple.c](/../sandbox/M1EMN/Exemples/bingham_simple.c)
and generalizes it for any Power Law fluid (using regularization method). Another
difference between the two is that this code calculates the second invariant of
deformation tensor at the face-centers of the cells instead of the cell centers.

## Mathematical Formulations

Unlike the Newtonian fluids, non-Newtonian fluids do not have a linear stress-strain rate relationship.
One way to represent the relationship is using the Generalized Newtonian fluid method:
$$ \tau = \tau_y + 2\mu_0D_{ij}^n $$
The fluid is such that
if $\|\tau\| \le  \tau_y$ then there is no motion $D_{ij}=0\:\forall\:(i,j)$
if the stress is high enough $\|\tau\| >  \tau_y$ then there is motion

*Note:* that $\|\tau\|$ is the modulus defined as the Euclidian norm  $\sqrt{\frac{1}{2}{\tau_{ij} \tau_{ij}}}$.
 It is not $\sqrt{\tau_{11}^2 + \tau_{12}^2}$ as in Balmorth et al. (2006), which is the Frobenius norm.

$D_{ij}$ is the shear strain rate tensor (or the deformation tensor)

$D_{ij}=(u_{i,j}+u_{j,i})/2$: the components in 2D:
$$D_{11}=\frac{\partial u}{\partial x}$$
$$D_{12} =\frac{1}{2}\left(\frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{21} =D_{12} =\frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{22}=\frac{\partial v}{\partial y}$$

In the Euclidian norm we have:
$$\|D\|=\sqrt{\frac{D_{ij}D_{ij}}{2}}$$
The second invariant defined by $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
is given by:
$$D_2^2= D_{ij}D_{ij}= \left( \frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2
 +  \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)^2$$
and we have obviously $\|D_{ij}\| = D_2/\sqrt{2}$

## Numerical regularization
$$ \tau_{ij} = \tau_y\left(\frac{D_{ij}}{\|D_{ij}\|}\right) + 2\mu_0\|D_{ij}\|^{n-1}D_{ij}^n $$
Factorising with $2D_{ij}$ to obtain a equivalent viscosity
$$\tau_{ij} = 2\left(\mu_0 \|D_{ij}\|^{n-1} + \frac{\tau_y}{2 \|D_{ij}\|}\right)D_{ij}$$
$$\tau_{ij} = 2 \mu_{eq}D_{ij}$$
$$\mu_{eq} = \mu_0\|D_{ij}\|^{n-1} + \frac{\tau_y}{2\|D_{ij}\|}$$
$\mu$ is the min of $\mu_{eq}$ and a large $\mu_{max}$ so that the viscosity does not blow up.
$$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
*Note:* We present here the formulation in Balmforth, he uses $\dot{\gamma}$ which is by his definition $\sqrt{\frac{1}{2}\dot{\gamma_{ij}}\dot{\gamma_{ij}}}$
and as $\dot{\gamma_{ij}}=2 D_{ij}$ then $\dot{\gamma}$ is $\sqrt{2}D_2$, that is why we have a $\sqrt{2}$ in the equations.





##  Exact solution in the proposed case


We look at an unidirectional flow, a pure shear flow  $u(y)$, $v=0$, so
$D_{11}=D_{22}=0$ and $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$.
$$ \tau_{12} = 2\mu D_{12}^n  + \tau_y =
  2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tau_y $$

Equilibrium between pressure gradient and viscosity (writting $\tau$ for a shorthand of $\tau_{12}$)
$$0=-\frac{\partial p}{\partial x} + \frac{\partial \tau}{\partial y}$$
as there is no stress at the free surface $y=h$, the stress is
$$ \tau = \left(-\frac{\partial p}{\partial x}\right)(h-y)$$
the stress $\tau$ increases from the free surface, as long as $\tau<\tau_y$,
we are under the threshold,
so shear is zero: $\frac{\partial u}{\partial y} =0$,  hence velocity is constant, say it is $U$.
Let us define
 $Y=h-\tau_y/(-\frac{\partial p}{\partial x})$, where  $\tau=\tau_y$.


So :
$$ \left\{\tau<\tau_y, \frac{\partial u}{\partial y} = 0,\:\&\:u=U\:\forall\:Y<y<h\right\} $$
Then going down:
$0<y<Y$ we have $\tau = 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tau_y$.
This gives:
$$\tau_y + 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n =  \left(-\frac{\partial p}{\partial x}\right)(Y-y)$$
After some straight forward manipulations:
$$\left(\frac{\partial u}{\partial y}\right) =  \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}(Y-y)^{1/n}
 = A_{p\mu}(Y-y)^{1/n}$$
and this allows to solve for the velocity profile
$$u = \frac{nA_{p\mu}}{n+1}\left(Y^{\frac{n+1}{n}} - \left(Y-y\right)^{\frac{n+1}{n}}\right)$$
which is indeed zero in $y=0$, and for  $y=Y$, we have the plug flow $u=U$ of value:
$$U= \frac{n}{n+1}A_{p\mu}Y^{\frac{n+1}{n}}$$
$$A_{p\mu} = \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}$$
For the present case $-\frac{\partial p}{\partial x} = 1, \mu = 1, h = 1$, which gives:
$A_{p\mu} = \frac{1}{2^{(1-n)/n}}, Y = 1 - \tau_y$

# Code
*/
#include "navier-stokes/centered.h"
char filename[80];
double tauy,mu_0,mumax;
double n;
int imax = 1e4;
#define dtmax (1e-3)
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  init_grid (1<<6);
  L0 = 1.0;
  origin (0.0, 0.0);
  DT = dtmax;
  stokes = true;
  TOLERANCE = 1e-5;

/**
 Values of yeild stress, viscosity, and coefficient.<br/>
 Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1<br/>
 Power law $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 0.5<br/>
 Herschel-Bulkley $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 0.5<br/>
 Bingham $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 1<br/>
*/

  mu_0 = 1.0;
  tauy= 0.0;
  n = 1.0;
  if (a >= 3){
    mu_0 = atof(arguments[2]);
  }
  if (a >= 4){
    tauy = atof(arguments[3]);
  }
  if (a >= 5){
    n = atof(arguments[4]);
  }

/**
  the regularisation value of viscosity
*/
  mumax=1000;

/**
 Right - left boundaries are periodic
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
 presure conditions are neumann 0.0
 */
  p[top] = neumann(0);
  pf[top] = neumann(0);
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);

  run();
}

// un is used to search for a stationary solution
scalar un[];
// muv will be used as the face vector for viscosity
face vector muv[];

/**
## Initialization event
*/
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  /**
    presure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0,0.0};
  a = mdpdx;
  /**
   Initialy at rest
  */
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
  foreach(){
    un[] = u.x[];
  }
  dump (file = "start");
}

/**
We look for a stationary solution. */
event logfile (i += 500; i <= imax) {
  double du = change (u.x, un);
  fprintf(ferr, "i = %d: err = %g\n", i, du);
  if (i > 0 && du < 1e-6){
    dump (file = filename);
    return 1; /* stop */
  }
  if (i==imax){
    dump (file = filename);
  }
}


event properties(i++) { // Overloading the properties event
  /**
  ## Implementation of generalized Newtonian viscosity

   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{N-1} + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$

   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$.

   The fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
  */
  double muTemp = mu_0;
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    if (D2 > 0.0) {
      double temp = tauy/(sqrt(2.0)*D2) + mu_0*exp((n-1.0)*log(D2/sqrt(2.0)));
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.0 || n < 1.0){
        muTemp = mumax;
      } else {
        muTemp = (n == 1.0 ? mu_0 : 0.0);
      }
    }
    muv.x[] = fm.x[]*(muTemp);
  }
  boundary ((scalar *){muv});
}

/**
## Running the code

Use the following `run.sh` script

~~~bash
#!/bin/bash
qcc -O2 -Wall Couette_NonNewtonian.c -o Couette_NonNewtonian -lm
./Couette_NonNewtonian lastNewt 1.0 0.0 1.0
./Couette_NonNewtonian lastShThn 1.0 0.0 0.5
./Couette_NonNewtonian lastHB 1.0 0.25 0.5
./Couette_NonNewtonian lastBing 1.0 0.25 1.0
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://www.dropbox.com/sh/8at1yk9vigovdg7/AAAr-Td7p106Kt_3cIK4mg_ia?dl=0)
<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/vmznrcg07xlxjd1/NonNewtonian_FACE.png?dl=0" width="50%">
  <figcaption><p align="center">Velocity and shear rate for Generalized Newtonian Fluids in plannar Couette flow.</figcaption>
</figure>

<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/mq1sj7px3pkwfqy/Bingham_FACE.png?dl=0" width="50%">
  <figcaption><p align="center">Velocity, shear rate and second invariant of deformation tensor $\|D_{ij}\|$.</figcaption>
</figure>

<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/zze7ajcj4g7lrd6/D2Bingham.png?dl=0" width="50%">
  <figcaption><p align="center">$\|D_{ij}\|$ contour for Bingham Fluid</figcaption>
</figure>

# Bibliography

* [Same example in Basilisk using the calculation of D2 at cell centers](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_simple.c)
  and its application to [1D Collapse](http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_collapse_noSV.c)

* [Related example in Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/couette.html)

* [Related example with augmented Lagrangian](http://basilisk.fr/sandbox/popinet/poiseuille-periodic.c)

*  K. F. Liu and C. C. Mei
 Liu, K.F. and Mei, C.C., 1990. Approximate equations for the slow spreading of a thin sheet of Bingham plastic fluid.
 Physics of Fluids A: Fluid Dynamics, 2(1), pp.30-36.; [doi: 10.1063/1.857821](https://doi.org/10.1063/1.857821)

* The Theoretical Formulations: Bird, R.B., 1987. Armstrong and RC Hassager, O.,“Dynamics of Polymeric Liquids”. v.1.

* Balmforth, N.J., Craster, R.V., Rust, A.C. and Sassi, R., 2006. Viscoplastic flow over an inclined surface. Journal of Non-Newtonian Fluid Mechanics, 139(1-2), pp.103-127.
  [doi: 10.1016/j.jnnfm.2006.07.010](https://doi.org/10.1016/j.jnnfm.2006.07.010)
*/
