/**
#  Planar Couette flow of Generalized Newtonian Fluid

This is the same case as done in [lid-bingham.c](/sandbox/popinet/lid-bingham.c)
The only difference between the two is that this code uses the method of regularization of stress-strain relationship
instead of the Augmented Lagrangian Method. The second invariant of stress tensor is calculated at the face-centers, similar to
[Couette_NonNewtonian.c](/sandbox/vatsal/GenaralizedNewtonian/Couette_NonNewtonian.c)

## Code
*/
#include "navier-stokes/centered.h"

/**
  Values of yield stress, viscosity, and coefficient.
  Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1
  Bingham $\mu_0 = 1.0$; $\tau_y > 0.0$ and n = 1
*/
char filename[80];
double tauy;
int counter;
#define mu_0 (1.0)
// The regularisation value of viscosity
double mumax = (1e3);

int imax = 1e5;
#define LEVEL 6
#define Maxdt (1e-4)
#define ERROR (1e-6)

scalar un[];
face vector muv[];

// initialization event
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
  foreach(){
    un[] = u.x[];
  }
  dump (file = "start");
}


// int main(int arg, char const *arguments[])
int main()
{
  // sprintf (filename, "%s", arguments[1]);
  // tauy = 1000/sqrt(2.0);
  init_grid (1<<LEVEL);
  L0 = 1.0;
  origin (-0.5, -0.5);
  DT = Maxdt;
  TOLERANCE = 1e-5;
  // CFL number
  CFL = 0.25;
  for (counter = 0; counter < 5; counter++){
    if (counter == 0){
      tauy = 0.0;
    } else {
      tauy = pow(10, counter-1)/sqrt(2);
    }
    if (counter>1){
      mumax = (1e4);
    }
    fprintf(ferr, "tauy = %g\n", tauy);
    sprintf (filename, "tau%d", counter);
    run();
  }
  // run();
}

// Top moving wall
u.t[top] = dirichlet(1);
/**
For the other no-slip boundaries this gives
*/
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

/**
We look for a stationary solution every 500 time steps
*/
event logfile (i=i+500; i <= imax) {
  double du = change (u.x, un);
  fprintf(ferr, "i = %d: dt = %g err = %g\n", i, dt, du);
  if (i > 0 && du < ERROR){
    dump (file = filename);
    return 1; /* stop */
  }
  if (i>imax-10){
    dump (file = filename);
  }
}

event properties(i++) {
  /**
  ## Implementation of generalized Newtonian viscosity
   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$
   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
  */
  double muTemp = mu_0;
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    if (D2 > 0.) {
      double temp = tauy/(sqrt(2.0)*D2) + mu_0;
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.0){
        muTemp = mumax;
      } else {
        muTemp = mu_0;
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
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 LidDrivenBingham.c -o LidDrivenBingham -lm
mpirun -np 6 ./LidDrivenBingham
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://www.dropbox.com/sh/n7d4bdss9zhpogw/AAAkkKuPUH--43IXtmcgs3yCa?dl=0)

<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/kbd5tugi888ip2l/LidDriven.png?dl=0?dl=0" width="100%">
  <figcaption><p align="center">Comparision of the present case with
  [/sandbox/M1EMN/Exemples/bingham_simple.c](/../sandbox/popinet/lid-bingham.c)
  and [Vola et al. (2003)](https://doi.org/10.1016/S0021-9991(03)00118-9).</figcaption>
</figure>

<figure>
<p align="center">
  <img src="https://dl.dropboxusercontent.com/s/mgrov5q5v4flvj2/tau3.png?dl=0" width="100%">
  <figcaption><p align="center">$\|D_{ij}\|$ contour and streamlines of the flow for case: $\tau_y = 100/\sqrt{2}$</figcaption>
</figure>

# Bibliography

* [Same example in Basilisk using the Augmented Lagrangian Method](http://basilisk.fr/sandbox/popinet/lid-bingham.c)

* Vola, D., Boscardin, L. and Latché, J.C., 2003. Laminar unsteady flows of Bingham fluids: a numerical strategy and some benchmark results. Journal of Computational Physics, 187(2), pp.441-456.
  [doi: 10.1016/S0021-9991(03)00118-9](https://doi.org/10.1016/S0021-9991(03)00118-9)
*/
