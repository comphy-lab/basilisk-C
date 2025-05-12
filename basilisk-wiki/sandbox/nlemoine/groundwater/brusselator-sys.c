/**
# Coupled reaction--diffusion equations

This is a test of a [Poisson-Helmoltz solver with coupled equations](http://basilisk.fr/sandbox/nlemoine/groundwater/diffusion-crossreact.h), it is based on the [Brusselator](/src/examples/brusselator.c) test case.

Two chemical compounds with concentrations $C_1$ and $C_2$ interact
according to the coupled reaction--diffusion equations:
$$
\partial_t C_1 = \nabla^2 C_1 + k(ka - (kb + 1)C_1 + C_1^2 C_2)
$$
$$
\partial_t C_2 = D \nabla^2 C_2  + k(kb C_1 - C_1^2 C_2)
$$

We will use a Cartesian (multi)grid, the generic time loop and the
time-implicit diffusion solver. */

#include "grid/multigrid.h"
#include "run.h"
//#include "diffusion-crossreact.h"
#include "diffusion.h"

/**
We need scalar fields for the concentrations. */

scalar C1[], C2[];

/**
We use the same parameters as [Pena and Perez-Garcia,
2001](/src/references.bib#pena2001) */

double k = 1., ka = 4.5, D = 8.;
double mu, kb;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
mgstats mgd,mgd1,mgd2;

/**
## Parameters

We change the size of the domain `L0` and set the tolerance of the
implicit diffusion solver. */

int main()
{
  init_grid (128);
  size (64);
  TOLERANCE = 1e-4;

  /**
  Here $\mu$ is the control parameter.  For $\mu > 0$ the system is
  supercritical (Hopf bifurcation). We test several values of $\mu$. */

  mu = 0.04; run();
  mu = 0.1;  run();
  mu = 0.98; run();
}

/**
## Initial conditions */

event init (i = 0)
{

  /**
  The marginal stability is obtained for `kb = kbcrit`. */

  double nu = sqrt(1./D);
  double kbcrit = sq(1. + ka*nu);
  kb = kbcrit*(1. + mu);

  /**
  The (unstable) stationary solution is $C_1 = ka$ and $C_2 = kb/ka$. It
  is perturbed by a random noise in [-0.01:0.01]. */

  foreach() {
    C1[] = ka ; 
    C2[] = kb/ka + 0.01*noise();
  }
}

/**
## Outputs

Here we create an mpeg animation of the $C_1$ concentration. The
`spread` parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movie (i = 1; i += 10)
{
  output_ppm (C1, linear = true, spread = 2, file = "f.mp4", n = 200);
  fprintf (stderr, "%d %g %g %d\n", i, t, dt, mgd.i);
}

/**
We make a PNG image of the final "pseudo-stationary" solution. */

event final (t = 3000)
{
  char name[80];
  sprintf (name, "mu-%g.png", mu);
  output_ppm (C1, file = name, n = 200, linear = true, spread = 2);
}

/**
## Time integration */

event integration (i++)
{

  /**
  We first set the timestep according to the timing of upcoming
  events. We choose a maximum timestep of 1 which ensures the stability
  of the reactive terms for this example. */

  dt = dtnext (1.);

  /**
  In order to define the coupling terms, we rearrange the evolution equations. We first define $\rho_1$ and $\rho_2$ such that:
  
$$
\partial_t C_1 = \nabla^2 C_1 + \underbrace{k(k_a - (k_b + 1)C_1 + C_1^2 C_2)}_{\normalsize\rho_1}
$$

$$
\partial_t C_2 = D \nabla^2 C_2  + \underbrace{k(k_b C_1 - C_1^2 C_2)}_{\normalsize\rho_2}
$$

and then we write:

$$
\rho_1^{(n+1)}\quad \approx\quad \rho_1^{(n)}\quad +\quad \left(\frac{\partial\rho_1}{\partial C_1}\right)^{(n)}\left(C_1^{(n+1)}-C_1^{(n)}\right)
\quad +\quad \left(\frac{\partial\rho_1}{\partial C_2}\right)^{(n)}\left(C_2^{(n+1)}-C_2^{(n)}\right)
$$

$$
\rho_2^{(n+1)}\quad \approx\quad \rho_2^{(n)}\quad +\quad \left(\frac{\partial\rho_2}{\partial C_1}\right)^{(n)}\left(C_1^{(n+1)}-C_1^{(n)}\right)
\quad +\quad \left(\frac{\partial\rho_2}{\partial C_2}\right)^{(n)}\left(C_2^{(n+1)}-C_2^{(n)}\right)
$$

The coupling coefficients are then
$$\beta_{11} = \frac{\partial\rho_1}{\partial C_1} = k\left(2\,C_1 C_2 -k_b-1\right) $$

$$\beta_{12} = \frac{\partial\rho_1}{\partial C_2} = k\,C_1^2 $$

$$\beta_{21} = \frac{\partial\rho_2}{\partial C_1} = k\left(k_b-2\,C_1 C_2\right) $$

$$\beta_{22} = \frac{\partial\rho_2}{\partial C_2} = -k\,C_1^2 $$

so that finally:

$$
\rho_1^{(n+1)}\quad=\quad\beta_{11}C_1^{(n+1)} + \beta_{12}C_2^{(n+1)} 
+ \underbrace{\left[\rho_1^{(n)}-\beta_{11}C_1^{(n)}-\beta_{12}C_2^{(n)}\right]}_{\normalsize r_1}
$$

$$
\rho_2^{(n+1)}\quad=\quad\beta_{21}C_1^{(n+1)} + \beta_{22}C_2^{(n+1)} 
+ \underbrace{\left[\rho_2^{(n)}-\beta_{21}C_1^{(n)}-\beta_{22}C_2^{(n)}\right]}_{\normalsize r_2}
$$
*/

  scalar r1[], r2[];
  scalar beta11[],beta12[],beta21[],beta22[];
  scalar theta1[],theta2[];
  face vector D1[],D2[];  

  foreach_face()
  {
    D1.x[] = 1.;
    D2.x[] = D;  
  }

  /* Cross-reaction matrix */
  
  double coupling = 0.;
  
  foreach()
  {
    r1[] = k * ( ka-(kb+1.)*C1[] + sq(C1[])*C2[] );
    beta11[] = k * ( 2.*C1[]*C2[] - kb - 1. );
    beta12[] = coupling * k * sq(C1[]);  
    r1[] -= ( beta11[]*C1[] + coupling * beta12[]*C2[] );
      
    r2[] = k * ( kb*C1[] - sq(C1[])*C2[] );
    beta21[] = coupling * k * ( kb - 2.*C1[]*C2[] );
    beta22[] = -k * sq(C1[]);  
    r2[] -= ( coupling * beta21[]*C1[] + beta22[]*C2[] );
      
    theta1[] = 1.;
    theta2[] = 1.;
  }
    
  face vector * Dl = {D1,D2};
  scalar * thetal = {theta1,theta2};
    
  scalar * MAT_BETA[]= {{beta11,beta12},
                        {beta21,beta22}};

  scalar * r = {r1,r2};
  scalar * C = {C1,C2};

  fprintf(stderr,"i = %d, done preparing inputs to sysdiffusion.\n",i);
  fflush(stderr);    
    
//  mgd = sysdiffusion(C, dt, Dl, r, MAT_BETA, thetal);
  mgd1 = diffusion (C1, dt, D1, r1, beta11,theta1);
  mgd2 = diffusion (C2, dt, D2, r2, beta22,theta2);
}

/**
![log](brusselator-sys/log)
*/

/**
## Results

We get the following stable [Turing
patterns](http://en.wikipedia.org/wiki/The_Chemical_Basis_of_Morphogenesis).

<center>
<table>
<tr>
<td>![](brusselator-sys/mu-0.04.png)</td>
<td>![](brusselator-sys/mu-0.1.png)</td>
 <td>![](brusselator-sys/mu-0.98.png)</td>
</tr>
<tr>
<td>$\mu=0.04$</td> 
<td>$\mu=0.1$ (stripes)</td> 
<td>$\mu=0.98$ (hexagons)</td>
</tr>
</table>
</center>

![Animation of the transitions](brusselator-sys/f.mp4)
*/
