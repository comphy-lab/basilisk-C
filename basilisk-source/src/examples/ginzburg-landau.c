/**
# The complex Ginzburg--Landau equation

The complex [Ginzburg--Landau
equation](http://codeinthehole.com/tutorial/index.html)
$$
\partial_t A = A + \left( 1 + i \alpha \right) \nabla^2 A - \left( 1 + i
\beta \right)  \left| A \right|^2 A
$$
with $A$ a complex number, is a classical model for phenomena
exhibiting [Hopf
bifurcations](http://en.wikipedia.org/wiki/Hopf_bifurcation) such as
Rayleigh-Bénard convection or superconductivity.

Posing $A_r=Re(A)$ and $A_i=Im(A)$ one gets the
coupled reaction--diffusion equations.
$$
\partial_t A_r = \nabla^2 A_r + A_r  \left( 1 - \left| A \right|^2 \right)
   - \alpha \nabla^2 A_i + \left| A \right|^2 \beta A_i
$$
$$
\partial_t A_i = \nabla^2 A_i + A_i  \left( 1 - \left| A \right|^2 \right)
   + \alpha \nabla^2 A_r - \left| A \right|^2 \beta A_r
$$

This system can be solved either with the reaction--diffusion solver
or with the generic linear system solver. */

#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"
#include "solve.h"

scalar Ar[], Ai[], A2[];
double alpha = 0., beta = 1.5;
double tend = 150;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion or linear solvers in `mg`. */

double dt;
mgstats mg;

/**
## Parameters

We change the size of the domain `L0` and set the spatial
resolution. */

int main() {
  size (100);
  init_grid (256);

  /**
  The maximum stable timestep depends on the numerical scheme
  used. The weakly-coupled scheme used on GPUs (see
  [below](#weakly-coupled-diffusion-solvers)) is much less stable than
  the generic coupled scheme. */
  
#if _GPU
  DT = 0.05;
#else
  DT = 0.15;
#endif
  
  run();

  /**
  We do a second run, but not on GPUs due to the limitation of the
  [weakly-coupled diffusion
  scheme](#weakly-coupled-diffusion-solvers). */
  
#if !_GPU
  alpha = -3.5, beta = 0.44;
  tend = 600.;
  size (300);
  run();
#endif // !_GPU
}

/**
## Initial conditions 

We use a white noise in $[-10^{-4}:10^{-4}]$ for both components. */

event init (i = 0) {
  foreach() {
    Ar[] = 1e-4*noise();
    Ai[] = 1e-4*noise();
  }
}

/**
## Time integration */

event integration (i++) {

  /**
  We first set the timestep according to the timing of upcoming
  events and the maximum timestep DT. */

  dt = dtnext (DT);
  
  /**
  We consider two different discretisations of the coupled system. The
  first one is the simplest and most accurate and robust, but does not
  work on GPUs for the moment. The second one gives an example of how to
  re-use the diffusion solvers to solve a weakly-coupled approximation:
  it also works on GPUs.

  ### Using the generic linear system solver */

#if !_GPU
  
  /**
  We store the values of $A_i(t)$ and $A_r(t)$. */
  
  scalar Art[], Ait[];
  foreach() {
    Art[] = Ar[];
    Ait[] = Ai[];
  }

  /**
  Since the generic linear system solver uses a multilevel strategy,
  we need to make sure that all fields (other than the unknowns)
  appearing in the expressions below are defined on all levels. */
  
  restriction ({Art, Ait});

  /**
  We solve for $(A_r,A_i)(t+dt)$ using a first-order, implicit-in-time scheme
  and store the multigrid statistics in `mg`. */
  
  #define laplacian(a) ((a[1,0] + a[-1,0] + a[0,1] + a[0,-1] - 4.*a[])/sq(Delta))
  
  msolve ({Ar, Ai}, {
      (Art[] - Ar[])/dt + laplacian(Ar) + Ar[]*(1. - A2[]) - alpha*laplacian(Ai) + A2[]*beta*Ai[],
      (Ait[] - Ai[])/dt + laplacian(Ai) + Ai[]*(1. - A2[]) + alpha*laplacian(Ar) - A2[]*beta*Ar[]
    }, &mg) {

    /**
    Since $|A(t+dt)|^2$ is a non-linear function of $(A_r,A_i)(t+dt)$
    we cannot discretise it directly, but we can update its value
    after each iteration of the multigrid solver, using the "ellipsis"
    argument of the `msolve()` macro. */
    
    foreach()
      A2[] = sq(Ar[]) + sq(Ai[]);
    restriction ({A2});
  }
  
  /**
  ### Weakly-coupled diffusion solvers */

#else // _GPU

  /**
  We compute $|A|^2$ and the other fields necessary as inputs to the 
  diffusion solver. */

  scalar r[], lambda[];
  foreach() {
    A2[] = sq(Ar[]) + sq(Ai[]);
    r[] = A2[]*beta*Ai[];
    lambda[] = 1. - A2[];
  }

  /**
  The scheme is less general than the generic one and works only for
  $\alpha = 0$. */

  assert (alpha == 0.);
  
  /**
  We use the diffusion solver (twice) to advance the system from $t$
  to $t+dt$. */

  diffusion (Ar, dt, r = r, beta = lambda);

  /**
  We use a time-split reaction term. Note that this seems to cause
  spurious oscillations in time in $|A|^2$ which are absent when using
  the coupled approach. */

  foreach() {
    r[] = - A2[]*beta*Ar[];
    lambda[] = 1. - A2[]; // this is necessary
  }
  mg = diffusion (Ai, dt, r = r, beta = lambda);
  
#endif // _GPU
}

/**
## Outputs

Here we create MP4 animations for both components. The `spread`
parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movies (t += tend/1000; t <= tend) {
  fprintf (stderr, "%g %g %g %d\n", t, dt, sqrt(normf(A2).max), mg.i);

  char name[80];
  sprintf (name, "Ai-%g-%g.mp4", alpha, beta);
  output_ppm (Ai, spread = 2, linear = true, file = name);
  sprintf (name, "A2-%g-%g.mp4", alpha, beta);
  output_ppm (A2, spread = 2, linear = true, file = name);
}

/**
For $\alpha=0$ and $\beta=1.5$ we get a typical "frozen state"
composed of "cellular structures" for $|A|^2$ and stationary spirals
for $A_i$.

<center><table>
<tr>
<td>![](ginzburg-landau/A2-0-1.5.mp4)(autoplay)</td>
<td>![](ginzburg-landau/Ai-0-1.5.mp4)(autoplay)</td>
</tr>
<tr><td>$|A|^2$</td>  <td>$A_i$</td></tr>
<caption>
Evolution of the norm and imaginary part for $\alpha=0$ and $\beta=1.5$.
</caption>
</table></center>

For $\alpha=-3.5$ and $\beta=0.44$ we observe the nucleation of
spiralling patterns that compete and finally invade the whole domain.

<center><table>
<tr>
<td>![](ginzburg-landau/A2--3.5-0.44.mp4)(autoplay)</td>
<td>![](ginzburg-landau/Ai--3.5-0.44.mp4)(autoplay)</td>
</tr>
<tr><td>$|A|^2$</td>  <td>$A_i$</td></tr>
<caption>
Evolution of the norm and imaginary part for $\alpha=-3.5$ and $\beta=0.44$.
</caption>
</table></center>

## See also

* [The Brusselator](brusselator.c). */
