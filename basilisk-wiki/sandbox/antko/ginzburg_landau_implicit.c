/**
# Implicit solving of the complex Ginzburg--Landau equation

We revisit the [Ginzburg-Landau example](/src/examples/ginzburg-landau.c) with the [coupled diffusion solver](./coupled_diffusion.h). 
We solve the following coupled evolution equations:
$$
\partial_t A_r = \nabla^2 A_r + A_r  \left( 1 - \left| A \right|^2 \right)
   - \alpha \nabla^2 A_i + \left| A \right|^2 \beta A_i
$$
$$
\partial_t A_i = \nabla^2 A_i + A_i  \left( 1 - \left| A \right|^2 \right)
   + \alpha \nabla^2 A_r - \left| A \right|^2 \beta A_r
$$
*/

#include "grid/multigrid.h"
#include "run.h"
#include "coupled_diffusion.h"

scalar Ar[], Ai[], A2[];
double alpha, beta;
double dt, dt_run = 0.05 [0,1], tend = 500. [0,1];
double domain_size;
mgstats mgd;
int nrun = 1;

/**
## Parameters
*/

int main() {
/** We first consider the set of parameters used in [ginzburg-landau.c](/src/examples/ginzburg-landau.c) */
  alpha = 0.;
  beta = 1.5;
  dt_run = 0.05;
  tend = 150.;
  size (100);
  init_grid (256);
  run();
/**
Here we recover the frozen state depicted in [Ginzburg-Landau example](/src/examples/ginzburg-landau.c) albeit with less emaning waves from the corners.

<center><table>
<tr>
<td>![](ginzburg_landau_implicit/A2_1.mp4)(autoplay)</td>
<td>![](ginzburg_landau_implicit/Ai_1.mp4)(autoplay)</td>
</tr>
<tr><td>$|A|^2$</td>  <td>$A_i$</td></tr>
<caption>
Evolution of the norm and imaginary part
</caption>
</table></center>
*/
/** In order to probe the stability of the implicit scheme, we rerun the same case but with a larger timestep. */
  nrun++;
  dt_run = 0.5;
  size (100);
  init_grid (256);
  run();
/**
Here is a comparison between the solutions obtained with two different timesteps.
<center><table>
<tr>
<td>![](ginzburg_landau_implicit/A2_1.mp4)(autoplay)</td>
<td>![](ginzburg_landau_implicit/A2_2.mp4)(autoplay)</td>
</tr>
<tr><td>$|A|^2$ for dt = 0.05</td>  <td>$|A|^2$ for dt = 0.5</td></tr>
</table></center>
*/
/** Finally, we explore another set of parameters (where the equations are fully coupled) prone to nucleation events. */
  nrun++;
  alpha = -3.5;
  beta = 0.44;
  dt_run = 0.05;
  tend = 500.;
  size (300);
  init_grid (256);
  run();
/**
For this different set of parameters, we observe the nucleation of spiralling patterns that compete and finally invade the whole domain.
<center><table>
<tr>
<td>![](ginzburg_landau_implicit/A2_3.mp4)(autoplay)</td>
<td>![](ginzburg_landau_implicit/Ai_3.mp4)(autoplay)</td>
</tr>
<tr><td>$|A|^2$</td>  <td>$A_i$</td></tr>
<caption>
Evolution of the norm and imaginary part
</caption>
</table></center>
*/
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

  dt = dtnext (dt_run);

  /**
  We compute $|A|^2$. */

  foreach()
    A2[] = sq(Ar[]) + sq(Ai[]);

  /**
  We use the coupled diffusion solver (once) to advance the system from $t$
  to $t+dt$. */

  scalar lambda1[], lambda2[], mu1[], mu2[];
  face vector DArAi[], DAiAr[];
  foreach() {
    lambda1[] = 1. - A2[];
    mu1[]     = A2[]*beta;
    lambda2[] = -A2[]*beta;
    mu2[]     = 1. - A2[];
  }
  foreach_face() {
    DArAi.x[] = -alpha;
    DAiAr.x[] = +alpha;
  }

  mgd = coupled_diffusion (Ar, Ai, dt, Dab = DArAi, Dba = DAiAr, 
                           alphaa = lambda1, betaa = mu1,
                           alphab = lambda2, betab = mu2);
}
/**
## Outputs

Here we create MP4 animations for both components. The `spread`
parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movies (i += 3; t <= tend) {
  fprintf (stderr, "%g %g\n", t, sqrt(normf(A2).max));
  char fnameAi[80];
  char fnameA2[80];
  sprintf(fnameAi, "Ai_%d.mp4", nrun);
  output_ppm (Ai, spread = 2, linear = true, file = fnameAi);
  sprintf(fnameA2, "A2_%d.mp4", nrun);
  output_ppm (A2, spread = 2, linear = true, file = fnameA2);
}


