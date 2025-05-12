/**
# A vortex instability

We compute the vorticity dynamics of an isolated vortex stucture. Such
structures may be unstable (Trieling et al. , 2010). We study a vortex
with the following angular (azimutal) velocity profile ($\mathbf{u}
= \{u_r, u_{\theta}\} =\{0, f(r)\}$), 

$$r < \pi :\ \ u_{\theta} =  \mathrm{sin}^3(r),$$
$$r \geq \pi :\ \ \  u_{\theta} =  0,$$

that is sufficiently smooth for future numerically-oriented studies.
*/
#include "navier-stokes/stream.h"

int maxlevel = 11;

int main() {
  L0 = 10.*pi;
  X0 = Y0 = -L0/2.;
  N = (1 << maxlevel);
  run();
}
/** 
## Setup

The solver requires the vorticity distribution to be initialized: 

$$r \leq \pi: \omega (r) = \frac{\mathrm{sin}^2(r) \times \left( \mathrm{sin}(r) + 3r\mathrm{cos}(r)\right)}{r}, $$
$$r > \pi: \omega (r) = 0.$$

Further, we ad a small perturbation to trigger a reproduceble mode-two
instability. 
 */
#define PERTURB (0.0005)
#define _X (x/(1 + PERTURB))
#define _Y (y*(1 + PERTURB))
#define R (sqrt(sq(_X) + sq(_Y)))
#define OMEGA ((R <= pi)*(sq(sin(R))*			\
			  (sin(R) + 3.*R*cos(R)))/R)

event init (t = 0) {
  TOLERANCE = 1e-4;
  foreach() {
    double omg = 0;
    foreach_child() foreach_child()
      omg += OMEGA;
    omega[] = omg/16.;
  } 
}
/**
The grid is adaptively refined and coarsend,
*/
#if TREE
event adapt (i++)
  adapt_wavelet ({psi}, (double[]){0.002}, maxlevel);
#endif

/**
And a movie is generated:
*/
event movie (t += 0.5; t <= 50) {
  output_ppm (omega, file = "mov.mp4", n = 512,
	      map = cool_warm, min = -1, max = 1);
#if TREE
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "level.mp4", n = 512,
	      min = 3, max = maxlevel);      
#endif
}

event stop (t = 50);
/**
## Results

A movie of the evolution of the vorticity distribution reveals the
dynamics:

![A mode-two instability can be observed](tripolar_vortex/mov.mp4)

The evolution fo the grid structure is as follows:

![The adaptive grid refinement](tripolar_vortex/level.mp4)

## Reference:

Trieling, R. R., Van Heijst, G. J. F., & Kizner, Z. (2010). Laboratory experiments on multipolar vortices in a rotating fluid. Physics of Fluids, 22(9), 094104.
 */