/**
# Errors for the Gaussian vortex pair experiment

Following [this case](/src/test/vortex.c): Two Gaussian vortices with
strength $\omega_0$ and size $\sigma$ are placed at a distance $d$
apart;

$$\omega (x,y) = \omega_0 \left(  e^{-\left(\frac{r_1}{\sigma}\right)^2} +
e^{-\left(\frac{r_2}{\sigma}\right)^2} \right) $$

with $r_1$ and $r_2$ the radial coordinates centered at each
vortex. We choose $\omega_0 = 1 [T^{-1}]$, $\sigma = \frac{1}{2}
[L]$ and $d = 4\sigma$.

![The dynamics of vortex merging](gauss_vor/omega.mp4) 

The aim of this page is to analyze the errors that arise due to the spatial discretization. As such we have a look at the errors in the
tendencyfield $(\frac{\partial \omega}{\partial t})$ at $t = 0$, 
using the $\omega - \psi$ Euler-equation solver.
 */
#include "grid/multigrid.h"
#include "navier-stokes/stream.h"

#define RAD(a) (sqrt(sq(x - a) + sq(y)))
/**
A superior solution is computed at a high reference level. The
tendency obtained from this solution is compared against those
obtained at various lower level grids, upto `maxlevel`. A movie is
generated at `movielevel`.
 */
double * tend;
int reflevel = 12, maxlevel = 10, movielevel = 8;

#define OFFSET(l) ((1 << (2*(l)))*((l) + 1))
#define _O (-GHOSTS)
#define C_IND(i,j,l) ((i+_O) + (1 << l)*(j+_O))
#define INDEX (OFFSET(level - 1) + C_IND(point.i, point.j, level)) 
/**
A function that stores the superior solution at levels in an array:
 */

void tendency_at_levels(scalar tendency) {
  tend = (double*) malloc (OFFSET(maxlevel)*sizeof(double));
  restriction ({tendency});
  foreach_cell()
    if (level <= maxlevel)
      tend[INDEX] = tendency[];
}

int main() {
  TOLERANCE = 1.e-7;
  L0 = 10;
  X0 = Y0 = -L0/2.;
  /**
First, the reference case is `run()`.
   */
  N = 1 << reflevel; 
  run();
  /**
Next, the convergence study is performed.
   */
  for (int l = 5; l <= maxlevel; l++) {
    N = 1 << l;
    run();
  }
  free (tend);
}
/**
For initialization, we approximate the cell-averaged value of
$\omega(t = 0, x, y)$ using 16 points.
 */
static inline double omega_init (Point point) {
  double a = 0;
  foreach_child()
    foreach_child()
      a += exp (-sq(RAD(1)*2.)) + exp (-sq(RAD(-1)*2.));
  return a/16.;
}

event init (t = 0) {
  foreach() 
    omega[] = omega_init (point);
  boundary ({omega});
}
/**
A single small timestep is used to estimate the tendency field.
 */
event stop (t = 0.0001) {
  assert (i == 1);
  scalar tendency[];
  foreach()
    tendency[] = (omega[] - omega_init (point))/dt;
  if (depth() == reflevel) {
    tendency_at_levels (tendency);
    output_ppm (tendency, file = "tend.png", n = 512, min = -0.1, max = 0.1);
    /**
The reference tendency field looks like this

![Reference tendency](gauss_vor/tend.png) 

If the computations are not performed at the reference level, we
compute the error field via the `tend` array.
     */
  } else {
    scalar error[];
    char str[99];
    sprintf (str, "error%d.png", depth());
    foreach()
      error[] = tendency[] - tend[INDEX]; 
    output_ppm (error, file = str, n = 512, min = -1e-4, max = 1e-4);
    /**
The error field is characterized by some non-obvious structure: 

![The error field for the level 8 computation](gauss_vor/error8.png)

We stop after one timestep or continue to generate the movie. 
     */
  }
  if (depth() != movielevel)
    return 1;
}

event movie (t += 1; t < 100) {
  if (depth() == movielevel)
    output_ppm (omega, file = "omega.mp4", min = -1, max = 1);
}
