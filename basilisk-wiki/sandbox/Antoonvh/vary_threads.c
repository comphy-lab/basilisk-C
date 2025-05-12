/**
# Varying the number of threads

The flow problem:
![Evolution of the vorticity field](vary_threads/omega.mp4)

The adaptive grid-level refinment

![quadtree level of refinment](vary_threads/level.mp4)

The number of threads varies over time as can be seen form the decomposition

![Thread ID distribution](vary_threads/tids.mp4)



 */
#include "navier-stokes/centered.h"

u.t[top] = dirichlet (0);

int main() {
  periodic (left);
  L0 = 15;
  X0 = -L0/2;
  Y0 = - 10;
  const face vector muc[] = {0.001, 0.001};
  mu = muc;
  N = 256;
  run();
}

event init (t = 0) {
  foreach() 
    u.y[] = exp(-sq(y) - sq(x));
}


int MAX_NUM, cells_per_thread = 500;

event change_thread_num (i++) {
  if (!MAX_NUM)
    MAX_NUM = omp_get_max_threads();
  if (grid->tn/omp_get_max_threads() > cells_per_thread &&
      omp_get_max_threads() < MAX_NUM) 
    omp_set_num_threads (omp_get_max_threads() + 1);
  else if (omp_get_max_threads() > 1 &&
	   grid->tn/(omp_get_max_threads() - 1) < cells_per_thread)
    omp_set_num_threads (omp_get_max_threads() - 1);
}

event output (t += 1) {
  printf ("%g %d %d\n",t, omp_get_max_threads(), grid->tn/omp_get_max_threads());
  scalar ti[];
  foreach()
    ti[] = tid();
  output_ppm (ti, file = "tids.mp4", min = 0, max = MAX_NUM, n = 400);
  vorticity (u, ti);
  output_ppm (ti, file = "omega.mp4", min = -.5, max = .5, n = 400, map = cool_warm);
  
  foreach()
    ti[] = level;
  output_ppm (ti, file = "level.mp4", min = 2, max = 9, n = 400);
}

event adapt (i++) {
  adapt_wavelet ({u.x, u.y}, (double[]){0.005, 0.005}, 9);
}

event stop (t = 150);
