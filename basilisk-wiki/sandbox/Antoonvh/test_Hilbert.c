/**
# Using the Hilbert iterator

![A simulation](test_Hilbert/movie.mp4)

![Grid-cell processor index](test_Hilbert/movieP.mp4)

~~~gnuplot The grid iterator and `pid()`
set key off
set size square
set cbrange [-0.5:2.5]
plot 'out' w l palette
~~~
 */
#include "quadtree_Hilbert.h"
#include "navier-stokes/centered.h"

double H_indexing (scalar H, bool leaves);
#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

u.t[top] = dirichlet (0);

const face vector nu[] = {0.001, 0.001};

int main() {
  mu = nu;
  L0 = 20;
  X0 = Y0 = -L0/2.;
  N = 256;
  run();
}

event init (t = 0) {
  double k = 3.83170597;
  scalar psi[];
  foreach()
    psi[] = ((RAD > 1)*((1/RAD))*ST +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -(psi[0, 1] - psi[0, -1])/(2*Delta);
    u.y[] = (psi[1] - psi[-1])/(2*Delta);
  }
}

event adapt (i++) 
  adapt_wavelet ((scalar*){u}, (double[]){0.05, 0.05}, 8);

event movie (t += 0.1){
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "movie.mp4", n = 400);
  foreach()
    omega[] = pid();
  output_ppm (omega, file = "movieP.mp4",
	      n = 400, min = 0, max = npe());
}

event stop (t = 20) {
  for (int piz = 0; piz < npe(); piz++) {
    if (pid() == piz)
      foreach()
	printf ("%g %g %d\n", x, y, pid());
    fflush (stdout);
#if _MPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif
  }
}
