/**
# [Witzig (1914)](#references) flow

This pages mimics [other
pages](http://basilisk.fr/_search?patterns=Womersley) on pulsatile
flow in a tube.

![The numerical solution and the analytical one](witzig/mov.mp4)

Search words: Womersley, womersley
 */

#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"
// The solution is written as a complex expression
#include <complex.h>

double Rad = 0.01;
double omega = 2*pi;
double dpdz = 5000;
double rho = 1000;
double mu; // = f(Wo);
double Wo = 5;
#define WOMERSLEY  (creal (prefactor*					\
			   (1 - (bessel0(cpow(I, 3./2.)*alpha*x/Rad, 100)/bessel0(cpow(I, 3./2.)*alpha, 100))) \
			   *cexp(I*omega*t)))

scalar u[];
u[right] = dirichlet (0);

int main () {
  mu = sq(Rad/Wo)*omega*rho;
  N = 32;
  L0 = Rad;
  DT = 2*pi/(omega * 40);
  run();
}

event dt_next(i++) {
  dt = dtnext (DT);
}

// complex Bessel function
double complex bessel0 (double complex arg, int iters) {
  double complex res = 0;
  for (int j = 0; j < iters; j++) {
    res += cpow (-1, j)/(sq(tgamma (j + 1)))*cpow(arg/2., 2*j);
  }
  return res;
}

event init (t = 0) {
  double alpha = Rad*sqrt(omega*rho/mu);
  double complex prefactor = I*dpdz/(rho*omega); 
  foreach() 
    u[] = WOMERSLEY;
}

event advance (i++, last) {
  face vector nuf[];
  foreach_face()
    nuf.x[] = x*mu;
  scalar source[], theta[];
  foreach() {
    source[] = -dpdz*cos(omega*(t + dt/2))*x;
    theta[] = rho*x;
  }
  diffusion (u, dt, nuf, r = source, theta = theta);
}

/**
Making a movie using Gnuplot and ffmpeg.
 */

FILE * gnuplotPipe;
event init (t = 0){
  double alpha = Rad*sqrt(omega*rho/mu);
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [0: %g]\n"
	  "set size square\n"
	  "set yr [-1.7: 1.5]\n"
	  "set key bottom left\n"
	  "set grid\n"
	  "set title 'Womersley = %g'\n"
	  "set xlabel 'r [m]'\n"
	  "set ylabel 'u [m/s]'\n", Rad, alpha);
}

int frame = 0;
event movie (i += 1) {
  double alpha = Rad*sqrt(omega*rho/mu);
  double complex prefactor = I*dpdz/(rho*omega); 
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot \
          '-' w l lw 5 t 'Witzig (1914)', \
          '-' pt 5 t 'Numerical'\n");
  foreach()
    fprintf(gnuplotPipe, "%g %g\n",x, WOMERSLEY);
  fprintf(gnuplotPipe, "%g  0\n", Rad);
  fprintf(gnuplotPipe, "e\n");
  foreach()
    fprintf(gnuplotPipe, "%g %g\n",x, u[]);
  fprintf(gnuplotPipe, "e\n");
  fflush (gnuplotPipe);
  frame++;
}

event stop (t = 4.*pi/omega, last) {
  pclose (gnuplotPipe);
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  return 1;
}

/**
## References

~~~bib
@book{witzig1914,
  title={{\"U}ber erzwungene Wellenbewegungen z{\"a}her, inkompressibler Fl{\"u}ssigkeiten in elastischen R{\"o}hren},
  author={Witzig, Konrad},
  year={1914},
  publisher={Universitat Bern.}
  }
~~~
 */
