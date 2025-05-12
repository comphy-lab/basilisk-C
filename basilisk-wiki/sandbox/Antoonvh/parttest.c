/**
# Particle advection test

On this page we test particle advection for a non-trivial, steady flow
field. We use the famous Lamb-Chaplygin vortex dipole soliton. It is
well known that this structure is characterized by an enclosed stream
line around the dipole's atmosphere. You may inspect its structure
[here](isoline.c).
*/
vector u[];
#include "view.h"
#include "poisson.h"
#define BVIEW 1
#include "particles.h"

scalar omega[];

int main() {
  L0 = 15;
  X0 = Y0 = -L0/2;
  DT = HUGE;
  foreach_dimension()
    periodic (left);
  init_grid (256);
  run();
  P_RK3 = true;
  run();
}

event init (t = 0) { //Init the dipole
  scalar psi[];
  double k = 3.83170597;
  TOLERANCE = 1e-6;
  /**
     The dipolar vortex is initialized on a grid with a maximum
     resolution that corresponds to $8192 \times 8192$ equidistant
     cells such that the error in the velocity estimation is not more
     than 0.25% of the free stream velocity ($U$).
  */
  int maxlevel = 13;
  while (depth() != maxlevel) {
    adapt_wavelet ({u.x, u.y}, (double[]){0.0025, 0.0025}, maxlevel);
    foreach() {
      double r = sqrt(sq(x) + sq(y));
      double s = x/r;
      omega[] =  ((r<1)*((-2*j1(k*r)*s/(k*j0(k)))))*sq(k);
    }
    poisson (psi, omega);
    foreach() {
      u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
      u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta) + 1.;
    }
  }
  DT = 0.01;
  dt = dtnext (DT);
  init_particles_2D_square_grid (40, 0, 0, 5.);// Add 40 x 40 = 1600 particles
}

event set_dtmax (i++)
  dt = dtnext (DT);
/**
##Output
   
We render a movie. 
*/

event movie (t += 0.1; t <= L0*4.) {
  view(width = 500, height = 500);
  squares("omega", min = -10, max = 10, map = cool_warm);
  scatter(loc = loc, s = 15);
  if (!P_RK3)
    save ("movie.mp4");
  else
    save ("movie_rk3.mp4");
    
}

/**
   ![Particle advection RK2](parttest/movie.mp4)
  
   ![Particle advection RK3](parttest/movie_rk3.mp4)
   
   Well done `particles.h`!
*/

event log_tracers (i += 5) {
  int ni = 0;
  char fname[99];
  if (!P_RK3)
    sprintf (fname, "particles");
  else
    sprintf (fname, "particles_rk3");
  static FILE * fp1 = fopen (fname, "w");
  for (int j = 0; j < n_part; j++) {
    if (sq(loc[j].x) + sq(loc[j].y) <= 1.)
      ni++;
  }
  fprintf(fp1, "%g\t%d\n", t, ni);
}

/**
  
   ~~~gnuplot Particles in the atmosphere stay there
   set xlabel 'time'
   set ylabel 'number of tracers'
   set key outside
   plot 'particles' u 1:2 w l lw 3 t 'RK2',\
        'particles_rk3' u 1:2 w l lw 3 t 'RK3'
     ~~~
*/
