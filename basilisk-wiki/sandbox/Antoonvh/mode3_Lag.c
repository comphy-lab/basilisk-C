/**
# Adaptive Hybrid Lagrangian-Eulerian $\omega-\psi$ dual-grid solver

![Works pretty good. Well done Lagrangian advection!](http://antoonvanhooft.nl/media/mode3_Lag.mp4)
 */
#define RKORDER 3

#include "view.h"
#define VIEW 1
#include "master-omgpsi_LAG.h"
#include "scatter2.h"
#include "../prouvost/AMR_tools/amr.h"


double tol = 1e-3;



int main() {
  //  maxlevel = 9;
  periodic(left);
  periodic_x = true;
  L0 = 8.;
  origin (-pi*4./3., -2*exp(1) + 1.2); // Not centered
  N = 1 << 7;
  run();
}

#define RAD (sqrt((sq(x) + sq(y))))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M))))

double P = 1e-2, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameter


event init (t = 0) {
  AMReps = tol;
  
  slave_init (L0, (coord){X0, Y0, 0}, 8, true, tol);
 
  double betav = 1./(1 - sq(b));
  p_omg = new_particles(0);
  int nc = 999, nf = 999;
  while (nc > 100 || nf > 100) {
    scalar omega_uf[]; //unfiltered omega 
    omega[bottom] = dirichlet (0);
    omega[top] = dirichlet (0);
    printf ("%d %ld\n", maxlevel, grid->tn);
    foreach() {
      omega[] = 0;
      double rp = RAD*RADP(P,m), omg = 0;
      if (rp <= 1.)
	omg = 2;
      else if (rp > 1 && rp <= b) 
	omg = 2*betav;
      omega_uf[] = omg;
    }
    const face vector alphaf[] = {-sq(0.4/(2.*pi)), -sq(0.4/(2.*pi))};
    poisson (omega, omega_uf, alphaf, unity);
    nf = adapt_metric({omega}).nf;
    nc = adapt_metric({omega}).nc;
  }
  printf ("%d %ld\n", maxlevel, grid->tn);
  foreach() {
    double o = omega[];
    foreach_child() {
      particle pn;
      pn.x = x;
      pn.y = y;
      pn.s = interpolate(omega, x, y);
      add_particle(pn, p_omg);
    }
  }
  velocity(p_omg);
}

event mov (t += 0.1, last) {
  view (width = 800, height = 800, bg = {0,0.5, 0.5});
  scatter_color (p_omg, s = 1, map = blue_white_red);
  save ("mode3.mp4");
  save ("mode3.png");
  cells();
  squares ("omega", map = blue_white_red);
  save ("grid.mp4");
  slave_level();
}

event adapt (i++ ) {
  compute_omega(p_omg, omega);
  //adapt_wavlet({omega},(double[]){tol}, maxlevel, 5);

  adapt_metric({omega});
  adapt_number();
}
event stop (t = 10) {
  return 1;
}
