/**
# Dipole-induced flow and a sharp-edged wall

Consider the case of a dipole-induced flow near a shared-edged
wall. What will happen to the fluid in the vicinity of the sharp edge?
Lets have a look.

![The evolution of the vorticity field (colors) and the tracers](sharp-edged/mov.mp4)

Obligatory refinement movie:

![Very refined](sharp-edged/cells.mp4)
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

int maxlevel = 11;
double cse = 1e-3;
double wedge_angle = 4.*pi/180.;

double xp = 0, yp = -5;

double Re = 2500;

face vector nu[];

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

Particles edge;

/**
I added the [wikipedia page on the Lamb-Chaplygin
dipole](https://en.wikipedia.org/wiki/Lamb%E2%80%93Chaplygin_dipole)
model some time ago. Better use it!
*/

#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)
void init_lamb (double sig, double xo, double yo){
  double k = 3.83170597;
  refine (sq(x - xo) + sq(y - yo) < 5 && level < maxlevel - 2);
  refine (sq(x - xo) + sq(y - yo) < 2 && level < maxlevel);
  scalar psi[], omg[];
  foreach()  {
    psi[] = 0;
    omg[] = -sig*((RAD<1)*((-2*j1(k*RAD)*ST/(k*j0(k)))))*sq(k);
  }
  poisson (psi, omg);
  foreach(){
    u.x[] += ((psi[0,1] - psi[0,-1])/(2*Delta));
    u.y[] += -(psi[1,0] - psi[-1,0])/(2*Delta);
  }
}

int main() {
  L0 = 40;
  X0 = Y0 = -L0/2.;
  mu = nu;
  for (xp = -2; xp < 2.1 ; xp += 2)
    run();
}

event init (t = 0) {
  // Add traces outside the edge, by removing those that are in.
  edge = init_tp_circle(n = 75, l = 0.5);
  remove_particles (edge, fabs(atan2(y, x)) < wedge_angle/2);
  // Iteratively refine for the sharp-edged geometry
  do {
    vertex scalar phi[];
    foreach_vertex() 
      phi[] = fabs (atan2(y, x));
    fractions (phi, cs, fs, wedge_angle/2.);
  } while (adapt_wavelet ({cs}, (double[]){cse}, maxlevel).nf > 5);
  // Add dipole
  init_lamb (1, xp, yp);
}

event properties (i++) {
  foreach_face() 
    nu.x[] = fm.x[]/Re;
}

double tend = 12;
event movies (t += 0.1; t <= tend) {
  scalar omega[];
  char text[99];
  sprintf (text, "Offset: x/R = %g", xp);
  view (fov = 8, ty = 0.05);
  vorticity (u, omega);
  squares ("omega", min = -10, max = 10, map = blue_white_red);
  draw_string (text, 1, 30, lw = 3);
  scatter (edge, s = 2);
  draw_vof ("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
  draw_vof ("cs", "fs", lw = 2);
  save ("mov.mp4");

  view (fov = 20 - 15*sq(sin(pi*t/(tend))));
  cells();
  save ("cells.mp4");
}

event adapt (i++) {
  adapt_wavelet ({cs, u}, (double[]){cse, 0.05, 0.05}, maxlevel);
}