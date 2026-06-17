/**
![Stokes' Paradox has its own Wikipedia entry, [see](https://en.wikipedia.org/wiki/Stokes%27_paradox) ](https://upload.wikimedia.org/wikipedia/commons/a/ad/Ggstokes.jpg)
   
# Stokes' Paradox

We want to model the potential and viscous (Stokes) flow around a
cylinder.

The potential-flow solution is readily found.

![Velocity magnitude field (color) and stream lines of the potential
 flow around a cylinder](Stokes_paradox/potential.png)

We compute the viscous evolution of the flow, hoping to find a steady state.

![Stokes' Paradox](Stokes_paradox/stokes.mp4)
 
For a steady, creeping flow, around a cylinder, the effect of the
viscous friction extends to infinity! It therefore quickly becomes
incompatible with our box boundary conditions. Infact, no steady
(non-trivial) solution exists for this case. This is the "Stokes
Paradox".
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

face vector nu[];

void show_flow(char * fname) {
  scalar psi[], omega[];
  vertex scalar psiv[];
  vorticity(u, omega);
  
  foreach() {
    psi[] = 0;
    if (cs[] > 0 && cs[] < 1) {
      coord n, p;
      embed_geometry (point, &p, &n);
      omega[] = embed_vorticity (point, u, p, n);
      }
      if (cs[] == 0)
	omega[] = 0;
    }
  psi[top] = dirichlet (Y0 + L0);
  psi[left] = dirichlet(y);
  psi[right] = dirichlet (y);
  psi[bottom] = dirichlet(Y0);
  psi[embed] = dirichlet (0);
  poisson (psi, omega, alpha = fs, tolerance = 1e-1);
  foreach_vertex() 
    psiv[] = (psi[0,0] + psi[-1,0] + psi[-1,-1] + psi[0,-1])/4.;
  scalar Up[];
  foreach() {
    Up[] = 0;
    foreach_dimension() 
      Up[] += sq(u.x[]);
    Up[] = Up[] > 0 ? sqrt(Up[]): 0;
  }
  view (tx = 0.1);
  draw_vof("cs", "fs", filled = -1, fc = {0.5, 0.5, 0.5});
  draw_vof("cs", "fs", lw = 2);

  squares ("Up", min = 0.8, max = 1.2, cbar = true, label = "U/Ui");
  
  translate (z = 0.02) 
    for (double val = -5; val <= 5; val += 1) 
      isoline ("psiv", val, lw = 2);
  
  save (fname);
}

int main() {
  L0 = 25;
  X0 = Y0 = -L0/2.;
  mu = nu;
  N = 256;
  run();
}

event init (t = 0) {
  /**
We start with finding the geometry;
  */
  
  refine (sq(x) + sq(y) < 2 && level < 9);
  refine (sq(x) + sq(y) < 2 && level < 10);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x) + sq(y)) - 1;
  fractions(phi, cs, fs);
  /**
Then the potential-flow solution;
   */
  scalar U[], zeros[];
  U[left] = dirichlet(L0);
  U[right] = dirichlet(0);
  foreach() 
    U[] = zeros[] = 0.;
  poisson (U, zeros, alpha = fs);
  foreach() 
    foreach_dimension() {
    u.x[] = (U[1] - U[-1])/(2*Delta);
  }
  show_flow("potential.png");
  
  /**
Finally, we attempt to find a Stokes-flow solution
   */
  stokes = true;
  DT = 0.05;
  u.n[left] = dirichlet(1);
  u.n[right] = neumann(0);
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  p[left] = neumann (0);
  p[right] = dirichlet (0);
}

event properties (i++) {
  foreach_face() 
    nu.x[] = fs.x[];
  DT *= 1.1;
}

event movie (t += 0.1) {
  show_flow ("stokes.mp4");
}

event stop (t = 10) {
  return 1;
}
