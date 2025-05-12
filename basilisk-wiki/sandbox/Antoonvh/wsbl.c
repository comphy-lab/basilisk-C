/**
# An embedded boundary in an SBL

This is a variant of [this example](wall.c)
 
![Movie](wsbl/mov.mp4)
 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "PointTriangle.h"


scalar b[], * tracers = {b};

double Re = 300;

b[bottom] = neumann(0.1*exp(-sq(x))*Re); // A local heat source
b[top] = dirichlet (y);

u.t[bottom] = dirichlet (0);
u.t[embed] = dirichlet (0);
u.n[embed] = dirichlet (0);

face vector av[], muc[];

double xw = -5, W = 0.5, H = 3; //Wall position, width, and height,

int maxlevel = 8;

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2;
  mu = muc;
  a = av;
  N = 128;
  run();
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar*){muc});
}

event init (t = 0) {
  vertex scalar phi[];
  coord wall1 = {xw, Y0}, wall2 = {xw, Y0 + H - W}, temp1[1];
  double temp2[1];
  foreach_vertex() {
    coord cc = {x, y};
    phi[] = sqrt(PointSegmentDistance (&cc, &wall1, &wall2,
				       temp1, temp2));
  }
  fractions (phi, cs, fs, W);
  foreach() {
    b[] = y;
    u.x[] = 1;
  }
  boundary (all);
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++) 
  diffusion (b, dt, mu);

event adapt (i++) 
  adapt_wavelet ({cs, b, u}, (double[]){0.01, 0.1, 0.05, 0.05},
		 maxlevel, 5);

event mov (t += 0.5) {
  draw_vof ("cs", "fs", filled = -1, fc = {0.8, 0.8, 0.8});
  draw_vof ("cs", "fs");
  squares ("b", linear = true, min = -0.1, max = 4.1); 
  mirror ({0,-1}) {
    cells();
  }
  save ("mov.mp4");
}

event stop (t = 100);