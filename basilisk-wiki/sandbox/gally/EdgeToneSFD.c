/** 
# Edge tone's baseflow calculated with the Selective Frequency Method.
*/

#include "navier-stokes/centered.h"
#include "utils.h"
#include "SFD.h"

bid inner;

double Reynolds = 150; 
int maxlevel = 9;
double b = 1.;  // jet diameter
double U0 = 1.; // for top hat
double W=10.;   // edge-distance

/**
We set the SFD parameters. */
double freq_SFD = 0.0373;
bool SFD_toggle;
event adapt_toggle (i++)
  SFD_toggle = (t >= 0.);

scalar omega[];
scalar base[];
scalar fluct[];
scalar un[];
face vector muv[];
double cf = 0.;
double du;
double Tend = 800;

int main(int argc, char * argv[])
{
  if (argc > 4) {
    maxlevel = atof (argv[1]);
    Reynolds = atof (argv[2]);
    cf = atof (argv[3]);
    Tend = atof (argv[4]);
  }
  display_control(Reynolds,0,1000);
  display_control(cf,0,1);
  display_control(maxlevel,0,15);
  L0 = 32; 
  origin (0., -L0/2.);
  N = 64;
  init_grid (N);
  mu = muv;
  
  run();
}

event init (i = 0) {  
  
  refine(level < 9  && (fabs(y) < Delta));
  mask (fabs(y) <= L0/(1 << 9) && x > W ? inner : none);
  foreach() {
    if (x < 0.1*W && fabs(y) < 0.1)
      u.y[] = 1e-1; // starting the instability
    base[] = omega[];
  }

/**
We set initial conditions. */  
  foreach()
    u.x[] = 0.;
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds; 
  vorticity (u, omega);
  refine (x < W/4. && fabs(y) < 1. && level < maxlevel);
}

/**
We set boundary conditions. */

u.n[left]  = dirichlet(fabs(y) > b/2. ? cf : 1.);
u.t[left]  = dirichlet(0.) ;

p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);

u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);

u.n[inner] = x < W ? neumann(0.) : dirichlet(0.);
u.t[inner] = x < W ? dirichlet(1.) : dirichlet(0.);


event logfile (i++, t<=Tend) {
  fprintf (stderr, "%d %g %g %g %g\n", i, t, dt, interpolate(u.y, 0.1*W, 0.), interpolate(u.y, 0.9*W, 0.));
  foreach()
    fluct[] = omega[] - base[];
}

event movie (t += 2; t <= Tend) {
  
  scalar m[];
  foreach()
    m[] = (fabs(y) <= L0/(1 << 9) && x > W) ? 0 : 1;
  
  output_ppm (omega, file="vorticity_baseflow.mp4", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-3, max=3, linear=true, mask=m);
  output_ppm (u.x, file="ux_baseflow.mp4", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-0.5, max=1.5, linear=true, mask=m);
}

event adapt (i++) 
  adapt_wavelet ({u}, (double[]){1e-3,1e-3}, maxlevel, 4);

/**

![Animation of the vorticity field](EdgeToneSFD/vorticity_baseflow.mp4)(loop)

![Animation of the horizontal velocity field](EdgeToneSFD/ux_baseflow.mp4)(loop)

*/