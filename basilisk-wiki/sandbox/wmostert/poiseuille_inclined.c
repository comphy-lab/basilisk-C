/**
   2D Poiseuille flow with a prescribed velocity in an angled channel. Problem file courtesy of J. Eggers. */

#include "navier-stokes/centered.h"

/**
   Define centre velocity U, domain sice L, channel half-length Lc2, channel half-width hw, inclination angle alpha, running time tfinal, and resolution level LEVEL. */
   
#define U (1.) 
#define L (8.) 
#define Lc2 (3.) 
#define hw (Lc2/8.0)
#define alpha (0.0) 
#define tfinal (5.) 
#define LEVEL (7)

/**
   Define velocities and coordinates in transformed (inclined) frame. */

#define ue(y) (U*(sq(hw)-sq(y))) 
#define yt(y) ((y-Lc2*sin(alpha))*cos(alpha)) // transversal coordinate in side channel
#define uc(y) (U*(sq(hw)-sq(yt(y)))/sq(hw)) // channel profile

int main() { // this is the actual program
  L0 = L;  // L0 is pre-defined
  /**
     Set the origin at the centre. */
  origin (-L0/2, -L0/2.);
  init_grid(1 << LEVEL);
  /**
     Viscosity is set to unity. */
  const face vector muc[] = {1.,1.}; 
  mu = muc;
  run(); 
}

/**
   Now we define the problem geometry through the use of boundary conditions and mask (the latter assigned on initialization. Although not strictly necessary for the problem, we would like user-defined left and right boundary conditions. (We could similarly use such user-defined BCs for the top and bottom, but do not do so in this example.) */

bid plug_right;
bid plug_left;

u.n[plug_left]  = dirichlet(uc(y));
p[plug_left]    = neumann(0.);
pf[plug_left]   = neumann(0.);

u.n[plug_right] = neumann(0.);
p[plug_right]   = dirichlet(0.);
pf[plug_right]  = dirichlet(0.);

u.t[top]   = dirichlet(0.);  // no slip
u.t[bottom]  = dirichlet(0.);  // no slip

/**
   Initial condition setting up the geometry, and a quiescent initial flow. */
event init (t = 0) {
  mask ((y > hw/cos(alpha) + x*tan(alpha)) ? top :  
	(y < -hw/cos(alpha) + x*tan(alpha)) ? bottom :  
	(x > Lc2*cos(alpha) - tan(alpha)*(y-Lc2*sin(alpha))) ? plug_right :
	(x < -Lc2*cos(alpha) - tan(alpha)*(y+Lc2*sin(alpha))) ? plug_left :
	none); // 

  foreach()
    u.x[] = 0.;
}


/**
   Run the output on every timestep. */
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
   Produce outputs. */
event movies (i += 4; t <= tfinal) {
  static FILE *fp = fopen("ux.ppm", "w");
  static FILE *fp1 = fopen("uy.ppm", "w");

  output_ppm (u.x, fp, n=1024, box= {{-L/2.,-L/2.},{L/2.,L/2.}},
	      linear = true, spread = 2);  
  output_ppm (u.y, fp1, n=1024, box= {{-L/2.,-L/2.},{L/2.,L/2.}},
	      linear = true, spread = 2);  
}

/**
   Adapt on the velocity. */
event adapt (i++) {
  adapt_wavelet ({u}, (double[]){3e-2,3e-2}, LEVEL, LEVEL-2);
}

event end (t=tfinal){
  fprintf (ferr, "Complete!\n");
  dump ("end");
}