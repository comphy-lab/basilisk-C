/**
# Flow around a 2D circular cylinder (embed) encapsulated in a fluid plastron (vof+contact angle)
Thanks to [M. Tavares](/sandbox/tavares/README) for sharing her case for [embed+vof+contact angle](/sandbox/tavares/test_cases/droplet-cylinder-embed.c)

Here consider a fixed solid cylinder (embed) at Re = 100, encapsulated in a fluid layer (plastron).
It's a multiphase flow. We've got variable viscosity (mu1,mu2), surface tension (sigma) and contact angles.

This case has no direct application, and represents a proof of concept.
*/

//#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tavares/contact-embed.h"
#include "view.h"

#define R0 0.5 // solid cylinder radius
#define R1 0.55 // plastron cylinder radius
#define xc 0.
#define yc 0. 
#define T 100.

double theta0, volume_vof_init;
int LEVEL = 8;

double Reynolds = 100.;

int main() {
  size (16.);

  /**
  We set the origin */

  origin (-L0/6., -L0/2.);

  init_grid (1 << LEVEL);


  /**
  We use a constant viscosity. */

	mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 
	/**
	We set plastron viscosity equal to 1/100 of the main fluid one.*/
  mu1 = 0.01*mu2; // plastron viscosity 
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

	/**
	Set a constant contact angle.*/
	const scalar c[] = 10.*pi/180.; // fixed contact angle...
	contact_angle = c;

	run();
}

/**
We set the boundary conditions, so to obtain a flow around a fixed cylinder. */
u.n[left]  = dirichlet(1.0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
Must impose no-slip on embedded boundaries!*/
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
  /**
  We define the solid cylinder (EMBED) and fluid cylinder (PLASTRON) 
  interface. */
  solid (cs, fs, (sq(x - xc) + sq(y - yc) - sq(R0)));

  fraction (f, - (sq(x - xc) + sq(y - yc) - sq(R1)));
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,1e-2}, LEVEL, LEVEL-4);
}

event movie(i+=10,t<=T){
  view(fov=5, tx = 0, ty = 0);
  draw_vof("cs", "fs",filled = -1);
  //draw_vof ("f", filled = 1, fc = {1,0,0});
	draw_vof ("f", lc = {1, 0, 0}, lw = 2);
	squares ("u.x", linear = true);
	// Draw grid only on upper part of flow
	cells (lc = {0.7, 0.7, 0.7});


  save("movie.mp4");
}
/**
![Some animation of plastron dynamics around a circular cylinder.](cylinder_plastron_online/movie.mp4)(width="60%")
*/

// fixme: Comparison to theory is missing (add soon) 
//
