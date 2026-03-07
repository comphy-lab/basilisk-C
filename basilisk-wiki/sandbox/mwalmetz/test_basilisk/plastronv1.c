/**

<span style="font-size: 200%;">**Code copié sur fpicella, annoté afin de mieux le comprendre**</span>
 
Pour le vrai code, cf. : [code cylindre plastron](/sandbox/fpicella/cylinder_plastron/cylinder_plastron_online.c)

# Flow around a 2D circular cylinder (embed) encapsulated in a fluid plastron (vof+contact angle)
Thanks to [M. Tavares](/sandbox/tavares/README) for sharing her case for [embed+vof+contact angle](/sandbox/tavares/test_cases/droplet-cylinder-embed.c)

Here consider a fixed solid cylinder (embed) at Re = 100, encapsulated in a fluid layer (plastron).
It's a multiphase flow. We've got variable viscosity (mu1,mu2), surface tension (sigma) and contact angles.

This case has no direct application, and represents a proof of concept.
*/

//On importe beaucoup de programmes exterieurs, a voir si je dois me plonger dedans, demander a Picella
//#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tavares/contact-embed.h"
#include "view.h"

//define dit au compilateur de remplacer R0 par 0.5 au dernier moment. Permet d'avoir un code plus lisible
#define R0 0.5 // solid cylinder radius
#define R1 0.55 // plastron cylinder radius
#define xc 0.
#define yc 0. 
#define T 100. //temps de simulation

double theta0, volume_vof_init; //theta0:   volume_vof_init:
int LEVEL = 8; //LEVEL:

double Reynolds = 100.;
// ordre de grandeur : si Re=100, qu'on a de l'eau et D0=1cm, V=1cm/s

int main() {
  size (16.); //equivalent a L0=16

  /**
  We set the origin */

  origin (-L0/6., -L0/2.); //on decale simplement l'origine

  init_grid (1 << LEVEL); //decalage binaire. Globalement forme juste une grille cartesienne de 2^LEVEL x 2^LEVEL mailles


  /**
  We use a constant viscosity. */

	mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity, autour de 10^-3 pour l'eau
	/**
	We set plastron viscosity equal to 1/100 of the main fluid one.*/
  mu1 = 0.01*mu2; // plastron viscosity, autour de 2.10^-5 pour l'air donc realiste
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