/**
Code by F.Picella ([original](https://basilisk.fr/sandbox/fpicella/cylinder_plastron/cylinder_plastron_online.c)) modified to suppress multiphase flow and model a simple 2D cylinder.

I just want to validate it by comparing it to a typical (Cd,Re) diagram.
*/


//#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//#include "tension.h"
//#include "tavares/contact-embed.h"
#include "view.h"

#define R0 0.5 // solid cylinder radius
#define xc 0.
#define yc 0. 
#define T 100.

int LEVEL = 10;

double Reynolds;

int main() {
  size (32.);

  /**
  We set the origin */

  origin (-8, -L0/2.);

  // boucle sur les LEVEL
  for ( Reynolds = 5 ; Reynolds < 180 ; Reynolds++){


    init_grid (1 << LEVEL);

  /**
  We use a constant viscosity. */

        mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 

        run();
   }
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

}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, LEVEL, LEVEL-3);
}

event compute_forces (i++, t<=T)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
    double Cx = (Fp.x+Fmu.x)*2.;
    double Cy = (Fp.y+Fmu.y)*2.;
    fprintf (stderr, "%+3.2e %06d %+6.5e %+6.5e %+6.5e %+6.5e \n",
             Reynolds, i, t, dt, Cx, Cy);
    fflush (stderr);
}

/**
If we need to visualise what we’re doing (it slows down the computing) :
*/

//
//event movie(i+=10,t<=T){
//  view(fov=5, tx = 0, ty = 0);
//  draw_vof("cs", "fs",filled = -1);
  //draw_vof ("f", filled = 1, fc = {1,0,0});
//        draw_vof ("f", lc = {1, 0, 0}, lw = 2);
//        squares ("u.x", linear = true);
        // Draw grid only on upper part of flow
//        cells (lc = {0.7, 0.7, 0.7});


//  save("movie.mp4");
//}

/**
Comparison to theory is coming on monday or tuesday : currently running ;)
*/

                                                            