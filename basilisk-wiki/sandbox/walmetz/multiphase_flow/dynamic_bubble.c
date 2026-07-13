/**
# Flow around a 2D circular cylinder (embed) encapsulated in a fluid plastron (vof+contact angle)

**Taken on F.Picella's sandbox but slightly modified** (see [the original one](https://basilisk.fr/sandbox/fpicella/cylinder_plastron/cylinder_plastron_online.c)

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
//#include "view.h"

#define R0 0.5 // solid cylinder radius
#define R1 0.55 // plastron cylinder radius
#define xc 0.
#define yc 0. 
#define T 200.
double theta0, volume_vof_init;
int LEVEL = 11;

double Reynolds = 40.;

int main() {
  size (32.);

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


event logfile (i++; t <= T){

  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq ((1.))*(2*R0));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((1.))*(2*R0));

  fprintf (stderr, "%f %02d %+6.5e %+6.5e %+6.5e  %+6.5e \n",
           Reynolds, i, t, dt,
           CD, CL);
  fflush (stderr);
}


//event movie(i+=10,t<=T){
//  view(fov=5, tx = 0, ty = 0);
//  draw_vof("cs", "fs",filled = -1);
//  //draw_vof ("f", filled = 1, fc = {1,0,0});
//        draw_vof ("f", lc = {1, 0, 0}, lw = 2);
//        squares ("u.x", linear = true);
//        // Draw grid only on upper part of flow
//        cells (lc = {0.7, 0.7, 0.7});
//
//
//  save("movie.mp4");
//}

/**
Although the code looks good and run, the method used for drag isn't working well on our simulation. Indeed, the drag is negative (physically impossible). An other method has been tried ([such as this one, p.33](https://www.phys.ens.fr/sites/default/files/2024-05/M_Rabaud_NotesCoursL3_FIP.pdf)) but drag is still negative...

Since my intership isn't long enough, I let it unresolved...

<img src="https://raw.githubusercontent.com/mwalmetz/intership_dalembert/main/cd_Re40.png" width="700" alt="Ma figure">
*/
