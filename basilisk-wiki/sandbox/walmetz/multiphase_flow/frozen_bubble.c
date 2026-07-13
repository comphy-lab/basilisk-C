/**
This code is about a solid cylinder inside a gas bubble. The solid cylinder don't have any influence, and is only a vestige. Indeed, this is a simplified [Picella's cylinder plastron](https://basilisk.fr/sandbox/fpicella/cylinder_plastron/cylinder_plastron_online.c), that freeze the dynamic bubble into a static bubble. It's a first validation case, for which [embed.h](https://basilisk.fr/src/embed.h) works. This has'nt been studied much more since it wasn't usefull. It only tells us that embed doesn't works for the dynamic bubble because of the triple point at the gaz separation.
*/

//#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//#include "tension.h"
//#include "tavares/contact-embed.h"
#include "view.h"


#define R0 0.5 // solid cylinder radius
#define R1 0.55 // plastron cylinder radius
#define xc 0.
#define yc 0. 
#define T 150.

int LEVEL = 10;

double Reynolds;

int main() {
  size (32.);

  /**
  We set the origin */

  origin (-8., -L0/2.);

  double list_mu[] = {0.01, 0.05, 0.1, 0.5, 1.};
  int size_mu = sizeof(list_mu) / sizeof(list_mu[0]);

  /** We want to analyze µ's effect */
  for (int i_mu=0 ; i_mu < size_mu ; i_mu++){

    double fact_mu = list_mu[i_mu];

    for (Reynolds=5 ; Reynolds<180; Reynolds+=5){

      init_grid (1 << LEVEL);

      mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 
      mu1 = fact_mu * mu2; // plastron viscosity 

      run();
      }
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

  fraction (f, - (sq(x - xc) + sq(y - yc) - sq(R1)));
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-3,3e-2,3e-2,1e-3}, LEVEL, LEVEL-3);
}


event logfile (i++; t <= T){

  /** We freeze the bubble*/    
  fraction (f, - (sq(x - xc) + sq(y - yc) - sq(R1)));

  /** Computation of $C_d$ and $C_l$*/
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq ((1.))*(2*R0));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((1.))*(2*R0));

  fprintf (stderr, "%3.2e %f %02d %+6.5e %+6.5e %+6.5e  %+6.5e \n",
           mu1/mu2, Reynolds, i, t, dt,
           CD, CL);
  fflush (stderr);
}



/** If we want the .mp4 (it slows the computation)*/
//event movie(i+=10,t<=T){
//  view(fov=5, tx = 0, ty = 0);
//  draw_vof("cs", "fs",filled = -1);
  //draw_vof ("f", filled = 1, fc = {1,0,0});
//        draw_vof ("f", lc = {1, 0, 0}, lw = 2);
//        squares ("u.x", linear = true);
//        // Draw grid only on upper part of flow
//        cells (lc = {0.7, 0.7, 0.7});
//  save("movie.mp4");
//}


/**
My results are compared to [Sheard et al. (JFM, 2005)](https://doi.org/10.1017%2FS0022112004002836), that compute a cylinder in a monophasic fluid flow (classical case).

It appears our code is valid for 1, since we are in the case of a classical solid cylinder in a monophasic fluid flow. Then as $µ_1/µ_2$ decreases, $C_D$ is decreasing too as we could expect. Since I did'nt found any theoric validation, we have to take a step back on these results.

Since the bubble ($µ_1$) is bigger than the cylinder, it modifies the Reynolds number and can affects the results. *I have to run a simulation for a 0.55 diameter cylinder to help to take a step back.*

<img src="https://raw.githubusercontent.com/mwalmetz/intership_dalembert/main/Cd_vs_Re_mu.png" width="900" alt="Ma figure">


*/