/** # Feedback Mechanism in a Jet Impinging a Flat Plate

This case study investigates the hydrodynamic feedback loop (edge-tone type mechanism) 
generated when an unstable fluid jet encounters a perpendicular flat plate. 
The simulation focuses on the self-sustained oscillations and vortex shedding 
occurring between the nozzle exit and the solid obstacle.


The development of this feedback loop highly depends on the Reynolds number ($Re$) 
and the coflow velocity ratio ($c_f$). These parameters govern the growth rate of 
the shear-layer instabilities. If inappropriate values are selected—such as a 
sub-critical Reynolds number or an excessive coflow that dampens the jet—the periodic 
vortex shedding will fail to establish, and the feedback phenomenon will not be observed.
*/

  



#include "navier-stokes/centered.h"
#include "utils.h"

bid inner;

double Reynolds = 150; 
int maxlevel = 8;
double b = 1.; // jet diameter
double U0 = 1.; // for top hat
double W=10.; // edge-distance

scalar omega[];
scalar base[];
scalar fluct[];
scalar un[];
face vector muv[];
double cf=0.1; 
double du;
double Tend = 1500;

/**You can pass command-line arguments to the program to easily run parametric studies and test different configurations.

*/

int main(int argc, char * argv[])
{
  if (argc > 4) {
    maxlevel = atof (argv[1]);
    Reynolds = atof (argv[2]);
    cf = atof (argv[3]);
    Tend = atof (argv[4]);
  }
  
  /**The interface includes interactive display controls to dynamically monitor and adjust key parameters during the simulation

*/
  display_control(Reynolds,0,1000);
  display_control(cf,0,1);
  display_control(maxlevel,0,15);


  
  /**The domain size is set to 32 and is vertically centered at the inlet.*/


  L0 = 32; 
  origin (0., -L0/2.);
  N = 64;
  mu = muv;
  init_grid (N);
  run();
}


/**The simulation can start from scratch or restore a previous snapshot. 
To accelerate convergence toward the feedback regime, a pre-computed 
baseflow checkpoint can be loaded (e.g., from `sandbox/Llourenco/feedbackmechanism/baseflow.c`). 
.
*/
event init (i = 0)
{
  if (!restore ("restart")) {
    //    fprintf (stderr, "Failed to restore\n");
    //   return 1;
  }
 /** We enforce a high resolution along the jet centerline and apply a geometric mask to define the solid flat plate located at $x > W$. */
  refine(level < 9  && (fabs(y) < Delta));
  mask (fabs(y) <= L0/(1 << 9) && x > W ? inner : none);
  /**Once this steady state is restored, a localized perturbation is applied 
to trigger the hydrodynamic instability*/
  foreach() {
    if (x < 0.1*W && fabs(y) < 0.1)
      u.y[] = 1e-1; // starting the instability
    base[] = omega[];
  }
}

/**We set a constant viscosity based on the jet width and the reference velocity
*/

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds; 
  vorticity (u, omega);
  refine (x < W/4. && fabs(y) < 1. && level < maxlevel);
}


/**We implement specific boundary conditions to generate the planar jet. 
A uniform coflow ($U_\infty = c_f$) is introduced outside the jet core ($|y| > b/2$) 
to confine the flow and ensure a proper interaction with the downstream flat plate.
*/

u.n[left]  = dirichlet(fabs(y) > b/2. ? cf : 1.);
u.t[left]  = dirichlet(0.) ;

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);

u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);

/**We define the no-slip condition on the flat plate (`inner` boundary). Note that 
the conditional statement ($x < W$) is redundant here since the mask already restricts this boundary to the region $x > W$. */

u.n[inner] = x < W ? neumann(0.) : dirichlet(0.);
u.t[inner] = x < W ? dirichlet(1.) : dirichlet(0.);

/**At each step, we log the simulation data and record the transverse velocity ($u_y$) at two virtual probes (near the nozzle and the plate) to get the time series.*/
event logfile (i++, t<=Tend) {
  fprintf (stderr, "%d %g %g %g %g\n", i, t, dt, interpolate(u.y, 0.1*W, 0.), interpolate(u.y, 0.9*W, 0.));
  /**We also update the vorticity fluctuation field ($fluct = \omega - base$).*/
  foreach()
    fluct[] = omega[] - base[];
}
/**
We produce animations of the vorticity and velocity field. */
event movie (t += 2; t <= Tend) {
  scalar m[];
  foreach()
    m[] = (fabs(y) <= L0/(1 << 9) && x > W) ? 0 : 1;
  output_ppm (omega, file="vorticity.mp4", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-3, max=3, linear=true, mask=m);
  output_ppm (u.x, file="ux.mp4", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-0.5, max=1.5, linear=true, mask=m);
}

event final_images (t = Tend) {
  scalar m[];
  foreach()
    m[] = (fabs(y) <= L0/(1 << 9) && x > W) ? 0 : 1;
  output_ppm (u.x, file="ux_final.png", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-0.5, max=1.5, linear=true, mask=m);
  output_ppm (omega, file="omega_final.png", n=512,
              box={{0,-W*0.75},{L0,W*0.75}},
              min=-3, max=3, linear=true, mask=m);
}


/**Dynamically adapts the mesh based on velocity errors to capture flow structures efficiently.
*/

event adapt (i++) 
  adapt_wavelet ({u}, (double[]){1e-3,1e-3}, maxlevel, 4);



/**
# Visualisations


![Animation of the vorticity field](timeSeries/vorticity.mp4)(loop)

![Animation of the velocity field (u.x)](timeSeries/ux.mp4)(loop)
*/