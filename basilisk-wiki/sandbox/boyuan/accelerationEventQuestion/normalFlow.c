/**
## Case set-up
  **The question proposed [here](http://basilisk.fr/sandbox/boyuan/accelerationEventQuestion/normalFlow.c#numerical-results) has had known reasons.** This may due to a possible bug in Basilisk. A careful investigation of the source code of [ns-centered](http://basilisk.fr/src/navier-stokes/centered.h) and [iforce](http://basilisk.fr/src/iforce.h) suggests that the acceleration term defined before in the centered NS solver is re-assigned to zero in the iforce subroutine. This is totally unreasonable and could lead to unexpected results. An effort to try to use body-force gravity term and surface tension is shown [here in my GitHub page](https://github.com/MGYBY/Basilisk_practice/tree/main/power-law/vof/normalFlow).

  Steady-uniform flow down an incline.
  
  * VOF Navier-Stokes in 2D with surface tension.
  * Gravity is modelled as body force.
  * Periodic boundary condition for left-and-right boundaries.
  * 32x32 cells in 2D domain. No adaptation.
  * x-axis aligns with the inclined plane and y-axis is perpendicular to the inclined plane. Gravity $\pmb{g}$ is decomposited to two components: $(g_x, g_y)$.
  * Initial condition: 3/4 of the vertical domain (32x24) is set up as fluid phase and 1/4 of the vertical domain (32x8) is set up as air phase; velocity field is zero $u_x=u_y=0$ initially.
  * The steady-state numerical solution is to be compared with analytical solutions.
*/
/**
The analytical solutions for such a problem are:
$$u_x = \frac{3}{2} \overline{u}(1-(1-\frac{y}{h})^2)$$
$$p=\rho_{fluid}g_y (h-y)$$
where $h$ is fluid depth, $g_y$ is the slope-normal component of gravity, $\overline{u}$ is depth-averaged velocity. The second equation suggests hydrostasy. 

This code tries to reproduce the steady-unform analytical solution. However, the numerical results are absurdly wrong. I also wrote Gerris codes for the same problem, but Gerris code could reasonably reproduce the analytical solution. The link for Gerris code for the same case could be found in the following link:
**[Gerris codes](https://mcgill-my.sharepoint.com/:f:/g/personal/boyuan_yu_mail_mcgill_ca/EnA6xFGANfhPn4fa1YI8VxYBcd9xWFGcPOabUCW5Yga18A?e=6PfJsT).**
*/

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
// alternatively, use momentum-conserving scheme
// #include "navier-stokes/conserving.h"
#include "tension.h"
// Gravity is modelled as body force instead of reduced gravity
// #include "reduced.h"
#include "view.h"
#include "tag.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/** Grid levels */

#define LEVEL 5

/** Problem-related parameters */

#define FR 1.50 // Froude number for the expected steady-state solution, but this dimensionless parameter is not used in the codes here.
#define RATIO 850.0/1.0 //density ratio, water to air
#define MURATIO 8.9e-4/17.4e-6 //dynamic viscosity ratio, water to air

// inclination angle of the channel. \sin\theta and \cos\theta
#define CHANNELSLOPE 0.015
#define CHANNELCOS (pow((1.0-pow(CHANNELSLOPE,2.0)),0.50))

#define MAXTIME 20.0 // Maximum runtime.

// expected analytical solutions for depth and depth-averaged velocity, but we start simulation with fluid at rest.
#define NORMALDEPTH 0.0020933486
#define NORMALVEL 0.214942405

double g_ = 9.81;

// square domain size
#define xextent_ (4.0/3.0*NORMALDEPTH)

/**
## Main body of the current codes
*/

/** ### Main */
int main()
{
  size (xextent_);

    rho1 = RATIO;
    rho2 = 1.0;

    mu1 = 8.9e-4;
    mu2 = mu1/MURATIO;

  // Surface tension seems not to change the solution too much, since there is very little interface curvature.
  f.sigma = 0.072;
  init_grid(1 << (LEVEL));
  // Acceleration using reduced gravity. But reduced gravity approach does not work for this case. 
//   G.y = (-CHANNELCOS)*g_;
//   G.x = (CHANNELSLOPE)*g_;

  /** Body-force gravity. This defines the acceleration vector $\pmb{a}$ in $\texttt{centered.h}$ file.*/
  const face vector gravity[] = {(CHANNELSLOPE)*g_, (-CHANNELCOS)*g_, 0.0};
  a = gravity;

  // periodic BC
  periodic (right);

  run();
}

/** ### Init event */
event init (i=0)
{
      // set 3/4 of the Ly domain with fluid and remaining 1/4 with air
      fraction (f, NORMALDEPTH-y);

      foreach() {
    	u.x[] = 0.0;
        u.y[] = 0.0;
//         p[] = y<=NORMALDEPTH ? (NORMALDEPTH-y)*(CHANNELCOS)*g_ : 0.0;
      }
      boundary ((scalar *){u});
}

/** ### Maximum run-time control and time & time-step logging */
event maxdt (t <= MAXTIME; t += 0.50);

event timingLog(i += 10) {
  fprintf (stderr, "%d %g %g \n", i, t, dt);
}

/** ### No AMR is used for now */
//-------------------ADAPTIVITY---------------------//
/*Adapt once on error in volume fraction, velocity field, and beach fraction*/
// event adapt(i++) {
//   //double uemax = 1e-5;
//
// //   double femax = 1e-3;
//   double uemax = NORMALVEL/120.0;
// //   adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
//   adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){uemax,uemax,uemax}, MAXLEVEL, MINLEVEL, 1);
// }

/**
### Dump output, Gfs file output, free-surface output*/

event snapshot (t += 1.0) {
  char nameOut[80];
  sprintf (nameOut, "dumpSnapshot-%g", t);
  dump(file=nameOut);
}

event outputGfsFiles (t += 1.0) {
    char name[80];
    sprintf(name, "out-%.3f.gfs", t);
    FILE *fp = fopen(name, "w");
    output_gfs(fp, translate = true);
    fclose (fp);
}

event outputCentVel(t += 1.0) {
  char resultname[40];
  sprintf( resultname, "centVel_%g.txt", t );
  FILE * fp = fopen(resultname, "w");
  for (double y = 0.; y < xextent_; y += xextent_/pow(2.,LEVEL))
        fprintf (fp, "%g %g %g %g \n", y, interpolate (u.x, xextent_/2, y), interpolate (u.y, xextent_/2, y), interpolate (f, xextent_/2, y));
  fclose (fp);
}

/**## Numerical Results */
/**
 <span style="color:red"> **The fluid stays stationary ($u_x(t)=u_y(t)=0$), which is apparently different from the analytical solutions and it is not reasonable at all. Why??** </span>*/
