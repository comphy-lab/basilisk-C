/**
# Collision of 3 droplets at the center of a periodic box

This test case aims to validate the extension of 
[the momentum conservation algorithm](/src/navier-stokes/centered.h) 
to the [multi-VOF method](no-coalescence.h). \
The idea is to use a large enough density contrast $\rho_2/\rho_1 = 1000$ 
to observe the kinetic energy transfer between phases when droplets are colliding. \
As it is a 2D test case, critical Weber is taken as 
$We = \frac{\rho_1*u_d^2*d}{\sigma} = 0.1$ to prevent large deformations.

*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase-multi.h"
#include "conserving-multi.h" // Remove this to go back to standard scheme
#include "tension.h"
#include "view.h"
#include "output.h"
#include "no-coalescence.h"
#include "tag.h"
#include "../misc/fractions_init.h"

/** Define mesh and time */
#define LEVEL 7
#define TMAX 1.
#define TMOVIES 0.01

/** Define physical parameters for the problem */
#define Ls 1.0 // Domain length for periodic box
#define u_d 1.0 // Droplet velocity
#define D1 (0.2*Ls) // Droplet diameter
#define rho_f 1.0 // Gas phase density
#define mu_f 0.0 // Gas phase dynamic viscosity
#define rho_d 1000 // Liquid phase density
#define mu_d 0.0 // Liquid phase dynamic viscosity
#define WeCrit 0.1 // Critical Weber
#define sig (rho_d*u_d*u_d*D1/WeCrit) // Surface tension

int main(int argc, char** argv)
{  
  L0 = Ls [0];
  DT = HUGE [0];
  init_grid(1<<(LEVEL));
  origin(-Ls/2.,-Ls/2.);
  rho1      = rho_d;
  mu1       = mu_d;
  rho2      = rho_f;
  mu2       = mu_f;
  f.sigma   = sig;
  TOLERANCE = 1e-5 [*];
  foreach_dimension()
    periodic(right);
  run();
}

/**Initializes the domain with three droplets at equidistance 
   and convergent velocities*/
event init (t=0)
{
  int ndroplets = 3;
  coord Positions[ndroplets];
  double Radii[ndroplets];
  for (int j = 0; j < ndroplets; j++){
    Radii[j] = 0.5*D1; 
  }
  Positions[0].x = 0;
  Positions[0].y = 0.5*sqrt(3)*D1;
  Positions[1].x = -D1;
  Positions[1].y = -0.5*sqrt(3)*D1;
  Positions[2].x = D1;
  Positions[2].y = -0.5*sqrt(3)*D1;
  foreach() {
    coord center = {x,y};
    f[] = refine_frac(center, Positions, Radii, ndroplets, Delta, 1);
    if (sq(x-Positions[0].x)+sq(y-Positions[0].y) < sq(0.5*D1)) {
      u.x[] = 0;
      u.y[] = -u_d;
    }
    if (sq(x-Positions[1].x)+sq(y-Positions[1].y) < sq(0.5*D1)) {
      u.x[] = cos(pi/6)*u_d;
      u.y[] = sin(pi/6)*u_d;
    }   
    if (sq(x-Positions[2].x)+sq(y-Positions[2].y) < sq(0.5*D1)) {
      u.x[] = -cos(pi/6)*u_d;
      u.y[] = sin(pi/6)*u_d;
    }
  }
}

/** Outputs a video of the vorticity with interface fragments.
*/

event movie (t = 0,t += TMOVIES; t <= TMAX)
{
  scalar omega[];
  vorticity (u, omega);
  squares(color = "omega");
  for (scalar s in interfaces){
    draw_vof(s.name, map = cool_warm, min = 0, max = 3, lw = 2.);
  }
  save ("vorticity.mp4");
}

/** Track the evolution of kinetic energy and surface energy.
*/

event logfile (i += 10) {
  double vol = 0., rke = 0., se = 0.;
  foreach (reduction(+:vol) reduction(+:rke) reduction(+:se)) {
    vol += dv();
    foreach_dimension() 
      rke += rho(f[])*sq(u.x[])*dv();
    for (scalar s in interfaces)
      if (s[] > EPS && s[] < 1 - EPS) {
        coord  normal = mycs(point, s), parea;
        double alpha = plane_alpha(s[], normal);
        double dS = pow(Delta, dimension - 1) * plane_area_center(normal, alpha, &parea);
        se += sig*dS;
      }
  }
  rke /= 2.*vol;
  se /= vol;
  fprintf (stderr,"%g %g %g\n",t, rke, se);
  fflush (stderr);
}

/**
![Collision of 3 droplets. The color field is the
 vorticity of the flow.](triCollision/vorticity.mp4)
 */