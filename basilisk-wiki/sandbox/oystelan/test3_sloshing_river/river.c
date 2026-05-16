/** Sloshing River
A testcase to investigate spurious currents caused by the gradient
limiter of thin cells with dry neighbour cells. The test case is a sloshing
river where which should not generate any vorticity/turbulence at the shorelines.
However they do when we use the original implementation, as can be seen here

![vertical velocity in a sloshing river (original monmod limiting of thin almost dry cells)](river/movie.mp4)

The reason for these currents are the gradient limiters (in this case minmod2) which for cells with
neighbours that is close to dry will result in unwanted spurious current into the domain.

The problem can easy enough be solved by introducing a "soft" cut-off in the advection scheme, 
giving much more predictable behaviour of the gradient limiters.

The result can be seen underneath.
![vertical velocity in a sloshing river (soft-cutoff of velocity)](river-new/movie.mp4)

*/


#include "grid/multigrid.h"
#if switch_solver
#include "../multilayer_stability_patch/hydro_patched.h"
#else
#include "layered/hydro.h"
#endif
#include "layered/implicit.h"


#if !NH // if NH is not defined:
scalar w, phi, omega;
#endif

/**
GPUs and bview are not compatible yet. */

#if !_GPU
# include "view.h"
#endif

#include "layered/remap.h"
#include "layered/perfs.h" 
#include <sys/stat.h>

#define T0 50.
#define L0 1.
#define LEVEL 7
FILE* fp_energy;

int main()
{
  periodic(right);
  init_grid (1 << LEVEL);
  size (1.);
  origin (-L0/2.,-L0/2.);
  G = 9.81; 
  nl = 11;
  TOLERANCE = 1E-5;
  //theta_H = 0.501;
  #if switch_solver
  h_lim = 1e-6; // soft cutoff to avoid rediculous velocities
  #endif
  run();
}

/** we define the direction of the river, and the velocity */
#define rot 0. // 0 = from left-right, pi/2 = top-bottom
#define U 0.02
#define zb(x,y) 0.05*(-cos((y*cos(rot)+x*sin(rot))*1.5*pi))

/** a ramp function to ensure a smooth start of the simulation */
double dt_max_rampup(double t, double dtinit, double rampuptime, double dtmax){
    return dtinit + dtmax*clamp(t/rampuptime,0,1.);
}
/** We initialize with water at a small tilt angle and with velocity in the x direction*/
event init (i = 0)
{
   #if !NH // if NH is not defined:
    assert (nl > 0);
    w = new scalar[nl];
    phi = new scalar[nl];
    omega = new scalar[nl];
    reset ({w, phi}, 0.);
  #endif
  CFL = 0.5;
  dt = 0.0001;
  foreach() {
    zb[] = zb(x,y);
  }
  // Set surface elevation at z = 0., and constant flow velocity
  foreach()
    foreach_layer(){
      h[] = max((-zb[]+y*sin(pi/180.))/nl,0.);
      u.x[] = h[] > 0. ? cos(rot)*U : 0.;
      u.y[] = 0.;
    }
  boundary(all);
}

/** update maximum allowed timestep */
event update_dt(i++){
    DT = dt_max_rampup(t, 0.00001, 1., 0.25);
}

event logfile (i += 10) {
  stats s = statsf (h);
#if IMPLICIT && !ML
  scalar u[];
  foreach()
    u[] = h[] > dry ? q.x[]/h[] : 0.;
  norm n = normf (u);
#else
  norm n = normf (u.x);
#endif
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %.4g %g %g %g %g %g\n", 
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt);
}

/**
And generate the movie of the free surface (this is quite
expensive). The movie is 45 seconds at 25 frames/second. 

Since GPUs and bview are not compatible yet, we replace 3D with 2D
outputs when running on GPUs. */

event movie (t += 0.05; t<= T0)
{
#if !_GPU  
  view (fov = 24, quat = {0, 0, 0, 1},
	tx = 0, ty = 0, width = 400, height = 400);
  char s[80];
  sprintf (s, "t = %.2f T0", t/T0);
  draw_string (s, size = 80);
  for (double x = -1; x <= 1; x++)
    translate (x)
      squares ("u10.y", linear = true, z = "eta", min = -0.1, max = 0.1);
  save ("movie.mp4");
#else // _GPU
  vector u = lookup_vector ("u29");
  output_ppm (u.x, file = "movie.mp4", linear = true, min = -0.15, max = 0.6, n = 512);
#endif // _GPU
}