/** Sloshing River
A testcase to investigate spurious currents caused by thin cells
when using different advection schemes

*/


#define output_vtk 1 // switch on/off vtk file output
#define solver_type 1 // 0 = original basilisk, 1 = new impl.



#include "grid/multigrid.h"
#if solver_type
#include "hydro_patched.h"
#else
#include "layered/hydro.h"
#endif
#include "layered/implicit.h"


#if !NH // if NH is not defined:
scalar w, phi, omega;
#endif

#include "layered/remap.h"
#include "layered/perfs.h"
#if output_vtk
#include "../output_vts_multilayer.h"
#endif
#include <sys/stat.h>


#define LEVEL 7
int num_omp = 12;
FILE* fp_energy;

int main()
{
  periodic(right);
  omp_set_num_threads(num_omp);
  init_grid (1 << LEVEL);
  size (1.);
  origin (-L0/2.,-L0/2.);
  G = 9.81; 
  nl = 11;
  TOLERANCE = 1E-5;
  #if solver_type
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

/**/
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
* output vts files which can be viewed in paraview
*/
#if output_vtk
//event output_domain_vtk(t += T0/4.; t <= 5*T0) {
event output_domain_vtk(t <= 100.; t += 0.05) {
    foreach_layer(){
       vorticity(u,omega);
    }
    struct stat st = { 0 };
    if (stat("vtk", &st) == -1) {
        mkdir("vtk", 0700);
    }
    static int j = 0;
    char name[100];
    sprintf(name, "%s/domain_%.6i.vts", "vtk", j++);
    //fprintf(stdout, "written to: %s\n", name);
    FILE* fp = fopen(name, "w");
    output_vts_bin_all_layers(fp,(scalar *){omega}, NULL, false, false);
    fclose(fp);
}
#endif
