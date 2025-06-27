/**
# Dual grid $\omega-\psi$ solver

This is used in combination with [master-omgpsi.h]().

The slave is concerned with finding the streamfunction ($\psi$) */

#include "poisson.h"
#include "utils.h"
#include "higher-order.h"
#include "../prouvost/AMR_tools/amr.h"
/** We need to copy data from the master and vice versa.*/

extern double master_value (const char * name, double xp = 0, double yp = 0,
			    int i, int j, int lev);

void copy_field_master (scalar s){
  foreach() {
    s[] = master_value (s.name, x, y, point.i, point.j, level);
  }
}

double slave_interpolate(const char * name, double xp = 0, double yp = 0,
			 double zp = 0, bool linear = false)
{
  if (!grid) {
    fprintf (stderr, "slave_interpolate: error: no grid! this may be a "
	     "master/slave synchronization issue\n");
    exit (1);
  }
  scalar s = lookup_field (name);
  if (s.i < 0) {
    fprintf (stderr, "slave_interpolate: error: unknown field '%s'\n", name);
    exit (1);
  }
  //printf ("%s\n", s.name);
  double v = interpolate (s, xp, yp, zp, linear);
  if (v == HUGE) 
    return 0;
  return v;
}

/**
The rest is for a function calling the poisson solver.
*/
scalar psi[]; //Reuse psi from previous step
face vector uf[]; //Always have uf on hand for for the master

psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);


bool adap = false;

int slave_init (double L = 8, coord O = {0,0,0}, int maxlvl = 7, bool adapt = false, double eps = 1e-2) {
  _init_solver();
  AMReps = eps;
  L0 = L;
  origin (O.x, O.y); 
  //maxlevel = maxlvl;
  if (grid)
    free_grid();
  init_grid (1 << maxlvl);
  psi.prolongation = refine_3rd;
  TOLERANCE = 1e-4;
  adap = adapt;
  return 0;
}

void compute_flow() {
  boundary({uf});
  foreach_face() {
    uf.x[] = (psi[0] - psi[-1])/Delta;
  }
}

void slave_solve_psi () {
  scalar omega[];
  copy_field_master (omega);
  poisson (psi, omega); 
  compute_flow();
  if (adap == true) {
    vector u[];
    foreach() {
      foreach_dimension()
	u.x[] = (uf.x[1] + uf.x[0])/2.;
    }
    multigrid_restriction((scalar*){u});
    boundary((scalar*){u});
    adapt_metric({u.x, u.y});
  }
}

//output function
#include "view.h"

void slave_level() {
  if (grid) {
    view (width = 800, height = 800);
    
    vertex scalar psiv[];
    foreach_vertex()
      psiv[] = (psi[] + psi[-1] + psi[-1,-1] + psi[0,-1])/4;
    for (double pv = -1; pv <= 1; pv += 0.1) {
      if (fabs(pv) > 1e-3)
	isoline("psiv", pv, lw = 2, lc = (float[]){0,0.5, 0.5});
    }
    cells();
    save ("slave.mp4");
    save ("slave.png");
  }
}


void slave_free() {
  free_grid();
}
