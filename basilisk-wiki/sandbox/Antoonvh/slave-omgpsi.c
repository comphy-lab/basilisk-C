/**
# Dual grid $\omega-\psi$ solver

This is used in combination with [master-omgpsi.h]().

The slave is concerned with finding the streamfunction ($\psi$) */

#include "poisson.h"
#include "utils.h"
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
  return interpolate (s, xp, yp, zp, linear);
}
  
/**
The rest is for a function calling the poisson solver.
 */
scalar psi[];

psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);

int slave_init () {
  _init_solver();
  L0 = 10;
  X0 = Y0 = -L0/2.;
  init_grid (N);
  return 0;
}

void slave_solve_psi () {
  if (!grid)
    slave_init();
  scalar omega[];
  copy_field_master (omega);
  poisson (psi, omega); // Fixme : omega does not really need to be stored.
  adapt_wavelet ({psi}, (double[]){0.01}, 8);
}

//output function
void slave_level() {
  if (!grid)
    slave_solve_psi();
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "slave-level.mp4", min = 2, max = 8, n = 300);
}


void slave_free() {
  free_grid();
}
