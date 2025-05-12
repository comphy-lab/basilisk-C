/**
# Dual grid $\omega-\psi$ solver

This is used in combination with [slave-omgpsi.c]().

The master is concerned with time integration of the vorticity ($\omega$) field.
*/
#pragma autolink slave-omgpsi.o
#include "utils.h"

double master_value (const char * name, double xp, double yp, int i, int j, int lev)
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
  Point point = {.i = i, .j = j, .level = lev}; // foreach_point() is not an alternative...
  Point p1 = locate (xp, yp);
  if (p1.level >= lev)
    return s[]; // more accurate and faster
  else // ?
    interpolate (s, xp, yp, 0, true);
}

extern double slave_interpolate(const char * name, double xp = 0, double yp = 0, double zp = 0, bool linear = false);

void copy_field_slave (scalar s){
  foreach() {
    double zp = 0;
    if (dimension == 3)
      zp = z;
    s[] = slave_interpolate (s.name, x, y, zp, true);
  }
}

extern void slave_solve_psi();
extern void slave_level();
extern int slave_init();
extern void slave_free();

#include "advection.h"
scalar omega[], * tracers = {omega}, psi[];

psi[right]  = dirichlet(0);
psi[left]   = dirichlet(0);
psi[top]    = dirichlet(0);
psi[bottom] = dirichlet(0);

event velocity (i++) {
  slave_solve_psi();
  copy_field_slave (psi);
  struct { double x, y; } f = {-1.[0],1.[0]};
  foreach_face()
    u.x[] = f.x*(psi[0,1] + psi[-1,1] - psi[0,-1] - psi[-1,-1])/(4.*Delta);
}

/**
We let the Poisson-solver grid know when the simulation is finished.
*/

event stop (t = end) {
  slave_free();
}
