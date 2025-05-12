/*
From Antoon's Sandbox
*/

#ifndef ADD_PART_MEM
#define ADD_PART_MEM coord u; coord u2; long unsigned int tag;
#endif

#define PA_INP particle * pp; double dt;
#define INP    (PA_inp){p(), &p(), 2.*dt}

#include "inertial-particles.h"

				    /**
For two-way couling, the forces on the particles are stored in a
field.
				    */

extern vector u;        //Fluid medium flow
extern face vector mu;  //Fluid medium dynamic viscosity
extern scalar rho;      //Fluid medium density
coord G;                //Gravity acceleration vector (should it be external) ?

coord p_acc (PA_inp inp) {
  particle pa = inp.p;
  double xp = pa.x; double yp = pa.y; double zp = pa.z;  

  double muc = is_constant (mu.x) ? constant(mu.x) : interpolate (mu.x, xp, yp, zp);
  double rp = pa.u2.y;
  double rhop = pa.u2.x;
  double tau  = rp ? (muc > 0 ? rhop * sq(2*rp)/(18*muc) : HUGE) : pa.u2.z ? pa.u2.z : HUGE;
  double itau = tau > 0 ? 1/tau : HUGE;
  double idt  = inp.dt > 0 ? 1./inp.dt: HUGE;
  foreach_dimension()
    inp.pp->u.x = (idt*pa.u.x + itau*interpolate(u.x, xp, yp ,zp))/(idt + itau);
  // Gravity
  double rhof = is_constant(rho) ? constant(rho) : interpolate (rho, xp, yp, zp);
  foreach_dimension() 
    inp.pp->u.x += rhop ? (G.x*(rhop - rhof)/rhop)/(idt + itau) : 0;
  // No additional acceleration
  return (coord){0., 0., 0.};
}

#if TWO_WAY
event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

void add_force (Particles P, vector F) {
  particle_boundary (P);
  double rhof = constant(rho);
  foreach_particle_in(P) {
    double muc = is_constant (mu.x) ? constant(mu.x) : interpolate (mu.x, x, y, z);
    double rp = p().u2.y;
    double pref = 6*pi*muc*rp;
    double xp = x, yp = y , zp = z;
    Point point = locate (x,y,z);
    double rhop = p().u2.x;
    foreach_dimension()
      F.x[] += pref*(p().u.x - interpolate(u.x, xp, yp, zp));
    pref = 4./3.*pi*cube(rp)*(rhop-rhof);
   
  }
}

event acceleration (i++) {
  face vector av = a;
  vector Fp[];
  foreach() 
    foreach_dimension()
      Fp.x[] = 0;
  foreach_P_in_list (inertial_particles)
    add_force (P, Fp);  
  
  double rhof = constant(rho);
  boundary ((scalar*){Fp});
  foreach_face() 
    av.x[] += face_value(Fp.x, 0)/(dv()*rhof);
}
#endif