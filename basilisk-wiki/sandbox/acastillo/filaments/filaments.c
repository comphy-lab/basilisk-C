#include "filaments.h"


double test_vorticity_filament0(coord pcar, int n_seg, double a, double* t0, coord* c, coord* tvec, coord* nvec, coord* bvec, int period){
  /* Compute the local coordinates required for the vorticity field.
  Each position P(x,y,z) is projected into the local Frenet-Serret frame to
  obtain a set of local coordinates, such that:
  i)  (P - X . T = 0
  ii) (P - X) . N = x_n
  ii) (P - X) . B = x_b
  This requires finding the value of X(t0) along each space curve that verifies
  i) through a minization process.
  Then, we use the local coordinates (x_n, x_b) to define a radial coordinate
  rho required to compute the vorticity of a Lamb-Oseen vortex as
  omega = Gamma/(pi a²) exp(-rho²/a²) . T
  where Gamma is the circulation and a the core size.
  */

  /* First, we approximate the minimal distance between the point P and each
  segment of the curve. */
  double dmin = 1e30, tmin=0;
  double rho_loc;
  for (int i = 0; i < n_seg; i++){
    if (vecdist2(pcar, c[i]) < dmin){
      dmin = vecdist2(pcar, c[i]);
      tmin = t0[i];
    }
  }

  if (period != 0)
    tmin = fmod(tmin + period*2*pi, period*2*pi);

  // If P is close to the vortex, we refine the initial guess
  coord ploc, ccar, frenet[3];
  double tq = frenet_projection_min(n_seg, a, t0, c, tvec, nvec, bvec, pcar, tmin);

  ccar = gsl_interp1d( n_seg, t0, c, tq);

  // Then, compute the local coordinates for the vortex
  frenet[0] = gsl_interp1d( n_seg, t0, tvec, tq);
  frenet[1] = gsl_interp1d( n_seg, t0, nvec, tq);
  frenet[2] = gsl_interp1d( n_seg, t0, bvec, tq);

  ploc.x = vecdot(vecdiff(pcar, ccar), frenet[0]);
  ploc.y = vecdot(vecdiff(pcar, ccar), frenet[1]);
  ploc.z = vecdot(vecdiff(pcar, ccar), frenet[2]);
  rho_loc = sqrt(vecdot(ploc, ploc));

  return rho_loc;
}


coord test_vorticity_filament1(coord pcar, int n_seg, double a, double* t0, coord* c, coord* tvec, coord* nvec, coord* bvec, int period){
  double dmin = 1e30, tmin=0;
  for (int i = 0; i < n_seg; i++){
    if (vecdist2(pcar, c[i]) < dmin){
      dmin = vecdist2(pcar, c[i]);
      tmin = t0[i];
    }
  }

  if (period != 0)
    tmin = fmod(tmin + period*2*pi, period*2*pi);

  // If P is close to the vortex, we refine the initial guess
  coord ccar, frenet[1];
  double tq = frenet_projection_min(n_seg, a, t0, c, tvec, nvec, bvec, pcar, tmin);

  ccar = gsl_interp1d( n_seg, t0, c, tq);

  // Then, compute the local coordinates for the vortex
  frenet[0] = gsl_interp1d( n_seg, t0, tvec, tq);

  return (coord) {frenet[0].x, frenet[0].y, frenet[0].z};
}































coord get_vorticity_filament(coord pcar, int n_seg, double a, double* t0, coord* c, coord* tvec, coord* nvec, coord* bvec, int period){
  /* Compute the local coordinates required for the vorticity field.
  Each position P(x,y,z) is projected into the local Frenet-Serret frame to
  obtain a set of local coordinates, such that:
  i)  (P - X . T = 0
  ii) (P - X) . N = x_n
  ii) (P - X) . B = x_b
  This requires finding the value of X(t0) along each space curve that verifies
  i) through a minization process.
  Then, we use the local coordinates (x_n, x_b) to define a radial coordinate
  rho required to compute the vorticity of a Lamb-Oseen vortex as
  omega = Gamma/(pi a²) exp(-rho²/a²) . T
  where Gamma is the circulation and a the core size.
  */

  /* First, we approximate the minimal distance between the point P and each
  segment of the curve. */
  coord omega;
  double dmin = 1e30, tmin=0;
  double rho_loc, omega_mag;
  for (int i = 0; i < n_seg; i++){
    if (vecdist2(pcar, c[i]) < dmin){
      dmin = vecdist2(pcar, c[i]);
      tmin = t0[i];
    }
  }

  if (period != 0)
    tmin = fmod(tmin + period*2*pi, period*2*pi);

  if (dmin < L0){
    // If P is close to the vortex, we refine the initial guess
    coord ploc, ccar, frenet[3];
    double tq = frenet_projection_min(n_seg, a, t0, c, tvec, nvec, bvec, pcar, tmin);

    ccar = gsl_interp1d( n_seg, t0, c, tq);

    // Then, compute the local coordinates for the vortex
    frenet[0] = gsl_interp1d( n_seg, t0, tvec, tq);
    frenet[1] = gsl_interp1d( n_seg, t0, nvec, tq);
    frenet[2] = gsl_interp1d( n_seg, t0, bvec, tq);

    ploc.x = vecdot(vecdiff(pcar, ccar), frenet[0]);
    ploc.y = vecdot(vecdiff(pcar, ccar), frenet[1]);
    ploc.z = vecdot(vecdiff(pcar, ccar), frenet[2]);
    rho_loc = sqrt(vecdot(ploc, ploc));

    // Last, we compute the vorticity for the vortex 1
    omega_mag = exp(-sq(rho_loc)/sq(a))/(pi*sq(a));
    omega = (coord) {omega_mag * frenet[0].x, omega_mag * frenet[0].y, omega_mag * frenet[0].z};
  }
  else {
    // Otherwise, if the point is too far, we set the vorticity to zero.
    foreach_dimension()
      omega.x = 0;
  }
  return omega;
}





coord get_vorticity_filament2(double Gamma, double Uc, coord pcar, int n_seg, double a, double* t0, coord* c, coord* tvec, coord* nvec, coord* bvec, int period){
  /* Compute the local coordinates required for the vorticity field.
  Each position P(x,y,z) is projected into the local Frenet-Serret frame to
  obtain a set of local coordinates, such that:
  i)  (P - X . T = 0
  ii) (P - X) . N = x_n
  ii) (P - X) . B = x_b
  This requires finding the value of X(t0) along each space curve that verifies
  i) through a minization process.
  Then, we use the local coordinates (x_n, x_b) to define a radial coordinate
  rho required to compute the vorticity of a Lamb-Oseen vortex as
  omega = Gamma/(pi a²) exp(-rho²/a²) . T + 2 Uc/a (rho/a) exp(-rho²/a²)
  where Gamma is the circulation and a the core size and Uc is the
  axial centerline velocity.
  */

  /* First, we approximate the minimal distance between the point P and each
  segment of the curve. */
  coord omega;
  double dmin = 1e30, tmin=0;
  double rho_loc, omega_mag;
  for (int i = 0; i < n_seg; i++){
    if (vecdist2(pcar, c[i]) < dmin){
      dmin = vecdist2(pcar, c[i]);
      tmin = t0[i];
    }
  }

  if (period != 0)
    tmin = fmod(tmin + period*2*pi, period*2*pi);

  if (dmin < L0){
    // If P is close to the vortex, we refine the initial guess
    coord ploc, ccar, frenet[3];
    double tq = frenet_projection_min(n_seg, a, t0, c, tvec, nvec, bvec, pcar, tmin);

    ccar = gsl_interp1d( n_seg, t0, c, tq);

    // Then, compute the local coordinates for the vortex
    frenet[0] = gsl_interp1d( n_seg, t0, tvec, tq);
    frenet[1] = gsl_interp1d( n_seg, t0, nvec, tq);
    frenet[2] = gsl_interp1d( n_seg, t0, bvec, tq);

    ploc.x = vecdot(vecdiff(pcar, ccar), frenet[0]);
    ploc.y = vecdot(vecdiff(pcar, ccar), frenet[1]);
    ploc.z = vecdot(vecdiff(pcar, ccar), frenet[2]);
    rho_loc = sqrt(vecdot(ploc, ploc));

    // Last, we compute the vorticity for the vortex 1
    omega_mag = exp(-sq(rho_loc)/sq(a))/(sq(a));
    omega.x = omega_mag * (Gamma/pi * frenet[0].x + (2.0 * Uc) * (-ploc.z * frenet[1].x + ploc.y * frenet[2].x));
    omega.y = omega_mag * (Gamma/pi * frenet[0].y + (2.0 * Uc) * (-ploc.z * frenet[1].y + ploc.y * frenet[2].y));
    omega.z = omega_mag * (Gamma/pi * frenet[0].z + (2.0 * Uc) * (-ploc.z * frenet[1].z + ploc.y * frenet[2].z));
  }
  else {
    // Otherwise, if the point is too far, we set the vorticity to zero.
    foreach_dimension()
      omega.x = 0;
  }
  return omega;
}
