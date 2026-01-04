/**
# Verification of Frenet-Serret Frame

This test specifically verifies the computation of the Frenet-Serret frame and local coordinates
for a ring vortex filament against analytical solutions.
*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"

int main() {
  L0 = 4.0;
  N = 32;
  init_grid(N);

  // --- 1. Setup Ring Filament ---
  int turns = 1; 
  int nseg_per_turn = 256; 
  int nseg = nseg_per_turn * turns + 1;
  double R = 1.0;
  double dphi = 2 * pi / ((double)nseg_per_turn);
  double phi_arr[nseg];
  double a1[nseg];
  coord C1[nseg];

  for (int i = 0; i < nseg; i++) {
    phi_arr[i] = dphi * (double)i - pi; 
    C1[i].x = R * cos(phi_arr[i]);
    C1[i].y = R * sin(phi_arr[i]);
    C1[i].z = 0.0;
    a1[i] = 0.10;
  }

  coord xshift = {0, 0, 0}, dxshift = {0, 0, 0};
  struct vortex_filament filament1;
  allocate_vortex_filament_members(&filament1, nseg);
  initialize_filaments(filament1, nseg, dphi, phi_arr, a1, C1, xshift, dxshift);

  // --- 2. Verification Points ---
  // We test at specific angles to be deterministic
  int n_tests = 10;
  double test_angles[] = {0.0, pi/2.0, -pi/2.0, pi/4.0, pi/3.0, 2*pi/3.0, pi/6.0, 5*pi/6.0, 7*pi/6.0, 11*pi/6.0};
  
  double max_err_T = 0, max_err_N = 0, max_err_B = 0;
  double max_err_rho = 0, max_err_theta = 0;

  fprintf(stderr, "--- Verifying Frenet-Serret Frame ---\n");

  for (int i = 0; i < n_tests; i++) {
    double alpha = test_angles[i];
    
    // Define a point P slightly off the curve at this angle
    // P = C(alpha) + delta_r * e_r + delta_z * e_z
    double delta_r = 0.2;
    double delta_z = 0.1;
    
    // Cylindrical coords of P
    double r_p = R - delta_r; // Inside the ring (N points inward)
    double z_p = delta_z;
    
    coord P = {r_p * cos(alpha), r_p * sin(alpha), z_p};
    
    // Analytical Expectations
    // Closest point on curve is exactly at phi = alpha
    coord T_exact = {-sin(alpha), cos(alpha), 0};
    coord N_exact = {-cos(alpha), -sin(alpha), 0}; // Inward normal for circle
    coord B_exact = {0, 0, 1};
    
    // Local coordinates
    // Vector from curve to P: P - C(alpha)
    // C(alpha) = (R cos alpha, R sin alpha, 0)
    // P - C = ((R-dr)cos - R cos, (R-dr)sin - R sin, dz - 0)
    //       = (-dr cos, -dr sin, dz)
    // Project onto Frame:
    // x_N = (P-C).N = (-dr cos)(-cos) + (-dr sin)(-sin) = dr (cos^2 + sin^2) = dr
    // x_B = (P-C).B = dz
    double x_N_exact = delta_r; 
    double x_B_exact = delta_z;
    
    double rho_exact = sqrt(x_N_exact*x_N_exact + x_B_exact*x_B_exact);
    double theta_local_exact = atan2(x_B_exact, x_N_exact);

    // Compute Numerical
    struct vortex_filament params1 = filament1;
    params1.pcar = P;
    struct local_filament vortex1 = get_local_coordinates(spatial_period=0, max_distance=4*L0, vortex=&params1);

    if (vortex1.near) {
      double err_T = sqrt(vecdist2(vortex1.Tvec, T_exact));
      double err_N = sqrt(vecdist2(vortex1.Nvec, N_exact));
      double err_B = sqrt(vecdist2(vortex1.Bvec, B_exact));
      
      double err_rho = fabs(vortex1.rho - rho_exact);
      double err_theta = fabs(vortex1.theta - theta_local_exact);
      if (err_theta > pi) err_theta = 2*pi - err_theta;

      fprintf(stderr, "Angle %.2f: Err T=%.1e N=%.1e B=%.1e Rho=%.1e Theta=%.1e\n", 
              alpha, err_T, err_N, err_B, err_rho, err_theta);

      if (err_T > max_err_T) max_err_T = err_T;
      if (err_N > max_err_N) max_err_N = err_N;
      if (err_B > max_err_B) max_err_B = err_B;
      if (err_rho > max_err_rho) max_err_rho = err_rho;
      if (err_theta > max_err_theta) max_err_theta = err_theta;
    } else {
        fprintf(stderr, "Angle %.2f: Point not found near filament!\n", alpha);
        return 1;
    }
  }

  fprintf(stderr, "Max Errors:\n T: %g\n N: %g\n B: %g\n rho: %g\n theta: %g\n", 
          max_err_T, max_err_N, max_err_B, max_err_rho, max_err_theta);

  if (max_err_T > 1e-4 || max_err_N > 1e-4 || max_err_B > 1e-4 || max_err_rho > 1e-4) {
    fprintf(stderr, "Frenet-Serret Verification FAILED\n");
    return 1;
  } else {
    fprintf(stderr, "Frenet-Serret Verification PASSED\n");
  }

  free_vortex_filament_members(&filament1);
  return 0;
}
