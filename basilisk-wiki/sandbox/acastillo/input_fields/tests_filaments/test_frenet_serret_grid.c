/**
# Verification of Frenet-Serret Frame on Grid

This test verifies the computation of the Frenet-Serret frame and local coordinates
for a ring vortex filament against analytical solutions for every point on the grid.
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

  // --- 2. Verification Loop ---
  fprintf(stderr, "--- Verifying Frenet-Serret Frame on Grid ---\n");

  double max_err_T = 0, max_err_N = 0, max_err_B = 0;
  double max_err_rho = 0, max_err_theta = 0;
  long count = 0;

  foreach() {
    // Avoid singularity at the center of the ring (z-axis) where theta is undefined/ambiguous for the frame
    // In our case, the closest point on the ring is well-defined everywhere except the z-axis (rho=0 in cylindrical)
    // However, for a ring in z=0, the closest point is unique everywhere except on the z-axis.
    double r_cyl = sqrt(x*x + y*y);
    if (r_cyl < 1e-6) continue; 

    // Analytical Expectations
    // For a point P(x,y,z), the closest point on the ring (R, 0, 0) is at angle atan2(y,x)
    double theta_cyl = atan2(y, x);
    
    // Closest point on curve
    coord P_curve = {R * cos(theta_cyl), R * sin(theta_cyl), 0};
    
    // Frenet Frame at P_curve
    // T is tangent to the circle: (-sin, cos, 0)
    coord T_exact = {-sin(theta_cyl), cos(theta_cyl), 0};
    // N is normal (inward for Frenet): (-cos, -sin, 0)
    coord N_exact = {-cos(theta_cyl), -sin(theta_cyl), 0}; 
    // B is binormal: (0, 0, 1)
    coord B_exact = {0, 0, 1};
    
    // Local coordinates
    // Vector from curve to P
    coord vec = {x - P_curve.x, y - P_curve.y, z - P_curve.z};
    
    // Project onto Frame
    double x_N = vecdot(vec, N_exact);
    double x_B = vecdot(vec, B_exact);
    
    double rho_exact = sqrt(x_N*x_N + x_B*x_B);
    double theta_local_exact = atan2(x_B, x_N);

    // Compute Numerical
    struct vortex_filament params1 = filament1;
    params1.pcar = (coord){x, y, z};
    struct local_filament vortex1 = get_local_coordinates(spatial_period=0, max_distance=4*L0, vortex=&params1);

    if (vortex1.near) {
      double err_T = sqrt(vecdist2(vortex1.Tvec, T_exact));
      double err_N = sqrt(vecdist2(vortex1.Nvec, N_exact));
      double err_B = sqrt(vecdist2(vortex1.Bvec, B_exact));
      
      double err_rho = fabs(vortex1.rho - rho_exact);
      double err_theta = fabs(vortex1.theta - theta_local_exact);
      if (err_theta > pi) err_theta = 2*pi - err_theta;

      if (err_T > max_err_T) max_err_T = err_T;
      if (err_N > max_err_N) max_err_N = err_N;
      if (err_B > max_err_B) max_err_B = err_B;
      if (err_rho > max_err_rho) max_err_rho = err_rho;
      if (err_theta > max_err_theta) max_err_theta = err_theta;
      
      count++;
    }
  }

  fprintf(stderr, "Checked %ld points.\n", count);
  fprintf(stderr, "Max Errors:\n T: %g\n N: %g\n B: %g\n rho: %g\n theta: %g\n", 
          max_err_T, max_err_N, max_err_B, max_err_rho, max_err_theta);

  // Thresholds might need to be slightly looser for grid points far away or due to discretization?
  // But analytical solution is exact for a ring.
  if (max_err_T > 1e-4 || max_err_N > 1e-4 || max_err_B > 1e-4 || max_err_rho > 1e-4) {
    fprintf(stderr, "Grid Verification FAILED\n");
    return 1;
  } else {
    fprintf(stderr, "Grid Verification PASSED\n");
  }

  free_vortex_filament_members(&filament1);
  return 0;
}
