/**
# Vorticity field for a Batchelor vortex 

In this example, we define a space-curve $\mathcal{C}(\xi,t)$ and compute a
local Frenet-Serret basis ($\bf\hat{T}, \hat{N}, \hat{B}$). For each point in
the computation domain, $\vec{x}$, we project into a curvilinear orthonormal
basis $({\bf\hat{e}}_\rho, {\bf\hat{e}}_\varphi, {\bf\hat{e}}_T)$, and write
the vorticity field for Batchelor vortex of unit Circulation:

$$
\begin{aligned}
\vec{\omega}(\vec{x},t) = \frac{1}{\pi a^2}e^{-\rho^2/a^2} {\bf\hat{T}}
\end{aligned}
$$
where $a$ is the size of the vortex core.

<table>
<tr>
<td><center>![Iso-surface of the vorticity magnitude](test_draw_filaments5/vorticity.png){ width="75%" }</center></td>
</tr>
</table>

*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"
#include "acastillo/input_fields/draw_filaments.h"

int minlevel = 4;
int maxlevel = 8;

int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << minlevel;
  init_grid(N);

  scalar v_mag[], v_ang[], omega1_mag[];
  vector omega[];
  foreach(){
    v_mag[] = 0;
    v_ang[] = 0;
    omega1_mag[] = 0;
    foreach_dimension(){
      omega.x[] = 0;
    }
  }

  int turns = 4;
  int nseg_per_turn = 128;
  int nseg = nseg_per_turn*turns+1;
  double R=1.0;
  double H=pi;
  double U_c = 0.0;
  double dphi = 2*pi/((double)nseg_per_turn);
  double phi[nseg];
  double a1[nseg], a2[nseg], a3[nseg];
  coord C1[nseg], C2[nseg], C3[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    phi[i] = dphi * (double)i - 2*pi;
    C1[i].x = R * cos(phi[i]);
    C1[i].y = R * sin(phi[i]);
    C1[i].z = (H/(2*pi)) * phi[i] - L0/2;
    a1[i] = 0.10;

    C2[i].x = R * cos(phi[i] + 2*pi/3);
    C2[i].y = R * sin(phi[i] + 2*pi/3);
    C2[i].z = (H/(2*pi)) * phi[i] - L0/2;
    a2[i] = 0.10;

    C3[i].x = R * cos(phi[i] + 4*pi/3);
    C3[i].y = R * sin(phi[i] + 4*pi/3);
    C3[i].z = (H/(2*pi)) * phi[i] - L0/2;
    a3[i] = 0.10;
  } 

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C1, a1);
    draw_tube_along_curve(nseg, C2, a2);
    draw_tube_along_curve(nseg, C3, a3);
    save ("prescribed_curve.png");
  }

  // Initialize the vortex filament
  coord xshift = {0, 0, turns*H}, dxshift = {0, 0, 0};
  struct vortex_filament filament1;
  allocate_vortex_filament_members(&filament1, nseg);
  initialize_filaments(filament1, nseg, dphi, phi, a1, C1, xshift, dxshift);

  struct vortex_filament filament2;
  allocate_vortex_filament_members(&filament2, nseg);
  initialize_filaments(filament2, nseg, dphi, phi, a2, C2, xshift, dxshift);

  struct vortex_filament filament3;
  allocate_vortex_filament_members(&filament3, nseg);
  initialize_filaments(filament3, nseg, dphi, phi, a3, C3, xshift, dxshift);
  
  // Display the curve and the Frenet-Serret frame
  {
    view (camera="iso");  
    draw_space_curve_with_vectors(filament1.nseg, filament1.C, filament1.Tvec, filament1.Nvec, filament1.Bvec, scale=0.25);   
    draw_space_curve_with_vectors(filament2.nseg, filament2.C, filament2.Tvec, filament2.Nvec, filament2.Bvec, scale=0.25);   
    draw_space_curve_with_vectors(filament3.nseg, filament3.C, filament3.Tvec, filament3.Nvec, filament3.Bvec, scale=0.25);   
    save ("prescribed_curve_with_vectors.png");
  }

  // We refine close to the curves
  scalar dmin[];
  for (int i = (maxlevel-minlevel-1); i >= 0; i--){
    foreach(){      
      struct vortex_filament params1;
      params1 = filament1;
      params1.pcar = (coord){x,y,z};

      struct vortex_filament params2;
      params2 = filament2;
      params2.pcar = (coord){x,y,z};

      struct vortex_filament params3;
      params3 = filament3;
      params3.pcar = (coord){x,y,z};

      dmin[] = 0;
      dmin[] =  (get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params1) < (i+1)*a1[0])*noise();
      dmin[] += (get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params2) < (i+1)*a2[0])*noise();
      dmin[] += (get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params3) < (i+1)*a3[0])*noise();
    }

    adapt_wavelet ((scalar*){dmin}, (double[]){1e-12}, maxlevel-i, minlevel);
    
    {
      cells(n = {1,0,0});     
      cells(n = {0,1,0});     
      cells(n = {0,0,1});     
      save ("cells.png"); 
      clear();
    }
  }

  // 4. Get the local coordinates in the Frenet-Serret Frame  
  foreach(){    
    struct vortex_filament params1;
    params1 = filament1;
    params1.pcar = (coord){x,y,z};
    struct local_filament vortex1 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params1);

    struct vortex_filament params2;
    params2 = filament2;
    params2.pcar = (coord){x,y,z};
    struct local_filament vortex2 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params2);
    
    struct vortex_filament params3;    
    params3 = filament3;
    params3.pcar = (coord){x,y,z};  
    struct local_filament vortex3 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params3);

    if (vortex1.near == 1){
      coord e_r, e_th;
      double phi = vortex1.theta - vortex1.theta0;
      foreach_dimension(){
        e_r.x  =  cos(phi) * vortex1.Nvec.x + sin(phi) * vortex1.Bvec.x;
        e_th.x = -sin(phi) * vortex1.Nvec.x + cos(phi) * vortex1.Bvec.x;
      }

      // 5. We use the coordinates to compute the vorticity field
      integration_results corr ;
      compute_A_with_derivatives(vortex1.rho, vortex1.a, U_c, &corr);

      batchelor_vortex sol;  
      sol = calculate_vortex_flow(&vortex1, U_c, &corr);     

      foreach_dimension(){
        omega.x[] += sol.w0_r * e_r.x + sol.w0_th * e_th.x + sol.w0_xi * vortex1.Tvec.x;
        omega.x[] += sol.w1_r * e_r.x + sol.w1_th * e_th.x + sol.w1_xi * vortex1.Tvec.x;
      }
    }
      
    
    if (vortex2.near == 1){
      coord e_r, e_th;
      double phi = vortex2.theta - vortex2.theta0;
      foreach_dimension(){
        e_r.x  =  cos(phi) * vortex2.Nvec.x + sin(phi) * vortex2.Bvec.x;
        e_th.x = -sin(phi) * vortex2.Nvec.x + cos(phi) * vortex2.Bvec.x;
      }

      // 5. We use the coordinates to compute the vorticity field
      integration_results corr ;
      compute_A_with_derivatives(vortex2.rho, vortex2.a, U_c, &corr);

      batchelor_vortex sol;  
      sol = calculate_vortex_flow(&vortex2, U_c, &corr);     

      foreach_dimension(){
        omega.x[] += sol.w0_r * e_r.x + sol.w0_th * e_th.x + sol.w0_xi * vortex2.Tvec.x;
        omega.x[] += sol.w1_r * e_r.x + sol.w1_th * e_th.x + sol.w1_xi * vortex2.Tvec.x;
      } 
    }
      
    
    if (vortex3.near == 1){
      coord e_r, e_th;
      double phi = vortex3.theta - vortex3.theta0;
      foreach_dimension(){
        e_r.x  =  cos(phi) * vortex3.Nvec.x + sin(phi) * vortex3.Bvec.x;
        e_th.x = -sin(phi) * vortex3.Nvec.x + cos(phi) * vortex3.Bvec.x;
      }

      // 5. We use the coordinates to compute the vorticity field
      integration_results corr ;
      compute_A_with_derivatives(vortex3.rho, vortex3.a, U_c, &corr);

      batchelor_vortex sol;  
      sol = calculate_vortex_flow(&vortex3, U_c, &corr);     

      foreach_dimension(){
        omega.x[] += sol.w0_r * e_r.x + sol.w0_th * e_th.x + sol.w0_xi * vortex3.Tvec.x;
        omega.x[] += sol.w1_r * e_r.x + sol.w1_th * e_th.x + sol.w1_xi * vortex3.Tvec.x;
      } 
    }
    
    foreach_dimension()
      omega1_mag[] = sq(omega.x[]);    
    omega1_mag[] = sqrt(omega1_mag[]);
  } 
  restriction ((scalar*){omega});
  
  {    
    isosurface ("omega1_mag",   1.00, color="omega1_mag");
    save("vorticity.png");
    clear();
  }
  
  free_vortex_filament_members(&filament1);
  free_vortex_filament_members(&filament2);
  free_vortex_filament_members(&filament3);
}
