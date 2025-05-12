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
#include "filaments.h"
#include "draw_filaments.h"

int minlevel = 4;
int maxlevel = 8;

int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << minlevel;
  init_grid(N);

  scalar omega1_mag[], omega2_mag[], omega3_mag[];
  vector omega[];
  foreach(){
    omega1_mag[] = 0;
    omega2_mag[] = 0;
    omega3_mag[] = 0;
    foreach_dimension(){
      omega.x[] = 0;
    }
  }

  int turns = 4;
  int nseg_per_turn = 128;
  int nseg = nseg_per_turn*turns+1;
  double R=1.0;
  double H=pi;
  double dtheta = 2*pi/((double)nseg_per_turn);
  double theta[nseg];
  double a1[nseg], a2[nseg], a3[nseg];
  coord C1[nseg], C2[nseg], C3[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i - 2*pi;
    C1[i].x = R * cos(theta[i]);
    C1[i].y = R * sin(theta[i]);
    C1[i].z = (H/(2*pi)) * theta[i] - L0/2;
    a1[i] = 0.020;

    C2[i].x = R * cos(theta[i] + 2*pi/3);
    C2[i].y = R * sin(theta[i] + 2*pi/3);
    C2[i].z = (H/(2*pi)) * theta[i] - L0/2;
    a2[i] = 0.030;

    C3[i].x = R * cos(theta[i] + 4*pi/3);
    C3[i].y = R * sin(theta[i] + 4*pi/3);
    C3[i].z = (H/(2*pi)) * theta[i] - L0/2;
    a3[i] = 0.020;
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
  allocate_vortex_filament(&filament1, nseg);
  initialize_filaments(filament1, nseg, dtheta, theta, a1, C1, xshift, dxshift);

  struct vortex_filament filament2;
  allocate_vortex_filament(&filament2, nseg);
  initialize_filaments(filament2, nseg, dtheta, theta, a2, C2, xshift, dxshift);

  struct vortex_filament filament3;
  allocate_vortex_filament(&filament3, nseg);
  initialize_filaments(filament3, nseg, dtheta, theta, a3, C3, xshift, dxshift);
  
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

    struct vortex_filament params2;
    params2 = filament2;
    params2.pcar = (coord){x,y,z};
    
    struct vortex_filament params3;    
    params3 = filament3;
    params3.pcar = (coord){x,y,z};  
    
    struct local_filament vortex1 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params1);
    struct local_filament vortex2 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params2);
    struct local_filament vortex3 = get_local_coordinates(spatial_period=0, max_distance=0.5, vortex=&params3);
    
    // 5. We use the coordinates to compute the vorticity field
    if (vortex1.near == 1)
      omega1_mag[] = exp(-sq(vortex1.rho/vortex1.a))/(pi*sq(vortex1.a));
    
    if (vortex2.near == 1)
      omega2_mag[] = exp(-sq(vortex2.rho/vortex2.a))/(pi*sq(vortex2.a));
    
    if (vortex3.near == 1)
      omega3_mag[] = exp(-sq(vortex3.rho/vortex3.a))/(pi*sq(vortex3.a));
    
    foreach_dimension()
      omega.x[] = omega1_mag[] * vortex1.Tvec.x 
      + omega2_mag[] * vortex2.Tvec.x    
      + omega3_mag[] * vortex3.Tvec.x;    
  }  
  restriction ((scalar*){omega});
  
  {    
    isosurface ("omega1_mag",   1.00, color="omega1_mag");
    isosurface ("omega2_mag",   1.00, color="omega2_mag");
    isosurface ("omega3_mag",   1.00, color="omega3_mag");
    save("vorticity.png");
    clear();
  }
  
  free_vortex_filament(&filament1);
  free_vortex_filament(&filament2);
  free_vortex_filament(&filament3);
}