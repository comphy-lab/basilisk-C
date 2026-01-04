/**
# Vorticity field for a Batchelor vortex 

In this example, we use the approach of [Blanco-Rodriguez et al. (2015)](#blanco2015internal) to compute the velocity and vorticity fields for a 
curved Batchelor vortex including first-order corrections from curvature.

<table>
<tr>
<td><center>![Axial vorticity at leading order](test_curvature_correction/w0_xi.png){ width="75%" }</center></td>
<td><center>![Axial vorticity at following order](test_curvature_correction/w1_xi.png){ width="75%" }</center></td>
</tr>
</table>

*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"
#include "acastillo/input_fields/draw_filaments.h"

int minlevel = 3;
int maxlevel = 6;

int main()
{
  L0 = 1.50;
  Z0 = -L0 / 2;
  N = 1 << minlevel;
  init_grid(N);

  scalar v_mag[], v_ang[];
  vector u[], omega[];
  vector u0[], omega0[];
  vector u1[], omega1[];
  foreach(){
    v_mag[] = 0;
    v_ang[] = 0;
    foreach_dimension(){
      u.x[] = 0;
      u0.x[] = 0;
      u1.x[] = 0;
      omega.x[] = 0;
      omega0.x[] = 0;
      omega1.x[] = 0;
    }
  }

  int turns = 4;
  int nseg_per_turn = 128;
  int nseg = nseg_per_turn*turns+1;
  double R=1.0;
  double H=0.0;
  double U_c = 1.0;
  double dphi = 2*pi/((double)nseg_per_turn);
  double phi[nseg]; 
  double a1[nseg];
  coord C1[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    phi[i] = dphi * (double)i - 2*pi;
    C1[i].x = R * cos(phi[i]);
    C1[i].y = R * sin(phi[i]);
    C1[i].z = 0.0;
    a1[i] = 0.10;
  } 

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C1, a1);
    save ("prescribed_curve.png");
  }

  // Initialize the vortex filament
  coord xshift = {0, 0, turns*H}, dxshift = {0, 0, 0};
  struct vortex_filament filament1;
  allocate_vortex_filament_members(&filament1, nseg);
  initialize_filaments(filament1, nseg, dphi, phi, a1, C1, xshift, dxshift);
  
  // Display the curve and the Frenet-Serret frame
  {
    view (camera="iso");  
    draw_space_curve_with_vectors(filament1.nseg, filament1.C, filament1.Tvec, filament1.Nvec, filament1.Bvec, scale=0.25);   
    save ("prescribed_curve_with_vectors.png");
  }

  FILE *fp = fopen("curve.txt", "w"); 
  for (int i = 0; i < nseg; i++){
    fprintf (fp, "%d %g %g %g %g %g %g %g %g %g \n", i, filament1.phi[i], 
         filament1.C[i].x, filament1.C[i].y, filament1.C[i].z, 
         filament1.sigma[i], filament1.kappa[i], filament1.tau[i], 
         filament1.s[i], filament1.theta0[i]);
  }
  
  // We refine close to the curve
  scalar dmin[];
  for (int i = (maxlevel-minlevel-1); i >= 0; i--){
    foreach(){      
      struct vortex_filament params1;
      params1 = filament1;
      params1.pcar = (coord){x,y,z};
      dmin[] = 0;
      dmin[] = (get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params1) < (i+1)*a1[0])*noise();    
    }
    
    adapt_wavelet ((scalar*){dmin}, (double[]){1e-12}, maxlevel-i, minlevel);
    
    {
      view (camera="iso"); 
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
    struct local_filament vortex1 = get_local_coordinates(spatial_period=1, max_distance=4*L0, vortex=&params1);
    
    if (vortex1.near == 1){

      integration_results corr ;
      compute_A_with_derivatives(vortex1.rho, vortex1.a, U_c, &corr);

      batchelor_vortex sol;  
      sol = calculate_vortex_flow(&vortex1, U_c, &corr);
      
      coord e_r, e_th;
      double phi = vortex1.theta - vortex1.theta0;
      foreach_dimension(){
        e_r.x  =  cos(phi) * vortex1.Nvec.x + sin(phi) * vortex1.Bvec.x;
        e_th.x = -sin(phi) * vortex1.Nvec.x + cos(phi) * vortex1.Bvec.x;
      }
      
      foreach_dimension(){
        u.x[] += sol.u0_r * e_r.x + sol.u0_th * e_th.x + sol.u0_xi * vortex1.Tvec.x;
        u.x[] += sol.u1_r * e_r.x + sol.u1_th * e_th.x + sol.u1_xi * vortex1.Tvec.x;

        omega.x[] += sol.w0_r * e_r.x + sol.w0_th * e_th.x + sol.w0_xi * vortex1.Tvec.x;
        omega.x[] += sol.w1_r * e_r.x + sol.w1_th * e_th.x + sol.w1_xi * vortex1.Tvec.x;
      }

      omega0.x[] = sol.w0_r;
      omega0.y[] = sol.w0_th;
      omega0.z[] = sol.w0_xi;

      omega1.x[] = sol.w1_r;
      omega1.y[] = sol.w1_th;
      omega1.z[] = sol.w1_xi;

      u0.x[] = sol.u0_r;
      u0.y[] = sol.u0_th;
      u0.z[] = sol.u0_xi;

      u1.x[] = sol.u1_r;
      u1.y[] = sol.u1_th;
      u1.z[] = sol.u1_xi;

      v_mag[] = vortex1.rho; 
      v_ang[] = phi; 
    }
  }
  
  
  
  {    
    view (camera="iso", fov=40);
    squares("u0.x", n = {1, 0, 0});
    squares("u0.x", n = {0, 1, 0});
    save("u0_r.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega0.x", n = {1, 0, 0});
    squares("omega0.x", n = {0, 1, 0});
    save("w0_r.png");
    clear();

    view (camera="iso", fov=40);
    squares("u1.x", n = {1, 0, 0});
    squares("u1.x", n = {0, 1, 0});
    save("u1_r.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega1.x", n = {1, 0, 0});
    squares("omega1.x", n = {0, 1, 0});
    save("w1_r.png");
    clear();
  }

  {    
    view (camera="iso", fov=40);
    squares("u0.y", n = {1, 0, 0});
    squares("u0.y", n = {0, 1, 0});
    save("u0_th.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega0.y", n = {1, 0, 0});
    squares("omega0.y", n = {0, 1, 0});
    save("w0_th.png");
    clear();

    view (camera="iso", fov=40);
    squares("u1.y", n = {1, 0, 0});
    squares("u1.y", n = {0, 1, 0});
    save("u1_th.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega1.y", n = {1, 0, 0});
    squares("omega1.y", n = {0, 1, 0});
    save("w1_th.png");
    clear();
  }

  {    
    view (camera="iso", fov=40);
    squares("u0.z", n = {1, 0, 0});
    squares("u0.z", n = {0, 1, 0});
    save("u0_xi.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega0.z", n = {1, 0, 0});
    squares("omega0.z", n = {0, 1, 0});
    save("w0_xi.png");
    clear();

    view (camera="iso", fov=40);
    squares("u1.z", n = {1, 0, 0});
    squares("u1.z", n = {0, 1, 0});
    save("u1_xi.png");
    clear();

    view (camera="iso", fov=40);
    squares("omega1.z", n = {1, 0, 0});
    squares("omega1.z", n = {0, 1, 0});
    save("w1_xi.png");
    clear();
  }
  
  free_vortex_filament_members(&filament1);
}

/**
# References

~~~bib

@article{blanco2015internal,
  title={Internal structure of vortex rings and helical vortices},
  author={Blanco-Rodr{\'\i}guez, Francisco J and Le Diz{\`e}s, St{\'e}phane and Sel{\c{c}}uk, Can and Delbende, Ivan and Rossi, Maurice},
  journal={Journal of Fluid Mechanics},
  volume={785},
  pages={219--247},
  year={2015},
  publisher={Cambridge University Press}
}

~~~
*/