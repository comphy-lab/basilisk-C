/**
# Vorticity field for a Batchelor vortex with jet component

In this example, we define a space-curve $\mathcal{C}(\xi,t)$ and compute a
local Frenet-Serret basis ($\bf\hat{T}, \hat{N}, \hat{B}$). For each point in
the computation domain, $\vec{x}$, we project into a curvilinear orthonormal
basis $({\bf\hat{e}}_\rho, {\bf\hat{e}}_\varphi, {\bf\hat{e}}_T)$, and write
the vorticity field for Batchelor vortex of unit Circulation:

$$
\begin{aligned}
\vec{\omega}(\vec{x},t) = \frac{\Gamma}{\pi a^2}e^{-\rho^2/a^2} {\bf\hat{T}} +
\frac{2 U_c}{a}\frac{\rho}{a}e^{-\rho^2/a^2} (-\sin{\phi}{\bf\hat{N}}+\cos{\phi}{\bf\hat{B}})
\end{aligned}
$$
where $a$ is the size of the vortex core.

<table>
<tr>
<td><center>![Iso-surface of the tangential vorticity](test_draw_filaments4/vorticity_tan.png){ width="75%" }</center></td>
<td><center>![Iso-surface of the angular vorticity](test_draw_filaments4/vorticity_ang.png){ width="75%" }</center></td>
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

  scalar omega1_mag[], omega2_mag[];
  vector omega1[], omega2[];
  foreach(){
    omega1_mag[] = 0;
    omega2_mag[] = 0;
    foreach_dimension(){
      omega1.x[] = 0;
      omega2.x[] = 0;
    }
  }

  int turns = 4;
  int nseg_per_turn = 128;
  int nseg = nseg_per_turn*turns+1;
  double R=1.0;
  double H=pi;
  double dtheta = 2*pi/((double)nseg_per_turn);
  double theta[nseg]; 
  double a1[nseg];
  coord C1[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i - 2*pi;
    C1[i].x = R * cos(theta[i]);
    C1[i].y = R * sin(theta[i]);
    C1[i].z = (H/(2*pi)) * theta[i] - L0/2;
    a1[i] = 0.020;
  } 

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C1, a1);
    save ("prescribed_curve.png");
  }

  // Initialize the vortex filament
  coord xshift = {0, 0, turns*H}, dxshift = {0, 0, 0};
  struct vortex_filament filament1;
  allocate_vortex_filament(&filament1, nseg);
  initialize_filaments(filament1, nseg, dtheta, theta, a1, C1, xshift, dxshift);

  // Display the curve and the Frenet-Serret frame
  {
    view (camera="iso");  
    draw_space_curve_with_vectors(filament1.nseg, filament1.C, filament1.Tvec, filament1.Nvec, filament1.Bvec, scale=0.25);   
    save ("prescribed_curve_with_vectors.png");
  }

  FILE *fp = fopen("curve.txt", "w"); 
  for (int i = 0; i < nseg; i++){
    fprintf (fp, "%d %g %g %g %g %g %g %g %g %g \n", i, filament1.theta[i], 
         filament1.C[i].x, filament1.C[i].y, filament1.C[i].z, 
         filament1.sigma[i], filament1.kappa[i], filament1.tau[i], 
         filament1.s[i], filament1.varphi0[i]);
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
      cells(n = {1,0,0});     
      cells(n = {0,1,0});     
      cells(n = {0,0,1});     
      save ("cells.png"); 
      clear();
    }
  }

  double Uc = 1.0;
  // 4. Get the local coordinates in the Frenet-Serret Frame  
  foreach(){
    struct vortex_filament params1;
    params1 = filament1;
    params1.pcar = (coord){x,y,z};    
    struct local_filament vortex1 = get_local_coordinates(spatial_period=0, max_distance=4*L0, vortex=&params1);
    
    if (vortex1.near == 1){
      
      // We may use the coordinates to compute the vorticity field
      omega1_mag[] = exp(-sq(vortex1.rho/vortex1.a))/(pi*sq(vortex1.a));
      omega2_mag[] = exp(-sq(vortex1.rho/vortex1.a))*(2.0*Uc*vortex1.rho/sq(vortex1.a));
      
      double phi = vortex1.phi + vortex1.varphi0;
      foreach_dimension(){
        omega1.x[] = omega1_mag[] * vortex1.Tvec.x;
        omega2.x[] = omega2_mag[] * (-vortex1.Nvec.x*sin(phi) + vortex1.Bvec.x*cos(phi));    
      }
    }
  }
  restriction ((scalar*){omega1, omega2});

  {
    vectors (u = "omega1", scale = 1);
    isosurface ("omega1_mag",   0.01, color="omega1_mag");
    isosurface ("omega1_mag",   0.10, color="omega1_mag");
    isosurface ("omega1_mag",   1.00, color="omega1_mag");
    save("vorticity_tan.png");
    clear();
  }

  {
    vectors (u = "omega2", scale = 1);
    isosurface ("omega2_mag",   0.01, color="omega2_mag");
    isosurface ("omega2_mag",   0.10, color="omega2_mag");
    isosurface ("omega2_mag",   1.00, color="omega2_mag");
    save("vorticity_ang.png");
    clear();
  }
  
  free_vortex_filament(&filament1);
}


/** 

## Outputs

~~~gnuplot Prescribed curbe $\mathcal{C(\xi,t)}$
set term pngcairo enhanced size 640,480 font ",12"
set output 'curve.png'
plot 'curve.txt' u 1:3 w lp title "C.x",\
              "" u 1:4 w lp title "C.y",\
              "" u 1:5 w lp title "C.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Strecthing coefficient $\sigma(\xi,t)$
set term pngcairo enhanced size 640,480 font ",12"
set output 'sigma.png'
plot 'curve.txt' u 1:6 w lp title "sigma"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Curvature $\kappa(\xi,t)$
set term pngcairo enhanced size 640,480 font ",12"
set output 'kappa.png'
plot 'curve.txt' u 1:7 w lp title "kappa"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Torsion $\tau(\xi,t)$
set term pngcairo enhanced size 640,480 font ",12"
set output 'tau.png'
plot 'curve.txt' u 1:8 w lp title "tau"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Arc-lenght parameter $\ell(\xi,t)$
set term pngcairo enhanced size 640,480 font ",12"
set output 'arclength.png'
plot 'curve.txt' u 1:9 w lp title "s"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Cumulative torsion $\varphi_0(\xi,t)$
set term pngcairo enhanced size 640,480 font ",12"
set output 'varphi.png'
plot 'curve.txt' u 1:10 w lp title "varphi0"
set xlabel "index"
set ylabel "value"
~~~

*/