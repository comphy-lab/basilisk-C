/**
# Motion of a elliptical vortex ring (using the filament framework)

In this example, we consider the motion of a elliptical vortex ring of radii
$R_1$ and $R_2$, and nominal core size $a$ using the [vortex filament
framework](../input_fields/filaments.h).

The filament approach is perfectly justified as long as the core size remains
small compared to other spatial scales (local curvature). We use the same
framework as in [Durán Venegas & Le Dizès (2019)](#duran2019),
[Castillo-Castellanos, Le Dizès & Durán Venegas (2021)](#castillo2021), and
[Castillo-Castellanos & Le Dizès (2019)](#castillo2022). All the vorticity is
considered as being concentrated along space-curves $\vec{x}\in\mathcal{C}$ 
which move as material lines in the fluid according to:
$$
\begin{aligned}
\frac{d\vec{x}_c}{dt} = \vec{U}_{ind} + \vec{U}_{\infty}
\end{aligned}
$$
where $\vec{x}_c$ is the position vector of vortex filament, $\vec{U}_{ind}$ the
velocity induced by the vortex filament and $\vec{U}_{\infty}$ an external
velocity field. The induced velocity $\vec{U}_{ind}$ is given by the 
[Biot-Savart](biot-savart.h) law. In this example, we set $\vec{U}_{\infty}$
as the advection velocity of a vortex ring with the same core size. 

*/

#include "grid/octree.h"
#include "run.h"
#include "view.h"
#include "../input_fields/filaments.h"
#include "../input_fields/draw_filaments.h"
#include "biot-savart.h"

/**
 In this example, $R_1=1.1$, $R_2=0.9$, $a=0.05$ and the vortex ring is
 discretized into $n_{seg}=64$ filaments.
*/
int nseg = 64;
double R1=1.10;
double R2=0.90;
double a=0.05;
double dtmax=0.01;
struct vortex_filament filament1;
struct vortex_filament filament1_prev;
coord Uinfty = {0., 0., 0.365405};

/**
 The main time loop is defined in [run.h]().
*/
int minlevel = 4;
int maxlevel = 8;
vector u[];
int main(){
  L0 = 10.0;
  X0 = Y0 = Z0 = -L0 / 2;
  
  N = 1 << minlevel;  
  init_grid(N);  
  run();
}

/**
## Initial conditions 
We consider the space-curve $\mathcal{C}(\xi,t)$ is parametrized as function of
$\theta(\xi,t)$. At time $t=0$,
$$
\begin{aligned}
x &= R_1\cos(\theta), \quad
\\
y &= R_2\sin(\theta), \quad
\\
z &= 0
\end{aligned}
$$
We also set the core size $a_{seg}$ such that the total vorticity is that of a
vortex ring of nominal core size $a_{ring}$, such that each segment has constant
volume
$$
\begin{aligned}
\pi a_{seg}^2 \ell_{seg} = a_{ring}^2 R \frac{2\pi^2}{n_{seg}} 
\end{aligned}
$$
This also means that vortex stretching with modify the local core size.

The curve $\mathcal{C}(\xi,t)$ will be stored as `struct vortex_filament` which 
must be released at the end of the simulation.
*/
event init (t = 0) {
  double dtheta = 2*pi/((double)nseg);
  double theta[nseg];
  double a1[nseg];
  double vol1[nseg];
  coord C1[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i;
    C1[i] = (coord) { R1 * cos(theta[i]), R2 * sin(theta[i]), 0.};
    vol1[i] = pi * sq(a) * ((R1 + R2)/2.) * dtheta;
    a1[i] = a;
  } 
  
  // We store the space-curve in a structure   
  allocate_vortex_filament_members(&filament1, nseg);
  allocate_vortex_filament_members(&filament1_prev, nseg);
  memcpy(filament1.theta, theta, nseg * sizeof(double));
  memcpy(filament1.C, C1, nseg * sizeof(coord));
  memcpy(filament1.a, a1, nseg * sizeof(double));
  memcpy(filament1.vol, vol1, nseg * sizeof(double));  

  local_induced_velocity(filament1);  
  for (int j = 0; j < nseg; j++) {
    filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
  }
  
    
  view (camera="iso");
  draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);
  save ("prescribed_curve.png");

  FILE *fp = fopen("curve.txt", "w"); 
  fclose(fp);  

  /**
  We also initialize the Cartesian grid close to the vortex filament and 
  initialize the velocity field to zero.
  */ 
  scalar dmin[];
  for (int i = (maxlevel-minlevel-1); i >= 0; i--){
    foreach(){      
      struct vortex_filament params1;
      params1 = filament1;
      params1.pcar = (coord){x,y,z};
      dmin[] = 0;
      dmin[] = (get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params1) < (i+1)*filament1.a[0])*noise();    
    }    
    adapt_wavelet ((scalar*){dmin}, (double[]){1e-12}, maxlevel-i, minlevel);    
  }  
  foreach(){
    foreach_dimension(){
      u.x[] = 0.;
    }
  }
}

event finalize(t = end){
  free_vortex_filament_members(&filament1);
  free_vortex_filament_members(&filament1_prev);
}

/** 
## A Simple Time Advancing Scheme

Time integration is done in two steps. First, we evaluate the (self-)induced
velocity at $\vec{x}_c$
*/
event evaluate_velocity (i++) {  

  memcpy(filament1_prev.Utotal, filament1.Utotal, nseg * sizeof(coord));

  local_induced_velocity(filament1);  
  for (int j = 0; j < nseg; j++) {
    filament1.Uauto[j] = nonlocal_induced_velocity(filament1.C[j], filament1);            
    foreach_dimension(){
      filament1.Utotal[j].x = filament1.Uauto[j].x + filament1.Ulocal[j].x - Uinfty.x;      
    }
  }
}

/**
Then, we advect the filament segments using an explicit Adams-Bashforth scheme
and update the core size
*/ 
event advance_filaments (i++, last) {      
  for (int j = 0; j < nseg; j++) {
    foreach_dimension(){
      filament1.C[j].x += dt*(3.*filament1.Utotal[j].x - filament1_prev.Utotal[j].x)/2.;
    }
  } 
  local_induced_velocity(filament1);  
  for (int j = 0; j < nseg; j++) {
    filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
  } 
  dt = dtnext (dtmax);
}

/** 
## Outputs

### Displacement of the vortex filament

The elliptical vortex ring is characterized by this flapping motion.

![Motion of the elliptical vortex ring](test_vortex_ellipse1/vortex_ellipse.mp4)

*/

event sequence (t += 0.10){
  view(camera="iso", fov=10);
  draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);  
  save("vortex_ellipse.mp4");
}

/**

### The solution outside the vortex tube

Vortex filaments are dealt as Lagrangian particles. We may also use an Eulerian
grid to evaluate the induced velocity at some point $\vec{x}\in\mathbb{R}^3$.

<center>
  <table>
  <tr>
  <td><center>![](test_vortex_ellipse1/axial_velocity_1.png){ width="100%" }</center></td>
  <td><center>![](test_vortex_ellipse1/axial_velocity_2.png){ width="100%" }</center></td>
  <td><center>![](test_vortex_ellipse1/axial_velocity_3.png){ width="100%" }</center></td>
  <td><center>![](test_vortex_ellipse1/axial_velocity_4.png){ width="100%" }</center></td>
  <td><center>![](test_vortex_ellipse1/axial_velocity_5.png){ width="100%" }</center></td>
  </tr>
  </table>
  <center>Axial velocity induced at the mid-plane taken at regular intervals </center>
</center>

*/

event slice (t += 1.0, t <= 5.0){  
  foreach(){
    coord p = {x,y,z};
    coord u_BS = nonlocal_induced_velocity(p, filament1);
    foreach_dimension(){
      u.x[] = u_BS.x - Uinfty.x;
    }
  }
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, maxlevel, minlevel);

  {
    char filename2[100]; 
    sprintf(filename2, "axial_velocity_%g.png", t);
    
    view(camera="right");  
    squares ("u.z", linear = true, spread = -1, n={1,0,0}, min=-0.33, max=0.33, map = cool_warm); 
    draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);  
    save(filename2);  
  }
}

/** 
### Time evolution

We may also follow the position of $\mathcal{C}$ and the induced velocities 
over time, 
*/

event final (t += 0.2, t <= 30.0){
  if (pid() == 0){
    FILE *fp = fopen("curve.txt", "a"); 
    write_filament_state(fp, &filament1);
    fclose(fp);
  }
}

/** 

~~~gnuplot Time evolution of the core size $a$ for a point initially located at $\theta=0$
set term pngcairo enhanced size 640,480 font ",12"
set output 'plot_core.png'
set xlabel "time"
set ylabel "core size"
plot 'curve.txt' u 1:18 every ::1::1 w lp title "\\theta=0"
~~~

~~~gnuplot Time evolution of the axial coordinate for a point initially located at $\theta=0$
set term pngcairo enhanced size 640,480 font ",12"
set output 'plot_position.png'
set xlabel "time"
set ylabel "axial coordinate"
plot 'curve.txt' u 1:5 every ::1::1 w lp title "\\theta=0"
~~~
*/

/**
# References

~~~bib

@article{duran2019, 
  title={Generalized helical vortex pairs},
  volume={865}, 
  DOI={10.1017/jfm.2019.65}, 
  journal={Journal of Fluid Mechanics},
  author={Durán Venegas, E. and Le Dizès, S.}, 
  year={2019}, 
  pages={523–545}
}

@article{castillo2021,
  title = {Closely spaced corotating helical vortices: General solutions},
  author = {Castillo-Castellanos, A. and Le Diz\`es, S. and Dur\'an Venegas, E.},
  journal = {Phys. Rev. Fluids},
  volume = {6},
  issue = {11},
  pages = {114701},
  numpages = {22},
  year = {2021},
  month = {Nov},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevFluids.6.114701},
  url = {https://link.aps.org/doi/10.1103/PhysRevFluids.6.114701}
}

@article{castillo2022, 
  title={Closely spaced co-rotating helical vortices: long-wave instability}, 
  volume={946},
  DOI={10.1017/jfm.2022.569}, 
  journal={Journal of Fluid Mechanics},
  author={Castillo-Castellanos, A. and Le Dizès, S.}, 
  year={2022}, 
  pages={A10}
}



~~~
*/