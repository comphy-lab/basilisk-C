/**
# Motion of two vortex rings (using the filament framework)

In this example, we consider the motion of two vortex rings of radius $R$ and
core size $a$ using the [vortex filament
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
as the advection velocity of one vortex ring with the same core size. 

*/

#include "grid/octree.h"
#include "run.h"
#include "view.h"
#include "../input_fields/filaments.h"
#include "../input_fields/draw_filaments.h"
#include "biot-savart.h"

/**
 In this example, $R_i=1.0$, $a_i=0.05$ and the vortex rings are discretized
 into $n_{seg}=64$ filaments.
*/
int nseg = 64;
double R = 1.0;
double d = 2.0;
double a=0.05;
double dtmax=0.005;
struct vortex_filament filament1;
struct vortex_filament filament1_prev;

struct vortex_filament filament2;
struct vortex_filament filament2_prev;
coord Uinfty = {0., 0., 0.365405};

/**
 The main time loop is defined in [run.h]().
*/
int minlevel = 4;
int maxlevel = 9;
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
We consider two space-curves, $\mathcal{C}_1(\xi,t)$ and $\mathcal{C}_2(\xi,t)$,
is parametrized as function of $\theta(\xi,t)$. At time $t=0$,
$$
\begin{aligned}
x_i &= R\cos(\theta), \quad
\\
y_i &= R\sin(\theta), \quad
\\
z_i &= \pm d/2
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


The curves will be stored as `struct vortex_filament` which must be released at
the end of the simulation.
*/
event init (t = 0) {
  double dtheta = 2*pi/((double)nseg);
  double theta[nseg];
  double a1[nseg];
  double a2[nseg];
  double vol1[nseg];
  double vol2[nseg];
  coord C1[nseg];
  coord C2[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i;
    C1[i] = (coord) { R * cos(theta[i]), R * sin(theta[i]), -d};
    C2[i] = (coord) { R * cos(theta[i]), R * sin(theta[i]), -2*d};
    
    vol1[i] = pi * sq(a) * R * dtheta;
    vol2[i] = pi * sq(a) * R * dtheta;
    
    a1[i] = a;
    a2[i] = a;
  }

  // We store the space-curve in a structure 
  allocate_vortex_filament_members(&filament1, nseg);
  allocate_vortex_filament_members(&filament1_prev, nseg);
  memcpy(filament1.theta, theta, nseg * sizeof(double));
  memcpy(filament1.C, C1, nseg * sizeof(coord));
  memcpy(filament1.a, a1, nseg * sizeof(double));  
  memcpy(filament1.vol, vol1, nseg * sizeof(double));

  local_induced_velocity(filament1);

  allocate_vortex_filament_members(&filament2, nseg);
  allocate_vortex_filament_members(&filament2_prev, nseg);
  memcpy(filament2.theta, theta, nseg * sizeof(double));
  memcpy(filament2.C, C2, nseg * sizeof(coord));
  memcpy(filament2.a, a2, nseg * sizeof(double));  
  memcpy(filament2.vol, vol2, nseg * sizeof(double));
    
  local_induced_velocity(filament2);

  for (int j = 0; j < nseg; j++) {
    filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
    filament2.a[j] = sqrt(filament2.vol[j]/(pi*filament2.s[j]));
  }

  view (camera="iso");  
  draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);
  draw_tube_along_curve(filament2.nseg, filament2.C, filament2.a);
  save ("prescribed_curve.png");

  FILE *fp = fopen("curve1.txt", "w"); 
  fclose(fp);  

  fp = fopen("curve2.txt", "w"); 
  fclose(fp);  

  /**
  We also initialize the Cartesian grid close to the vortex filaments and 
  initialize the velocity field to zero.
  */ 
  scalar dmin[];
  for (int i = (maxlevel-minlevel-1); i >= 0; i--){
    foreach(){      
      struct vortex_filament params1;
      params1 = filament1;
      params1.pcar = (coord){x,y,z};
      double dmin1 = get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params1);

      struct vortex_filament params2;
      params2 = filament2;
      params2.pcar = (coord){x,y,z};
      double dmin2 = get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params2);

      dmin[] = 0;
      dmin[] = ( max(dmin1, dmin2) < (i+1)*filament1.a[0])*noise();    
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

  free_vortex_filament_members(&filament2);
  free_vortex_filament_members(&filament2_prev);
}

/** 
## A Simple Time Advancing Scheme

Time integration is done in two steps. First, we evaluate the (self-)induced
velocity at $\vec{x}_c$, then the mutually induced velocities
*/
event evaluate_velocity (i++) {  
  
  memcpy(filament1_prev.Utotal, filament1.Utotal, nseg * sizeof(coord));
  memcpy(filament2_prev.Utotal, filament2.Utotal, nseg * sizeof(coord));

  local_induced_velocity(filament1);
  local_induced_velocity(filament2);  
  for (int j = 0; j < nseg; j++) {
    filament1.Uauto[j] = nonlocal_induced_velocity(filament1.C[j], filament1);
    filament2.Uauto[j] = nonlocal_induced_velocity(filament2.C[j], filament2);

    filament1.Umutual[j] = nonlocal_induced_velocity(filament1.C[j], filament2);
    filament2.Umutual[j] = nonlocal_induced_velocity(filament2.C[j], filament1);
    foreach_dimension(){
      filament1.Utotal[j].x = filament1.Uauto[j].x + filament1.Umutual[j].x + filament1.Ulocal[j].x - Uinfty.x;
      filament2.Utotal[j].x = filament2.Uauto[j].x + filament2.Umutual[j].x + filament2.Ulocal[j].x - Uinfty.x;      
    }
  }
}

/**
Then, we advect the filament segments using an explicit Adams-Bashforth scheme
and update the core sizes
*/ 
event advance_filaments (i++, last) {      
  for (int j = 0; j < nseg; j++) {
    foreach_dimension(){
      filament1.C[j].x += dt*(3.*filament1.Utotal[j].x - filament1_prev.Utotal[j].x)/2.;
      filament2.C[j].x += dt*(3.*filament2.Utotal[j].x - filament2_prev.Utotal[j].x)/2.;
    }
  } 
  local_induced_velocity(filament1);  
  local_induced_velocity(filament2);  
  for (int j = 0; j < nseg; j++) {
    filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
    filament2.a[j] = sqrt(filament2.vol[j]/(pi*filament2.s[j]));
  } 
  dt = dtnext (dtmax);
}

/** 
## Outputs

### Displacement of the vortex filaments

The two vortex rings are characterized by this leapfrogging motion.

![Motion of the vortex rings](test_2vortex_rings1/2vortex_rings.mp4)

*/

event sequence (t += 0.20){
  view(camera="iso", fov=20);
  draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);  
  draw_tube_along_curve(filament2.nseg, filament2.C, filament2.a);  
  save("2vortex_rings.mp4");
}

/**

### The solution outside the vortex tube

Vortex filaments are dealt as Lagrangian particles. We may also use an Eulerian
grid to evaluate the induced velocity at some point $\vec{x}\in\mathbb{R}^3$.

<center>
  <table>
  <tr>
  <td><center>![](test_2vortex_rings1/axial_velocity_0.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_5.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_10.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_15.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_20.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_25.png){ width="100%" }</center></td>
  </tr>
  <tr>
  <td><center>![](test_2vortex_rings1/axial_velocity_30.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_35.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_40.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_45.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_50.png){ width="100%" }</center></td>
  <td><center>![](test_2vortex_rings1/axial_velocity_55.png){ width="100%" }</center></td>
  </tr>
  </table>
  <center>Axial velocity induced at the mid-plane taken at regular intervals </center>
</center>

*/

event slice (t += 5.0){  
  foreach(){
    coord p = {x,y,z};
    coord u1_BS = nonlocal_induced_velocity(p, filament1);
    coord u2_BS = nonlocal_induced_velocity(p, filament2);
    foreach_dimension(){
      u.x[] = u1_BS.x + u2_BS.x - Uinfty.x;
    }
  }
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, maxlevel, minlevel);

  {
    char filename2[100]; 
    sprintf(filename2, "axial_velocity_%g.png", t);
    
    view(camera="left");  
    squares ("u.z", linear = true, spread = -1, n={1,0,0}, min=-1.00, max=1.00, map = cool_warm); 
    draw_tube_along_curve(filament1.nseg, filament1.C, filament1.a);  
    draw_tube_along_curve(filament2.nseg, filament2.C, filament2.a);  
    save(filename2);  
  }
}

/** 
### Time evolution

We may also follow the position of $\mathcal{C}_1$ and $\mathcal{C}_2$ over
time, as well as the evolution of the core sizes
*/

event final (t += 0.2, t <= 60.0){
  if (pid() == 0){
    FILE *fp = fopen("curve1.txt", "a"); 
    write_filament_state(fp, &filament1);
    fclose(fp);

    fp = fopen("curve2.txt", "a"); 
    write_filament_state(fp, &filament2);
    fclose(fp);
  }
}

/** 



~~~gnuplot Time evolution of the axial coordinates for points initially located at $\theta=0$
set term pngcairo enhanced size 640,480 font ",12"
set output 'plot_axial_position.png'
set xlabel "time"
set ylabel "axial coordinate"
plot 'curve1.txt' u 1:5 every ::1::1 w lp title "ring 1", \
     'curve2.txt' u 1:5 every ::1::1 w lp title "ring 2"
~~~

~~~gnuplot Time evolution of the radial coordinates for points initially located at $\theta=0$
set term pngcairo enhanced size 640,480 font ",12"
set output 'plot_radial_position.png'
set xlabel "time"
set ylabel "radial coordinate"
plot 'curve1.txt' u 1:(sqrt($3*$3 + $4*$4)) every ::1::1 w lp title "ring 1", \
     'curve2.txt' u 1:(sqrt($3*$3 + $4*$4)) every ::1::1 w lp title "ring 2"
~~~

~~~gnuplot Time evolution of the core size $a$ for points initially located at $\theta=0$
set term pngcairo enhanced size 640,480 font ",12"
set output 'plot_core.png'
set xlabel "time"
set ylabel "core size"
plot 'curve1.txt' u 1:18 every ::1::1 w lp title "ring 1", \
     'curve2.txt' u 1:18 every ::1::1 w lp title "ring 2"
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