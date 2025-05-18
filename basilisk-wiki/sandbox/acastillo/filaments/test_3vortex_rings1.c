/**
# Motion of three thin vortex rings (using the vortex filament framework)

In this example, we consider the motion of three vortex rings of radius $R$ and
core size $a$ using the [vortex filament framework](../input_fields/filaments.h).

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
[Biot-Savart](biot-savart.h) law. 

In this example, we're going to create a simple solver for the equation above
and compute the external field in the Eulerian grid. We're also going to compare
the motion of the vortex rings to the case of a single vortex ring. To this end,
we'll use the translational velocity of the circular vortex ring as
$\vec{U}_{\infty}$. The three vortex rings are characterized by a more complex
leapfrogging motion.

![Motion of the vortex rings](test_3vortex_rings1/3vortex_rings.mp4)

*/

#include "grid/octree.h"
#include "run.h"
#include "view.h"
#include "../input_fields/filaments.h"
#include "../input_fields/draw_filaments.h"
#include "biot-savart.h"

/**
In this example, $R_i=1.0$, $a_i=0.05$ and the vortex rings are discretized
 into $n_{seg}=64$ filaments. Also, the separation distance $d=2.0$ 
*/
int nseg = 64;
double R = 1.0;
double d = 2.0;
double a = 0.05;
double dtmax = 0.001;
struct vortex_filament ring1;
struct vortex_filament ring2;
struct vortex_filament ring3;

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
We consider three space-curves, $\mathcal{C}_i(\xi,t) (i=1,2,3)$, parametrized
as function of $\theta(\xi,t)$. At time $t=0$,
$$
\begin{aligned}
x_i &= R\cos(\theta), \quad
\\
y_i &= R\sin(\theta), \quad
\\
z_i &= z_0 + i d
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
  double a3[nseg];
  double vol1[nseg];
  double vol2[nseg];
  double vol3[nseg];
  coord C1[nseg];
  coord C2[nseg];
  coord C3[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i;
    C1[i] = (coord) { R * cos(theta[i]), R * sin(theta[i]), -2*d};
    C2[i] = (coord) { R * cos(theta[i]), R * sin(theta[i]), -1*d};
    C3[i] = (coord) { R * cos(theta[i]), R * sin(theta[i]),  0.};
    
    vol1[i] = pi * sq(a) * R * dtheta;
    vol2[i] = pi * sq(a) * R * dtheta;
    vol3[i] = pi * sq(a) * R * dtheta;
    
    a1[i] = a;
    a2[i] = a;
    a3[i] = a;
  }

  // We store the space-curve in a structure 
  allocate_vortex_filament_members(&ring1, nseg);
  memcpy(ring1.theta, theta, nseg * sizeof(double));
  memcpy(ring1.C, C1, nseg * sizeof(coord));
  memcpy(ring1.a, a1, nseg * sizeof(double));  
  memcpy(ring1.vol, vol1, nseg * sizeof(double));
  local_induced_velocity(ring1);

  allocate_vortex_filament_members(&ring2, nseg);
  memcpy(ring2.theta, theta, nseg * sizeof(double));
  memcpy(ring2.C, C2, nseg * sizeof(coord));
  memcpy(ring2.a, a2, nseg * sizeof(double));  
  memcpy(ring2.vol, vol2, nseg * sizeof(double));    
  local_induced_velocity(ring2);

  allocate_vortex_filament_members(&ring3, nseg);
  memcpy(ring3.theta, theta, nseg * sizeof(double));
  memcpy(ring3.C, C3, nseg * sizeof(coord));
  memcpy(ring3.a, a3, nseg * sizeof(double));  
  memcpy(ring3.vol, vol3, nseg * sizeof(double));    
  local_induced_velocity(ring3);

  for (int j = 0; j < nseg; j++) {
    ring1.a[j] = sqrt(ring1.vol[j]/(pi*ring1.s[j]));
    ring2.a[j] = sqrt(ring2.vol[j]/(pi*ring2.s[j]));
    ring3.a[j] = sqrt(ring3.vol[j]/(pi*ring3.s[j]));
  }

  view (camera="iso");  
  draw_tube_along_curve(ring1.nseg, ring1.C, ring1.a);
  draw_tube_along_curve(ring2.nseg, ring2.C, ring2.a);
  draw_tube_along_curve(ring3.nseg, ring3.C, ring3.a);
  save ("prescribed_curve.png");

  FILE *fp = fopen("curve1.txt", "w"); 
  fclose(fp);  

  fp = fopen("curve2.txt", "w"); 
  fclose(fp);  

  fp = fopen("curve3.txt", "w"); 
  fclose(fp);  

  /**
  We also initialize the Cartesian grid close to the vortex filaments and 
  initialize the velocity field to zero.
  */ 
  scalar dmin[];
  for (int i = (maxlevel-minlevel-1); i >= 0; i--){
    foreach(){      
      struct vortex_filament params1;
      params1 = ring1;
      params1.pcar = (coord){x,y,z};
      double dmin1 = get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params1);

      struct vortex_filament params2;
      params2 = ring2;
      params2.pcar = (coord){x,y,z};
      double dmin2 = get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params2);

      struct vortex_filament params3;
      params3 = ring3;
      params3.pcar = (coord){x,y,z};
      double dmin3 = get_min_distance(spatial_period=0, max_distance=4*L0, vortex=&params3);

      dmin[] = 0;
      dmin[] = ( max(max(dmin1, dmin2), dmin3) < (i+1)*ring1.a[0])*noise();    
    }    
    adapt_wavelet ((scalar*){dmin}, (double[]){1e-12}, maxlevel-i, minlevel);    
  }  
  foreach(){
    foreach_dimension(){
      u.x[] = 0.;
    }
  }
}

/**
We release the vortex filaments at the end of the simulation. 
*/
event finalize(t = end){
  free_vortex_filament_members(&ring1);
  free_vortex_filament_members(&ring2);
  free_vortex_filament_members(&ring3);
}

/** 
## A Simple Time Advancing Scheme

Time integration is done in two steps. First, we evaluate the (self-)induced
velocity at $\vec{x}_c$, then the mutually induced velocities
*/
event evaluate_velocity (i++) {  
  
  memcpy(ring1.Uprev, ring1.U, nseg * sizeof(coord));
  memcpy(ring2.Uprev, ring2.U, nseg * sizeof(coord));
  memcpy(ring3.Uprev, ring3.U, nseg * sizeof(coord));

  local_induced_velocity(ring1);
  local_induced_velocity(ring2);  
  local_induced_velocity(ring3);  
  for (int j = 0; j < nseg; j++) {
    ring1.Uauto[j] = nonlocal_induced_velocity(ring1.C[j], ring1);
    ring2.Uauto[j] = nonlocal_induced_velocity(ring2.C[j], ring2);
    ring3.Uauto[j] = nonlocal_induced_velocity(ring3.C[j], ring3);

    coord U12 = nonlocal_induced_velocity(ring1.C[j], ring2);
    coord U13 = nonlocal_induced_velocity(ring1.C[j], ring3);
    coord U21 = nonlocal_induced_velocity(ring2.C[j], ring1);
    coord U23 = nonlocal_induced_velocity(ring2.C[j], ring3);
    coord U31 = nonlocal_induced_velocity(ring3.C[j], ring1);
    coord U32 = nonlocal_induced_velocity(ring3.C[j], ring2);
    
    foreach_dimension(){
      ring1.Umutual[j].x = U12.x + U13.x;
      ring2.Umutual[j].x = U21.x + U23.x;
      ring3.Umutual[j].x = U31.x + U32.x;

      ring1.U[j].x = ring1.Uauto[j].x + ring1.Umutual[j].x + ring1.Ulocal[j].x - Uinfty.x;
      ring2.U[j].x = ring2.Uauto[j].x + ring2.Umutual[j].x + ring2.Ulocal[j].x - Uinfty.x;      
      ring3.U[j].x = ring3.Uauto[j].x + ring3.Umutual[j].x + ring3.Ulocal[j].x - Uinfty.x;      
    }
  }
}

/**
Then, we advect the filament segments using an explicit Adams-Bashforth scheme
*/ 
event advance_filaments (i++) {      
  for (int j = 0; j < nseg; j++) {
    foreach_dimension(){
      ring1.C[j].x += dt*(3.*ring1.U[j].x - ring1.Uprev[j].x)/2.;
      ring2.C[j].x += dt*(3.*ring2.U[j].x - ring2.Uprev[j].x)/2.;
      ring3.C[j].x += dt*(3.*ring3.U[j].x - ring3.Uprev[j].x)/2.;
    }
  } 
}

/** 
Finally, we compute the new arch-lenghts and update the core sizes to preserve 
the total vorticity
*/
event advance_filaments (i++, last) {      
  local_induced_velocity(ring1);  
  local_induced_velocity(ring2);
  local_induced_velocity(ring3);  
  for (int j = 0; j < nseg; j++) {
    ring1.a[j] = sqrt(ring1.vol[j]/(pi*ring1.s[j]));
    ring2.a[j] = sqrt(ring2.vol[j]/(pi*ring2.s[j]));
    ring3.a[j] = sqrt(ring3.vol[j]/(pi*ring3.s[j]));
  } 
  dt = dtnext (dtmax);
}

/** 
## Outputs

### Displacement of the vortex filaments
*/

event sequence (t += 0.20){
  view(camera="iso", fov=20);
  draw_tube_along_curve(ring1.nseg, ring1.C, ring1.a);  
  draw_tube_along_curve(ring2.nseg, ring2.C, ring2.a);  
  draw_tube_along_curve(ring3.nseg, ring3.C, ring3.a);  
  save("3vortex_rings.mp4");
}

/**

### The solution outside the vortex tube

Vortex filaments are treated as Lagrangian particles and may exist outside the
mesh. We may also evaluate the induced velocity at any point 
$\vec{x}\in\mathbb{R}^3$.
<center>
  <table>
  <tr>
  <td>![](test_3vortex_rings1/axial_velocity_10.png){ width="100%" }</td>
  <td>![](test_3vortex_rings1/axial_velocity_15.png){ width="100%" }</td>
  <td>![](test_3vortex_rings1/axial_velocity_20.png){ width="100%" }</td>
  </tr>
  <tr>
  <td>![](test_3vortex_rings1/axial_velocity_25.png){ width="100%" }</td>
  <td>![](test_3vortex_rings1/axial_velocity_30.png){ width="100%" }</td>
  <td>![](test_3vortex_rings1/axial_velocity_35.png){ width="100%" }</td>
  </tr>
  </table>
  <center>Axial velocity induced at the mid-plane taken at regular intervals </center>
</center>

*/

event slice (t += 5.0){  
  foreach(){
    coord p = {x,y,z};
    coord u1_BS = nonlocal_induced_velocity(p, ring1);
    coord u2_BS = nonlocal_induced_velocity(p, ring2);
    coord u3_BS = nonlocal_induced_velocity(p, ring3);
    foreach_dimension(){
      u.x[] = u1_BS.x + u2_BS.x + u3_BS.x - Uinfty.x;
    }
  }
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, maxlevel, minlevel);

  {
    char filename2[100]; 
    sprintf(filename2, "axial_velocity_%g.png", t);
    
    view(camera="left");  
    squares ("u.z", linear = true, spread = -1, n={1,0,0}, min=-1.00, max=1.00, map = cool_warm); 
    draw_tube_along_curve(ring1.nseg, ring1.C, ring1.a);  
    draw_tube_along_curve(ring2.nseg, ring2.C, ring2.a);  
    draw_tube_along_curve(ring3.nseg, ring3.C, ring3.a);  
    save(filename2);  
  }
}

/** 
### Time evolution

We may also follow the position of $\mathcal{C}_i$ and the induced velocities
over time, 
*/

event final (t += 0.2, t <= 75.0){
  if (pid() == 0){
    FILE *fp = fopen("curve1.txt", "a"); 
    write_filament_state(fp, &ring1);
    fclose(fp);

    fp = fopen("curve2.txt", "a"); 
    write_filament_state(fp, &ring2);
    fclose(fp);

    fp = fopen("curve3.txt", "a"); 
    write_filament_state(fp, &ring3);
    fclose(fp);
  }
}

/** 

- Leapfrogging can be seen by following the axial coordinates. The outer ring
  becomes slower, while the inner rings shoots forward. Together, the
  three rings move faster than a single ring.

~~~pythonplot Evolution of the axial coordinate
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
ax = plt.figure()
data1 = np.loadtxt('curve1.txt', delimiter=' ', usecols=[0,16])
data2 = np.loadtxt('curve2.txt', delimiter=' ', usecols=[0,16])
data3 = np.loadtxt('curve3.txt', delimiter=' ', usecols=[0,16])
plt.plot(data1[:,0], data1[:,1], label='Vortex 1')
plt.plot(data2[:,0], data2[:,1], label='Vortex 2')
plt.plot(data3[:,0], data3[:,1], label='Vortex 3')

plt.xlabel(r'Time $t$')
plt.ylabel(r'Translational velocity $U_z$')
plt.xlim([0,75])
plt.legend()
plt.tight_layout()
plt.savefig('plot_Uz_vs_t.svg')
~~~

~~~pythonplot Evolution of the axial coordinate
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
ax = plt.figure()
data1 = np.loadtxt('curve1.txt', delimiter=' ', usecols=[0,4])
data2 = np.loadtxt('curve2.txt', delimiter=' ', usecols=[0,4])
data3 = np.loadtxt('curve3.txt', delimiter=' ', usecols=[0,4])
plt.plot(data1[:,0], data1[:,1], label='Vortex 1')
plt.plot(data2[:,0], data2[:,1], label='Vortex 2')
plt.plot(data3[:,0], data3[:,1], label='Vortex 3')

plt.xlabel(r'Time $t$')
plt.ylabel(r'Axial coordinate $z$')
plt.xlim([0,75])
plt.legend()
plt.tight_layout()
plt.savefig('plot_z_vs_t.svg')
~~~

- Leapfrogging is also reflected on the radial coordinate

~~~pythonplot Evolution of the radial coordinate
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
ax = plt.figure()
data1 = np.loadtxt('curve1.txt', delimiter=' ', usecols=[0,2,3])
data2 = np.loadtxt('curve2.txt', delimiter=' ', usecols=[0,2,3])
data3 = np.loadtxt('curve3.txt', delimiter=' ', usecols=[0,2,3])
plt.plot(data1[:,0], np.sqrt(data1[:,1]**2+data1[:,2]**2), label='Vortex 1')
plt.plot(data2[:,0], np.sqrt(data2[:,1]**2+data2[:,2]**2), label='Vortex 2')
plt.plot(data3[:,0], np.sqrt(data3[:,1]**2+data3[:,2]**2), label='Vortex 3')

plt.xlabel(r'Time $t$')
plt.ylabel(r'Radial coordinate $r$')
plt.xlim([0,75])
plt.legend()
plt.tight_layout()
plt.savefig('plot_r_vs_t.svg')
~~~

- Strecthing of the cores is also visible,

~~~pythonplot Evolution of the core size
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
ax = plt.figure()
data1 = np.loadtxt('curve1.txt', delimiter=' ', usecols=[0,17])
data2 = np.loadtxt('curve2.txt', delimiter=' ', usecols=[0,17])
data3 = np.loadtxt('curve3.txt', delimiter=' ', usecols=[0,17])
plt.plot(data1[:,0], data1[:,1], label='Vortex 1')
plt.plot(data2[:,0], data2[:,1], label='Vortex 2')
plt.plot(data3[:,0], data3[:,1], label='Vortex 3')

plt.xlabel(r'Time $t$')
plt.ylabel(r'Core size $a$')
plt.xlim([0,75])
plt.legend()
plt.tight_layout()
plt.savefig('plot_a_vs_t.svg')
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