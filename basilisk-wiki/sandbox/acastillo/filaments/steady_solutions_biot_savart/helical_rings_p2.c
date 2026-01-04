/**
# Quasi-steady solutions for a vortex filament in shaped like a torus unknot, $\mathcal{U}_{2,1}$

Vortex rings and vortex helices are well-known steady solutions of the
Biot-Savart law. Vortex filaments in the shape of a torus knots and unknots may
also have steady solutions. In this example, we consider a vortex filament in
the shape of a torus unknot, $\mathcal{U}_{2,1}$, which is can be seen as a
helical coil wrapped around a torus with a single poloidal revolution and two
toroidal revolutions. We can also think of this curve as a thin vortex ring
superimposed with a Kelvin wave of azimuthal wavenumber $m=2$ and
finite-amplitude. Here, we consider 3 cases with different amplitudes of the
Kelvin wave, $A=0.1R$, $A=0.2R$ and $A=0.3R$, where $R$ is the radius of the
undeformed vortex ring. The shape of the vortex filaments are shown below.

![Steady vortex filaments in the shape of a torus unknot, $\mathcal{U}_{2,1}$](./helical_rings_p2/prescribed_curve1.png){width=75%}

The vortex filaments are initialized from a file containing the coordinates of
the filament centerline, the core size and the volume of each segment. 

The vortices move as material lines in the fluid according to:
$$
\begin{aligned}
\frac{d\vec{x}_c}{dt} = \vec{U}_{ind} + \vec{U}_{\infty}
\end{aligned}
$$
where $\vec{x}_c$ is the position vector of vortex filament, $\vec{U}_{ind}$ the
velocity induced by the vortex filament and $\vec{U}_{\infty}$ an external
velocity field. The induced velocity $\vec{U}_{ind}$ is given by the 
[Biot-Savart](biot-savart.h) law. 

In this example, the three vortex filaments are **independent** but are shown
together for visualization purposes. The vortex filaments rotate with a constant
angular velocity around the vertical axis and translate with a constant
velocity along the vertical axis. The motion of the vortex filaments is shown
below.

<div style="text-align: center;">
  <video controls style="width: 75%; max-width: 100%; height: auto; margin: 0 auto;">
    <source src="./helical_rings_p2/top_view.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>
</div>

<div style="text-align: center;">
  <video controls style="width: 75%; max-width: 100%; height: auto; margin: 0 auto;">
    <source src="./helical_rings_p2/side_view.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>
</div>
*/

#include "grid/octree.h"
#include "run.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"
#include "acastillo/input_fields/draw_filaments.h"
#include "acastillo/filaments/biot-savart.h"

/**
 In this example, we use $a$ as the characteristic length scale and rescale all
 other lengths by $1/a$. This ways, the velocity induced by a circular vortex ring of radius $R$ and circulation $\Gamma=2\pi$ remains of order unity.
*/

#define SCALE (0.05)

int nseg ;
double R = 1.0/SCALE;
double a = 0.05;
double dtmax = 0.05;
double tend = 2000.0;
struct vortex_filament helical_ring1;
struct vortex_filament helical_ring2;
struct vortex_filament helical_ring3;

coord Uinfty = {0., 0., 0.};

/**
 The main time loop is defined in [run.h]().
*/
int minlevel = 4;
int maxlevel = 4;
vector u[];
int main(){
  L0 = 9.0/SCALE;
  X0 = Y0 = Z0 = -L0 / 2;
  
  N = 1 << minlevel;  
  init_grid(N);  
  run();
}

/**
## Initial conditions 
We read the vortex filaments from a file. The files must have the following
format:
```
nseg
a
phi_0  x_0  y_0  z_0  dl_0
phi_1  x_1  y_1  z_1  dl_1
...
phi_(nseg-1)  x_(nseg-1)  y_(nseg-1)  z_(nseg-1)  dl_(nseg-1)
```
where `nseg` is the number of segments, `a` the core size (assumed constant),
`phi_i` the azimuthal angle of the segment, `(x_i,y_i,z_i)` the coordinates of
the center of the segment and `dl_i` the length of the segment. 
Viewed from the top, the vortex filaments look like an ellipse, but are actually
3D curves. 

~~~pythonplot Evolution of the axial coordinate of the filaments 
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
ax = plt.figure()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)
data = np.loadtxt('./helical_rings_p2.filament1.txt', delimiter=' ', usecols=[0,1,2,3], skiprows=2)
ax1.plot(data[:,1], data[:,2], label=r'$A=0.1R$')
ax2.plot(data[:,3], data[:,2], label=r'$A=0.1R$')

data = np.loadtxt('./helical_rings_p2.filament2.txt', delimiter=' ', usecols=[0,1,2,3], skiprows=2)
ax1.plot(data[:,1], data[:,2], label=r'$A=0.2R$')
ax2.plot(data[:,3], data[:,2], label=r'$A=0.2R$')

data = np.loadtxt('./helical_rings_p2.filament3.txt', delimiter=' ', usecols=[0,1,2,3], skiprows=2);
ax1.plot(data[:,1], data[:,2], label=r'$A=0.3R$');
ax2.plot(data[:,3], data[:,2], label=r'$A=0.3R$')

phi = np.linspace(0,2*np.pi)
ax1.plot(np.cos(phi), np.sin(phi), '--', color='grey', lw=0.5)

ax1.legend(loc='best')
ax1.axis('image')
ax1.set_xlabel(r'Coordinate $x$')
ax1.set_ylabel(r'Coordinate $y$')

ax2.axis('image')
ax2.set_xlabel(r'Coordinate $z$')
plt.tight_layout()
plt.savefig('plot_init.svg')
~~~

*/

void read_filament_from_file(const char *filename, struct vortex_filament *filament, double R=1.0, double xshift=0.0, double yshift=0.0, double zshift=0.0, double scale=1.0) {
  // Debug: Print current working directory and list files
  fprintf(stderr, "Attempting to open: %s\n", filename);
  fprintf(stderr, "Current working directory contents:\n");
  system("pwd > log2 2>&1");
  system("ls -la >> log2 2>&1");

  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Error: could not open %s\n", filename);
    fprintf(stderr, "Errno: %s\n", strerror(errno));
    exit(1);
  }  

  int nseg;
  double a;
  fscanf(fp, "%d", &nseg);
  fscanf(fp, "%lf", &a);
  
  double phi[nseg];
  double a1[nseg];
  double vol1[nseg];
  coord C1[nseg];

  for (int i = 0; i < nseg; i++) {
    double x, y, z, dl;
    if (fscanf(fp, "%lf %lf %lf %lf %lf", &phi[i], &x, &y, &z, &dl) != 5) {
      fprintf(stderr, "Error: invalid format in %s at line %d\n", filename, i+2);
      exit(1);
    }
    C1[i] = (coord){x * scale + xshift, y * scale + yshift, z * scale + zshift};
    a1[i] = a * scale;
    vol1[i] = pi * sq(a) * dl * cube(scale);
  }
  fclose(fp);

  allocate_vortex_filament_members(filament, nseg);
  memcpy(filament->phi, phi, nseg * sizeof(double));
  memcpy(filament->C, C1, nseg * sizeof(coord));
  memcpy(filament->a, a1, nseg * sizeof(double));  
  memcpy(filament->vol, vol1, nseg * sizeof(double));
}

event init (t = 0) {

  read_filament_from_file("./helical_rings_p2.filament1.txt", &helical_ring1, R, xshift=-L0/3, zshift=-L0/2, scale=1/SCALE);
  read_filament_from_file("./helical_rings_p2.filament2.txt", &helical_ring2, R, xshift= 0.0,  zshift=-L0/2, scale=1/SCALE);
  read_filament_from_file("./helical_rings_p2.filament3.txt", &helical_ring3, R, xshift= L0/3, zshift=-L0/2, scale=1/SCALE);

  local_induced_velocity(&helical_ring1);  
  local_induced_velocity(&helical_ring2);
  local_induced_velocity(&helical_ring3);
  
  view (camera="front", fov=6, width=2560, height=800);
  draw_tube_along_curve(helical_ring1.nseg, helical_ring1.C, helical_ring1.a, scale=1);
  draw_tube_along_curve(helical_ring2.nseg, helical_ring2.C, helical_ring2.a, scale=1);
  draw_tube_along_curve(helical_ring3.nseg, helical_ring3.C, helical_ring3.a, scale=1);
  save ("prescribed_curve1.png");
  clear();
  
  view (camera="top", fov=25, width=1920, height=1920);
  draw_tube_along_curve(helical_ring1.nseg, helical_ring1.C, helical_ring1.a, scale=1);
  draw_tube_along_curve(helical_ring2.nseg, helical_ring2.C, helical_ring2.a, scale=1);
  draw_tube_along_curve(helical_ring3.nseg, helical_ring3.C, helical_ring3.a, scale=1);
  save ("prescribed_curve2.png");
  clear();
}

/**
We release the vortex filaments at the end of the simulation. 
*/
event finalize(t = end){
  free_vortex_filament_members(&helical_ring1);
  free_vortex_filament_members(&helical_ring2);
  free_vortex_filament_members(&helical_ring3);
}

/** 
## A Simple Time Advancing Scheme

Time integration is done in two steps. First, we evaluate the (self-)induced
velocity at $\vec{x}_c$
*/
event evaluate_velocity (i++) {

  memcpy(helical_ring1.Uprev, helical_ring1.U, helical_ring1.nseg * sizeof(coord));
  local_induced_velocity(&helical_ring1);
  for (int j = 0; j < helical_ring1.nseg; j++) {
    helical_ring1.Uauto[j] = nonlocal_induced_velocity(helical_ring1.C[j], &helical_ring1);
    foreach_dimension(){
      helical_ring1.U[j].x = helical_ring1.Uauto[j].x + helical_ring1.Ulocal[j].x - Uinfty.x;      
    }
  }

  memcpy(helical_ring2.Uprev, helical_ring2.U, helical_ring2.nseg * sizeof(coord));
  local_induced_velocity(&helical_ring2);
  for (int j = 0; j < helical_ring2.nseg; j++) {
    helical_ring2.Uauto[j] = nonlocal_induced_velocity(helical_ring2.C[j], &helical_ring2);
    foreach_dimension(){
      helical_ring2.U[j].x = helical_ring2.Uauto[j].x + helical_ring2.Ulocal[j].x - Uinfty.x;      
    }
  }

  memcpy(helical_ring3.Uprev, helical_ring3.U, helical_ring3.nseg * sizeof(coord));
  local_induced_velocity(&helical_ring3);
  for (int j = 0; j < helical_ring3.nseg; j++) {
    helical_ring3.Uauto[j] = nonlocal_induced_velocity(helical_ring3.C[j], &helical_ring3);
    foreach_dimension(){
      helical_ring3.U[j].x = helical_ring3.Uauto[j].x + helical_ring3.Ulocal[j].x - Uinfty.x;      
    }
  }
}

/**
Then, we advect the filament segments using an explicit Adams-Bashforth scheme
*/ 
event advance_filaments (i++) {      
  for (int j = 0; j < helical_ring1.nseg; j++) {
    foreach_dimension(){
      helical_ring1.C[j].x += dt*(3.*helical_ring1.U[j].x - helical_ring1.Uprev[j].x)/2.;
    }
  } 

  for (int j = 0; j < helical_ring2.nseg; j++) {
    foreach_dimension(){
      helical_ring2.C[j].x += dt*(3.*helical_ring2.U[j].x - helical_ring2.Uprev[j].x)/2.;
    }
  }
  
  for (int j = 0; j < helical_ring3.nseg; j++) {
    foreach_dimension(){
      helical_ring3.C[j].x += dt*(3.*helical_ring3.U[j].x - helical_ring3.Uprev[j].x)/2.;
    }
  }
}

/** 
Finally, we advance by one time-step
*/
event advance_filaments (i++, last) {      
  dt = dtnext (dtmax);
}

/** 
## Outputs

- Viewed from the top, we can see the elliptical ring seems to be rotating with 
a constant angular velocity around the vertical axis and translating with a
constant velocity along the vertical axis. However, these velocities depend on
the amplitude of the Kelvin wave.

~~~pythonplot Evolution of the x- and y-coordinates of the first filament segment
import numpy as np
import matplotlib.pyplot as plt
	
a = 0.05
n = 128
ax = plt.figure()
data = np.loadtxt('curve1.txt', delimiter=' ', usecols=[0,2,3,4])
data[:,1] = data[:,1] + 9/a/3
r = np.sqrt(data[:,1]**2 + data[:,2]**2)
plt.plot(data[::n,3]*a, r[::n]*a, label='$A=0.1R$')

data = np.loadtxt('curve2.txt', delimiter=' ', usecols=[0,2,3,4])
r = np.sqrt(data[:,1]**2 + data[:,2]**2)
plt.plot(data[::n,3]*a, r[::n]*a, label='$A=0.2R$')

data = np.loadtxt('curve3.txt', delimiter=' ', usecols=[0,2,3,4])
data[:,1] = data[:,1] - 9/a/3
r = np.sqrt(data[:,1]**2 + data[:,2]**2)
plt.plot(data[::n,3]*a, r[::n]*a, label='$A=0.3R$')


plt.xlabel(r'Coordinate $z/R$')
plt.ylabel(r'Coordinate $r/R$')
plt.legend()
plt.tight_layout()
plt.savefig('plot_xy_vs_t.svg')
~~~

*/


event movie (t += 5.00){
  view (camera="front", fov=6, width=2560, height=800);
  draw_tube_along_curve(helical_ring1.nseg, helical_ring1.C, helical_ring1.a, scale=1);
  draw_tube_along_curve(helical_ring2.nseg, helical_ring2.C, helical_ring2.a, scale=1);
  draw_tube_along_curve(helical_ring3.nseg, helical_ring3.C, helical_ring3.a, scale=1);
  save("top_view.mp4");

  view (camera="top", fov=25, width=1920, height=1920);
  draw_tube_along_curve(helical_ring1.nseg, helical_ring1.C, helical_ring1.a, scale=1);
  draw_tube_along_curve(helical_ring2.nseg, helical_ring2.C, helical_ring2.a, scale=1);
  draw_tube_along_curve(helical_ring3.nseg, helical_ring3.C, helical_ring3.a, scale=1);
  save("side_view.mp4");
}

event final (t += 5.00, t <= tend){
  if (pid() == 0){
    FILE *fp;
    fp = fopen("curve1.txt", "a"); 
    write_filament_state(fp, &helical_ring1); 
    fclose(fp);

    fp = fopen("curve2.txt", "a"); 
    write_filament_state(fp, &helical_ring2); 
    fclose(fp);

    fp = fopen("curve3.txt", "a"); 
    write_filament_state(fp, &helical_ring3); 
    fclose(fp);
  }
}


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
