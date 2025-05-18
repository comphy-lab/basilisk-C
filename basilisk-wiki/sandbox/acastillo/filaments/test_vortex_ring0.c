/**
# Motion of a thin vortex ring (using the vortex filament framework)

In this example, we consider the motion of a vortex ring of radius $R$ and core
size $a$ using the [vortex filament framework](../input_fields/filaments.h).

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

- In this test, we verify that the translational velocity $u_z$ is independent
  of the discretization:

~~~pythonplot Contributions to the total axial velocity from different sources as function of $n_{seg}$
import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()
data = np.loadtxt('curve.txt', delimiter=' ', usecols=[0,2,3,4])
plt.plot(data[:,0], data[:,1], ls='--', lw=1.4, label='local contributions');
plt.plot(data[:,0], data[:,2], ls='-.', lw=1.6, label='non-local contributions');
plt.plot(data[:,0], data[:,3], ls='-',  lw=2,   label='total velocity');
plt.legend(loc='best');
plt.xlabel(r'$n_{seg}$')
plt.ylabel(r'$U_z$')
plt.xlim([8,128])
plt.tight_layout()
plt.savefig('plot_Uz_vs_nseg.svg')
~~~

*/

#include "grid/octree.h"
#include "run.h"
#include "view.h"
#include "../input_fields/filaments.h"
#include "../input_fields/draw_filaments.h"
#include "biot-savart.h"

int minlevel = 4;
int maxlevel = 8;
vector u[];
int main(){
  L0 = 10.0;
  X0 = Y0 = Z0 = -L0 / 2;
  
  N = 1 << minlevel;  
  init_grid(N);  

  if (pid() == 0){
    FILE *fp = fopen("curve.txt", "w"); 
    fclose(fp);
  }  
  
  double R = 1.0;
  double a=0.01;
  struct vortex_filament filament1;
  coord Uinfty = {0., 0., 0.};  
  for (int nseg = 8; nseg < 128; nseg+=1){

    double dtheta = 2*pi/((double)nseg);
    double theta[nseg];
    double a1[nseg];
    double vol1[nseg];
    coord C1[nseg];

    // Define a curve 
    for (int j = 0; j < nseg; j++){
      theta[j] = dtheta * (double)j;
      C1[j] = (coord) { R * cos(theta[j]), R * sin(theta[j]), 0.};
      vol1[j] = pi * sq(a) * R * dtheta;
      a1[j] = a;
    } 
    
    // We store the space-curve in a structure   
    allocate_vortex_filament_members(&filament1, nseg);
    memcpy(filament1.theta, theta, nseg * sizeof(double));
    memcpy(filament1.C, C1, nseg * sizeof(coord));
    memcpy(filament1.a, a1, nseg * sizeof(double));  
    memcpy(filament1.vol, vol1, nseg * sizeof(double));

    local_induced_velocity(filament1);  
    for (int j = 0; j < nseg; j++) {
      filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
    }      
    local_induced_velocity(filament1);  
    for (int j = 0; j < nseg; j++) {
      filament1.Uauto[j] = nonlocal_induced_velocity(filament1.C[j], filament1);            
      foreach_dimension(){
        filament1.U[j].x = filament1.Uauto[j].x + filament1.Ulocal[j].x - Uinfty.x;      
      }
    }
    if (pid() == 0){
      FILE *fp = fopen("curve.txt", "a");
      for (int j = 0; j < nseg; j++) {
        fprintf(fp, "%d %d ", nseg, j);        
        fprintf(fp, "%g ", filament1.Ulocal[j].z);
        fprintf(fp, "%g ", filament1.Uauto[j].z);
        fprintf(fp, "%g \n", filament1.U[j].z);
      }
      fputs("\n", fp);  
      fclose(fp);
    }
    free_vortex_filament_members(&filament1);
  }
  
/**
- We also verify that the translational velocity $U_z$ has a log-dependency on
  the core size $a$, according to:
  $$
  U_z = \frac{\Gamma}{4\pi}\left(\log\left(\frac{8}{a}\right)-0.5\right)
  $$
  which is all contained inside the local contributions.

~~~pythonplot Contributions to the total axial velocity from different sources as function of $a$
import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()
data = np.loadtxt('curve2.txt', delimiter=' ', usecols=[2,3,4,5])
plt.semilogx(data[:,0], data[:,1], ls='--', lw=1.4, label='local contributions');
plt.semilogx(data[:,0], data[:,2], ls='-.', lw=1.6, label='non-local contributions');
plt.semilogx(data[:,0], data[:,3], ls='-',  lw=2,   label='total velocity');

def func(x):
  return (1/(4*np.pi))*(np.log(8/x)-0.5)
plt.semilogx(data[:,0], func(data[:,0]), '-k', label=r'$(\Gamma/4\pi)(\log(8/a)-0.5$')

plt.legend(loc='best');
plt.xlabel(r'$a$')
plt.ylabel(r'$U_z$')
plt.xlim([5e-4,5e-1])
plt.tight_layout()
plt.savefig('plot_Uz_vs_a.svg')
~~~
*/
  if (pid() == 0){
    FILE *fp = fopen("curve2.txt", "w"); 
    fclose(fp);
  }  

  for (int i = 1; i < 1000; i++){
    int nseg = 64;
    double a = 0.0005*(double)i;
    double dtheta = 2*pi/((double)nseg);
    double theta[nseg];
    double a1[nseg];
    double vol1[nseg];
    coord C1[nseg];

    // Define a curve 
    for (int j = 0; j < nseg; j++){
      theta[j] = dtheta * (double)j;
      C1[j] = (coord) { R * cos(theta[j]), R * sin(theta[j]), 0.};
      vol1[j] = pi * sq(a) * R * dtheta;
      a1[j] = a;
    } 
    
    // We store the space-curve in a structure   
    allocate_vortex_filament_members(&filament1, nseg);
    memcpy(filament1.theta, theta, nseg * sizeof(double));
    memcpy(filament1.C, C1, nseg * sizeof(coord));
    memcpy(filament1.a, a1, nseg * sizeof(double));  
    memcpy(filament1.vol, vol1, nseg * sizeof(double));
    
    local_induced_velocity(filament1);  
    for (int j = 0; j < nseg; j++) {
      filament1.a[j] = sqrt(filament1.vol[j]/(pi*filament1.s[j]));
    }      
    local_induced_velocity(filament1);  
    for (int j = 0; j < nseg; j++) {
      filament1.Uauto[j] = nonlocal_induced_velocity(filament1.C[j], filament1);            
      foreach_dimension(){
        filament1.U[j].x = filament1.Uauto[j].x + filament1.Ulocal[j].x - Uinfty.x;      
      }
    }
    if (pid() == 0){
      FILE *fp = fopen("curve2.txt", "a");
      for (int j = 0; j < nseg; j++) {
        fprintf(fp, "%d ", i);        
        fprintf(fp, "%d ", j);        
        fprintf(fp, "%g ", filament1.a[j]);        
        fprintf(fp, "%g ", filament1.Ulocal[j].z);
        fprintf(fp, "%g ", filament1.Uauto[j].z);
        fprintf(fp, "%g \n", filament1.U[j].z);
      }
      fputs("\n", fp);  
      fclose(fp);
    }
    free_vortex_filament_members(&filament1);
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