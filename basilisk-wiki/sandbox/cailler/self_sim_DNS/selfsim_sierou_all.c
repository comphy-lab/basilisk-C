/** 

# Axisymmetric Self-Similar DNS of a Recoiling Cone (no dipolar flow)

Reproduction in the self-similar domain of the [Sierou \& Lister, (2004)](#sierou2004) 
liquid cone recoil driven by surface tension **without** any dipolar far-field flow 
$(\overline{\mu}_0 = 0)$, for the following initial cone angles:

$$
\theta_0 \text{ [rad] } = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2
$$

This is an extension to the **axisymmetric** case of 
[this simulation](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_keller_all.c). 
Therefore, we strongly suggest to the reader to take a look to the following 
*Wiki* pages:

* The simpler 2D--case (*Keller \& Miksis* problem) is exposed 
[here](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c);
* The self-similar solver used in this file is explained 
[there](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README).

~~~pythonplot Self-similar solver *VS* Sierou \& Lister
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
sierou_0p2rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_0.2rad.dat", 
                             unpack=False)
sierou_0p4rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_0.4rad.dat", 
                             unpack=False)
sierou_0p6rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_0.6rad.dat", 
                             unpack=False)
sierou_0p8rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_0.8rad.dat", 
                             unpack=False)
sierou_1p0rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_1.0rad.dat", 
                             unpack=False)
sierou_1p2rad = np.loadtxt ("../../sierou_lister/sierou_data/raw_sierou_fig2_1.2rad.dat", 
                             unpack=False)

path_shape_0p2rad = "../../sierou_lister/sierou_data/basilisk_selfsim_sierou_shape_th11.4592_t5_N11.dat"
path_shape_0p4rad = "selfsim_sierou_shape_th22.9183_t5_N9.dat"
path_shape_0p6rad = "selfsim_sierou_shape_th34.3775_t5_N9.dat"
path_shape_0p8rad = "selfsim_sierou_shape_th45.8366_t10_N8.dat"
path_shape_1p0rad = "selfsim_sierou_shape_th57.2958_t10_N8.dat"
path_shape_1p2rad = "selfsim_sierou_shape_th68.7549_t10_N8.dat"

facets_0p2rad = np.loadtxt (path_shape_0p2rad, unpack=False)
facets_0p4rad = np.loadtxt (path_shape_0p4rad, unpack=False)
facets_0p6rad = np.loadtxt (path_shape_0p6rad, unpack=False)
facets_0p8rad = np.loadtxt (path_shape_0p8rad, unpack=False)
facets_1p0rad = np.loadtxt (path_shape_1p0rad, unpack=False)
facets_1p2rad = np.loadtxt (path_shape_1p2rad, unpack=False)

# # Get segments of the output_facets() function of Basilisk C:
N_seg_0p2rad = int (0.5*facets_0p2rad.shape[0])
N_seg_0p4rad  = int (0.5*facets_0p4rad.shape[0])
N_seg_0p6rad  = int (0.5*facets_0p6rad.shape[0])
N_seg_0p8rad  = int (0.5*facets_0p8rad.shape[0])
N_seg_1p0rad = int (0.5*facets_1p0rad.shape[0])
N_seg_1p2rad = int (0.5*facets_1p2rad.shape[0])

segments_0p2rad  = np.split (facets_0p2rad , indices_or_sections=N_seg_0p2rad )
segments_0p4rad  = np.split (facets_0p4rad , indices_or_sections=N_seg_0p4rad )
segments_0p6rad  = np.split (facets_0p6rad , indices_or_sections=N_seg_0p6rad )
segments_0p8rad  = np.split (facets_0p8rad , indices_or_sections=N_seg_0p8rad )
segments_1p0rad  = np.split (facets_1p0rad , indices_or_sections=N_seg_1p0rad )
segments_1p2rad  = np.split (facets_1p2rad , indices_or_sections=N_seg_1p2rad )


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 12
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,5)
ax.set_ylim(0,5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')

# Basilisk results with:
# ----------------------
# adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, maxlvl=9, minlvl=6);
# unrefine(x < -7./24.*SIZE); 
# unrefine( ((f[] == 1.) || (f[] == 0)) && (y < 10.*Delta) ); 
# t_end = 1e1

for segment in segments_0p2rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[0]) 
for segment in segments_0p4rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[2]) 
for segment in segments_0p6rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[4]) 
for segment in segments_0p8rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[6]) 
for segment in segments_1p0rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[8]) 
for segment in segments_1p2rad:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[10]) 

# Sierou & Lister (2004) results:
# -------------------------------
ax.plot(sierou_0p2rad[:,0], sierou_0p2rad[:,1], '--', lw=0.8, color="darkred")
ax.plot(sierou_0p4rad[:,0], sierou_0p4rad[:,1], '--', lw=0.8, color="darkred")
ax.plot(sierou_0p6rad[:,0], sierou_0p6rad[:,1], '--', lw=0.8, color="darkred")
ax.plot(sierou_0p8rad[:,0], sierou_0p8rad[:,1], '--', lw=0.8, color="darkred")
ax.plot(sierou_1p0rad[:,0], sierou_1p0rad[:,1], '--', lw=0.8, color="darkred")
ax.plot(sierou_1p2rad[:,0], sierou_1p2rad[:,1], '--', lw=0.8, color="darkred", 
        label=r'Sierou & Lister')

ax.legend(frameon=False, loc=4)

plt.text(4.2, 1., r'$\theta_0 = 0.2$', fontsize=14,
         ha='center', va='center', rotation=11.4592, color=mapcolors[0])
plt.text(4., 1.9, r'$\theta_0 = 0.4$', fontsize=14,
         ha='center', va='center', rotation=22.9183, color=mapcolors[2])
plt.text(3.8, 2.8, r'$\theta_0 = 0.6$', fontsize=14,
         ha='center', va='center', rotation=34.3775, color=mapcolors[4])
plt.text(3.2, 3.55, r'$\theta_0 = 0.8$', fontsize=14,
         ha='center', va='center', rotation=45.8366, color=mapcolors[6])
plt.text(2.5, 4.15, r'$\theta_0 = 1.0$', fontsize=14,
         ha='center', va='center', rotation=57.2958, color=mapcolors[8])
plt.text(1.5, 4.3, r'$\theta_0 = 1.2$', fontsize=14,
         ha='center', va='center', rotation=68.7549, color=mapcolors[10])


plt.savefig('selfsim_solver_sierou_fig2_repro_N8-9.svg') 
~~~

As shown on the above figure: 

* as in the 2D--case, all the numerical results with the `Basilisk` self-similar 
solver match those from [Sierou \& Lister, (2004)](#sierou2004), which 
**validates the self-similar solver for the 3D--AXI configuration**;
* one can observe that surface wave amplitudes are lower than in the 2D--case, 
because of the wave energy spread on a $\xi$--increasing surface in 3D: 
surface waves amplitudes scale as $\mathcal{O}\left(\overline{R}^{-5}\right)$ 
in 3D--AXI, against $\mathcal{O}\left(\overline{R}^{-7/2}\right)$ in 2D; 
* no visible effect of the viscosity concerning the shapes of the interfaces, 
compared to the potential results of the reference paper, as in 2D. 

*/




/** 
## Code

### General Parameters 

The numerical configuration proceeds from previous preliminary studies done 
both in the [physical domain](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c), 
where it appeared that a box of size $[12 \times 12]$ would be relevant to 
avoid border effects and visualizing a conical shape interface at boundaries, 
and in the [self-similar domain](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_keller_all.c) 
for the 2D--case, in order to know the optimized resolution and running time 
parameters for converging to a steady (scale invariant) state with the 
self-similar solver. 

*/


#define LEVEL 7 
#define SIZE 12e0
#define X_OFFSET (1./3.)*SIZE

#include "axi.h"
#include "selfsim_centered_keller.h"
#include "contact.h"
#include "selfsim_two-phase_keller.h"
#include "tension.h"


// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 

scalar f0[]; // auxiliary volume fraction field used for BC 
vector n_front[]; // normal vector of the interface
scalar alpha_front[]; // intercept linked to n_front[] during VoF-reconstruction 

/** We test for the following initial cone angles:
$\theta_0 \text{ [rad] } = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2.$ 

For large angles, we do not need to use a high level of refinement, since 
capillary waves are almost non existent, therefore $N_{max} = 8$ is sufficient. 

For $\theta_0 = 0.4$ or $0.6 \text{ [rad]}$, results are slightly improved 
by increasing by one level of refinement. In that case, we cut in half the 
total running time to speed up the results (at $\tau_{end} = 5$, the steady 
state is already reached in the self-similar space).

Due to the very low cone angle for $\theta_0 = 0.2 \, \text{rad} \approx 11.5$°, 
along with the increased need for curvature computation, 
we need to go as far as $N_{max} = 11$ to avoid non convergence of the 
*Poisson* solver. On a *MacBook Pro M1*, the computation for this single 
angle lasted *47 minutes*... 
This is the reason why we do not run this case on the server as it would exceed 
the time limitation, and directly provide the results obtained in local for 
this level of refinement. 

However, the non convergence for $N_{max}$ lower than $11$ only gives *warnings* 
in the `log` file, and no direct computation errors: you are free to try with 
$\theta_0 = 0.2 \,\, / \,\, N_{max} = 9$ for instance. 
*/

int k; // counter
double theta_0[] = {0.4, 0.6, 0.8, 1.0, 1.2};
double lvl_max[] = {9, 9, 8, 8, 8};
double t_end[] = {5., 5., 10., 10., 10.};


/** 
### Boundary Conditions
*/
#include "selfsim_f_BC_keller_all.h"

/** 
Contrary to the 2D--case of [Keller \& Miksis, (1983)](#keller1983), 
a *capillary flow* does exist for the 3D--AXI case of [Sierou \& Lister, (2004)](#sierou2004): 
the azimuthal curvature of the cone is not zero, hence there is a *Laplace* 
pressure gradient generating this capillary flow. 

For the special case of a surface tension recoiling liquid cone with no external 
fluid and no forced dipolar far-field flow ($\overline{\mu}_0 = 0$), 
the authors have shown that in self-similar spherical-polar coordinates, 
the velocity far from the apex is simply: 

$$
\mathbf{\overline{u}}_{\infty} 
= \dfrac{\cot \theta_0}{\overline{R}^2} \, \mathbf{e}_{\overline{R}} 
= \dfrac{\cot \theta_0}{(\xi^2 + \eta^2)^{3/2}} \, \boldsymbol{\xi} 
$$

where:

$$
\overline{R} = \sqrt{\xi^2 + \eta^2}, \quad 
\tan \theta = \eta/\xi, \quad  
\xi = \left(\dfrac{\rho}{\sigma t^2} \right) \, z, \quad 
\eta = \left(\dfrac{\rho}{\sigma t^2} \right) \, r
$$

with the classical notations $\sigma$ for surface tension and $\rho$ for 
the liquid density.

As a result, **boundary conditions have to be modified** according to the 
below figure:

![Initial configuration for the 3D--AXI recoil of a cone without dipolar flow](BC_sierou_no_dipolar_flow.png){width="30%"}

It can also be noticed that *symmetry* conditions are applied onto the gas phase, 
since we do not expect any physical phenomenon at stake in this region. 

*/


  /* Normal vectors of the interface and related intercepts BC */

// !!! Normal vectors are renormalized with the L1-norm !!!
// pay attention to the "foreach_dimension()" declaration for the normal components
n_front.x[right] = -sin(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])));
n_front.y[right] = cos(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])));
alpha_front[right] = plane_alpha (f[ghost], (coord){
    -sin(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k]))), 
    cos(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])))
});

n_front.y[top] = -sin(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])));
n_front.x[top] = cos(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])));
alpha_front[top] = plane_alpha (f[ghost], (coord){
   -sin(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k]))), 
    cos(theta_0[k])/(fabs(sin(theta_0[k])) + fabs(cos(theta_0[k])))
});

  /* Volume Fractions BC */
f[right] = f0[1,0];
f[top] = f0[0,1];

  /* Contact Angles */
h.t[bottom] = contact_angle (pi/2.);
h.t[right] = contact_angle (pi/2. - theta_0[k]); 
h.t[top] = contact_angle (pi - theta_0[k]); 

  /* Velocity BC */

// Symmetry conditions on the axis:
u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0.);


// Capillary flow due to the curvature gradients in AXI:
u.t[right] = dirichlet ( y/(pow(x*x + y*y, 3./2.) * (tan(theta_0[k]))) );
u.n[right] = dirichlet ( x/(pow(x*x + y*y, 3./2.) * (tan(theta_0[k]))) );
u.n[top] = dirichlet ( y/(pow(x*x + y*y, 3./2.) * (tan(theta_0[k]))) );
u.t[top] = dirichlet ( x/(pow(x*x + y*y, 3./2.) * (tan(theta_0[k]))) );

/** 
### Generic Events
*/

int main() {

  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 1e-3, mu2 = 1e-5; 

  for (k = 0; k < 5; k++){    
    size(SIZE) ;
    init_grid(1 << LEVEL) ;
    origin(-X_OFFSET, 0) ;
    run() ;
  }
}


event init (t = 0){
  // Initial free-surface definition:
  #if TREE
    do{
    fraction(f0, - (y - x*tan(theta_0[k])) ) ; 
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, lvl_max[k]).nf); 

    // Initialization of the ghot cells'auxiliary field: 
    foreach_boundary(right)
      f0[1,0] = f_BC_right(f0[], f0[0,1], f0[0,-1]);
    foreach_boundary(top)
      f0[0,1] = f_BC_top(f0[], f0[1,0], f0[-1,0]);

    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
  #else
    fraction(f0, - (y - x*tan(theta_0[k])) ) ; 

    foreach_boundary(right)
      f0[1,0] = f_BC_right (f0[], f0[0,1], f0[0,-1]);
    foreach_boundary(top)
      f0[0,1] = f_BC_top(f0[], f0[1,0], f0[-1,0]);
  #endif
  // ---------------------------
  foreach(){
    f[] = f0[]; // do not forget to define the main volume fraction field
  }

  // computation of curvature for visualization:
  boundary({f});
  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
  curvature (f, curv_viz, 1, add = false);
}


event vorti (i++){
  // cf. bug report: https://groups.google.com/g/basilisk-fr/c/ok-OhtzO1Pk
  foreach()
    omega[] = y*(u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
  boundary({omega});
  curvature (f, curv_viz, 1, add = false);
}


event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, lvl_max[k], 6);
    unrefine(x < -7./24.*SIZE); // for backflow conditions
    unrefine( ((f[] == 1.) || (f[] == 0)) && (y < 10.*Delta) ); // for spurious vorticity on the axis
  #endif 

  foreach_boundary(right)
    f0[1,0] = f_BC_right (f0[], f0[0,1], f0[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top(f0[], f0[1,0], f0[-1,0]);

  #if TREE
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
  #endif

  boundary({f0});
  boundary({f});
  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
}


event end (t = t_end[k]){}


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


/** 
### Outputs
*/

event profiles (t = t_end[k]){                   
  char filename[200] ;
  sprintf(filename, "selfsim_sierou_shape_th%g_t%g_N%g.dat", 
          theta_0[k]*180/pi, t, lvl_max[k]) ;
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}






/**
## References
~~~bib
@article{keller1983,
 ISSN = {00361399},
 URL = {http://www.jstor.org/stable/2101434},
 author = {Joseph B. Keller and Michael J. Miksis},
 journal = {SIAM Journal on Applied Mathematics},
 number = {2},
 pages = {268--277},
 publisher = {Society for Industrial and Applied Mathematics},
 title = {Surface Tension Driven Flows},
 urldate = {2025-02-24},
 volume = {43},
 year = {1983}
}

@article{sierou2004,
  author = {Sierou, A. and Lister, J. R.},
  title = {Self-similar recoil of inviscid drops},
  journal = {Physics of Fluids},
  volume = {16},
  number = {5},
  pages = {1379-1394},
  year = {2004},
  month = {05},
  issn = {1070-6631},
  doi = {10.1063/1.1689031},
  url = {https://doi.org/10.1063/1.1689031},
  eprint = {https://pubs.aip.org/aip/pof/article-pdf/16/5/1379/19153172/1379\_1\_online.pdf},
}
~~~
*/