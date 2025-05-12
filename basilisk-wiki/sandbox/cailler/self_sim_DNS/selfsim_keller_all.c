/** 
# Self-Similar DNS for the *Keller \& Miksis (1983)* problem: *different angles*

Reproduction in the self-similar domain of the [Keller \& Miksis, (1983)](#keller1983) 
wedge recoil driven by surface tension, for all the following initial wedge angles:
$\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$°. 

* The *Keller \& Miksis* problem is exposed [here](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c);
* The self-similar solver is explained [there](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README);
* This page is simply an extension of the [convergence study](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_keller_conv.c) 
presented for the special case of $\theta_0 = 45$°, and where additional details 
are given upon the parameters used.


~~~pythonplot Self-similar solver *VS* Keller \& Miksis
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
keller_27p5deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_27.5.dat", 
                             unpack=False)
keller_32p5deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_32.5.dat", 
                             unpack=False)
keller_45deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_45.dat", 
                             unpack=False)
keller_65deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_65.dat", 
                             unpack=False)
keller_80deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_80.dat", 
                             unpack=False)

path_shape_27p5deg = "selfsim_keller_shape_th27.5_t5_N9.dat"
path_shape_32p5deg = "selfsim_keller_shape_th32.5_t5_N9.dat"
path_shape_45deg = "selfsim_keller_shape_th45_t10_N8.dat"
path_shape_65deg = "selfsim_keller_shape_th65_t10_N8.dat"
path_shape_80deg = "selfsim_keller_shape_th80_t10_N8.dat"

facets_27p5deg = np.loadtxt (path_shape_27p5deg, unpack=False)
facets_32p5deg = np.loadtxt (path_shape_32p5deg, unpack=False)
facets_45deg = np.loadtxt (path_shape_45deg, unpack=False)
facets_65deg = np.loadtxt (path_shape_65deg, unpack=False)
facets_80deg = np.loadtxt (path_shape_80deg, unpack=False)

# Get segments of the output_facets() function of Basilisk C:
N_seg_27p5deg = int (0.5*facets_27p5deg.shape[0])
N_seg_32p5deg = int (0.5*facets_32p5deg.shape[0])
N_seg_45deg = int (0.5*facets_45deg.shape[0])
N_seg_65deg = int (0.5*facets_65deg.shape[0])
N_seg_80deg = int (0.5*facets_80deg.shape[0])

segments_27p5deg = np.split (facets_27p5deg, indices_or_sections=N_seg_27p5deg)
segments_32p5deg = np.split (facets_32p5deg, indices_or_sections=N_seg_32p5deg)
segments_45deg = np.split (facets_45deg, indices_or_sections=N_seg_45deg)
segments_65deg = np.split (facets_65deg, indices_or_sections=N_seg_65deg)
segments_80deg = np.split (facets_80deg, indices_or_sections=N_seg_80deg)


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 10
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,4)
ax.set_ylim(0,4)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')

# Basilisk results with:
# ----------------------
# adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, maxlvl=8, minlvl=6);
# t_end[k] = 10
for segment in segments_27p5deg:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[0]) 
for segment in segments_32p5deg:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[2]) 
for segment in segments_45deg:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[4]) 
for segment in segments_65deg:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[6]) 
for segment in segments_80deg:
    ax.plot(segment[:, 0], segment[:, 1], lw=1.2, color = mapcolors[8]) 

# Keller & Miksis (1983) results:
# -------------------------------
ax.plot(keller_27p5deg[:,0], keller_27p5deg[:,1], '--', lw=0.8, color="darkred")
ax.plot(keller_32p5deg[:,0], keller_32p5deg[:,1], '--', lw=0.8, color="darkred")
ax.plot(keller_45deg[:,0], keller_45deg[:,1], '--', lw=0.8, color="darkred")
ax.plot(keller_65deg[:,0], keller_65deg[:,1], '--', lw=0.8, color="darkred")
ax.plot(keller_80deg[:,0], keller_80deg[:,1], '--', lw=0.8, color="darkred", 
        label=r'Keller & Miksis')

ax.legend(frameon=False, loc=4)
plt.text(3.6, 1.6, r'$\theta_0 = 27.5$°', fontsize=14,
         ha='center', va='center', rotation=27.5, color=mapcolors[0])
plt.text(3.5, 2.43, r'$\theta_0 = 32.5$°', fontsize=14,
         ha='center', va='center', rotation=32.5, color=mapcolors[2])
plt.text(3, 3.2, r'$\theta_0 = 45$°', fontsize=14,
         ha='center', va='center', rotation=45, color=mapcolors[4])
plt.text(1.46, 3.45, r'$\theta_0 = 65$°', fontsize=14,
         ha='center', va='center', rotation=65, color=mapcolors[6])
plt.text(.47, 3.5, r'$\theta_0 = 80$°', fontsize=14,
         ha='center', va='center', rotation=80, color=mapcolors[8])


plt.savefig('selfsim_solver_keller_all_deg_N8.svg') 
~~~

As shown on the above figure: 

* all the curves obtained with the self-similar solver in 2D 
reproduce correctly the numerical results of *Keller \& Miksis*: **this validates 
our self-similar solver in 2D**;
* smaller angles need better refinement to be satisfying 
[not possible on the `Basilisk` server due to running time constraints]; 
* no visible effect of the viscosity concerning the shapes of the interfaces, 
compared to the potential results of the reference paper. 

*/


/** 
## Code

### General Parameters 

From the convergence study it appears that $\tau_{end} = 10$ with a maximum 
level of refinement $N_{max} = 8$ should be enough to converge towards the 
results of *Keller \& Miksis*, with respect to the server running time constraints. 

**However**, for smaller angles, the high curvature variations need a better 
refinement. Considering the time constraints with the online server for the 
computations, we suggest to go up to $N_{max} = 9$ for $\theta_0 = 27.5$° and 
$32.5$°, with a running time divided by 2 ($\tau_{end} = 5$ in that case). 
The results are slightly improved, though it is not satisfying compared to 
the results one could obtain by using $N_{max} = 10$ to match with the [simulations 
done in the physical domain](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_all.c#reproduction-of-the-keller-miksis-results).
*/


#define LEVEL 7 
#define SIZE 12e0
#define X_OFFSET (1./3.)*SIZE


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

/** We test for the following initial wedge angles:
$\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$° */
int k; // counter
double theta_0[] = {27.5*pi/180., 32.5*pi/180., 45.*pi/180., 65.*pi/180, 80.*pi/180};
double lvl_max[] = {9, 9, 8, 8, 8};
double t_end[] = {5., 5., 10., 10., 10.};

/** 
### Boundary Conditions
*/

#include "selfsim_f_BC_keller_all.h"

/** 
For a clear depiction of the BCs, see 
[this representation](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_all.c#boundary-conditions).
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

// No flow on the opposite side of the apex:
u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);


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
  vorticity (u, omega);
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
  sprintf(filename, "selfsim_keller_shape_th%g_t%g_N%g.dat", 
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
~~~
*/