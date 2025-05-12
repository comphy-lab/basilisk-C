/** 
# 2D surface tension driven recoil of a liquid wedge: *different initial angles*

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

All the modifications to the native libraries used here 
are well explained in [this test case](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/).
Do not refrain to take a look to better understand the "tricks" 
used in this code for managing non-closed, moving interfaces.</div>


## Introduction

The numerical results of [Keller \& Miksis, (1983)](#keller1983) 
for a 2D liquid wedge recoil driven by surface tension are 
reproduced for all the initial wedge angles 
tested by the authors: $\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$°. 

Please, refer to [this page](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c) 
for a general presentation of the physical problem, and its related 
numerical convergence study.


## Reproduction of the *Keller \& Miksis* results

~~~pythonplot Basilisk *VS* Keller \& Miksis
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

# Use LaTeX fonts in the plot:
plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=20) # size for the title
plt.rc('axes', labelsize=17) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 


t = 1e-2
t_scaled = np.power(t,-2./3.)

#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
keller_27p5deg = np.loadtxt ("../keller_data/raw_keller_fig2_27.5.dat", 
                             unpack=False)
keller_32p5deg = np.loadtxt ("../keller_data/raw_keller_fig2_32.5.dat", 
                             unpack=False)
keller_45deg = np.loadtxt ("../keller_data/raw_keller_fig2_45.dat", 
                             unpack=False)
keller_65deg = np.loadtxt ("../keller_data/raw_keller_fig2_65.dat", 
                             unpack=False)
keller_80deg = np.loadtxt ("../keller_data/raw_keller_fig2_80.dat", 
                             unpack=False)

path_shape_27p5deg = "keller_shape_th27.5_t0.01_s1.dat"
path_shape_32p5deg = "keller_shape_th32.5_t0.01_s1.dat"
path_shape_45deg = "keller_shape_th45_t0.01_s1.dat"
path_shape_65deg = "keller_shape_th65_t0.01_s1.dat"
path_shape_80deg = "keller_shape_th80_t0.01_s1.dat"

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
# adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, maxlvl=11, minlvl=6);
# t_end = 1e-2
for segment in segments_27p5deg:
    ax.plot(segment[:, 0]*t_scaled, segment[:, 1]*t_scaled, 
            lw=1.2, color = mapcolors[0]) 
for segment in segments_32p5deg:
    ax.plot(segment[:, 0]*t_scaled, segment[:, 1]*t_scaled, 
            lw=1.2, color = mapcolors[2]) 
for segment in segments_45deg:
    ax.plot(segment[:, 0]*t_scaled, segment[:, 1]*t_scaled, 
            lw=1.2, color = mapcolors[4]) 
for segment in segments_65deg:
    ax.plot(segment[:, 0]*t_scaled, segment[:, 1]*t_scaled, 
            lw=1.2, color = mapcolors[6]) 
for segment in segments_80deg:
    ax.plot(segment[:, 0]*t_scaled, segment[:, 1]*t_scaled, 
            lw=1.2, color = mapcolors[8]) 

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


plt.savefig('keller_fig2_repro_N11_t0.01.svg') 
~~~
*/


/**
## Code 

### General Parameters

From the [convergence study](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c), 
its appears that, for a box of size unity in the 
non-dimensionalized physical space, $\tilde{t}_{end} = 10^{-2}$ 
corresponds to the time before the reflection of capillary 
waves on the boundaries. In that case, a maximum level of 
$N_{max} = 11$ is sufficient to capture the smallest structures 
with enough precision for having a sharp interface.
*/


#define LEVEL 7 
#define MAXLEVEL 11
#define SIZE 1e0
#define T_END 1e-2
#define X_OFFSET (1./3.)*SIZE

#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase_keller_fig2.h"
#include "tension.h"

// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 

scalar f0[]; // auxiliary volume fraction field used for BC 
vector n_front[]; // normal vector of the interface
scalar alpha_front[]; // intercept linked to n_front[] during VoF-reconstruction 

scalar u_L2[]; // for movie


/** 
 As said earlier, the following initial angles are tested:
 $\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$°. 
 */

int k; // counter
double theta_0[] = {27.5*pi/180., 32.5*pi/180., 
  45.*pi/180., 65.*pi/180, 
  80.*pi/180};
  
  
/**
### Boundary Conditions
*/

#include "f_BC_keller_fig2_all.h"

/**
We recall the following boundary conditions:

<img src="img_keller/keller_init_scheme.png" alt="Initial configuration with Boundary Conditions" height="350" class="center"/>

Also, as soon as the simulation begins, the contact angle on the 
symmetry axis becomes $90$°, and therefore does not remain at 
the initial $\theta_0$ value, but the far-field angle 
$\beta_0$ does not change during the simulation.

> **N.B.:** 

>  + all the "right" boundary conditions are designed for 
    $\theta_0 = 27.5$°, $32.5$°, $45$°;

>  + all the "top" boundary conditions are designed for 
    $\theta_0 = 65$°, $80$°;
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

// No flow on the right and top borders:
u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

// Outflow on the left gas-border (to set a pressure reference):
p[left] = dirichlet(0.);
pf[left] = dirichlet(0.);
u.n[left] = neumann(0.);


/**
### Generic Events
*/

int main() {
  
  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 0, mu2 = 0; // disable viscosity
  
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
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, MAXLEVEL).nf); 

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
  foreach()
    u_L2[] = sqrt(sq(u.x[]) + sq(u.y[])); 
}


event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, MAXLEVEL, LEVEL-1);
    unrefine (f[] == 0);
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


event end (t = T_END){
  dt = 1e-10;
}


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


/**
### Outputs

We want to plot the interface profile for different times of the 
surface tension driven recoil of the wedge:
*/


event profiles (t = T_END){               
  char filename[200] ;
  sprintf(filename, "keller_shape_th%g_t%g_s1.dat", theta_0[k]*180/pi, t) ;
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}


event vel_p_maps (t = T_END) 
{
  scalar x_interf[];
  scalar y_interf[];
  char fileup[200];
  char fileinterf[200];
  sprintf (fileup, "keller_u_p_th%g_t%g_s1.dat", theta_0[k]*180/pi, t);
  sprintf (fileinterf, "keller_interf_th%g_t%g_s1.dat", theta_0[k]*180/pi, t);
  FILE * fup = fopen (fileup, "w");
  FILE * finterf = fopen (fileinterf, "w");
  position (f, x_interf, {1,0}, add=false);
  position (f, y_interf, {0,1}, add=false);
  foreach (){
    fprintf(fup, "%g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[], omega[]);
    if ( interfacial (point, f) ) {
      if (f[] > 0.1 && f[] < 0.9){
        fprintf(finterf, "%g %g\n", x_interf[], y_interf[]);
      }
    }
  } 
  fclose(fup);
  fclose(finterf);
}





/**
## Additional Remarks 

+ Slight discrepancies between our results and those of *Keller \& Miksis* exist, 
that can be attributed to the differences of numerical methods used;

+ The authors used for instance the *trapezoidal rule*, which is a quite 
inaccurate method (order 1) compared to the *Simpson's rule* (order 3) or 
*Gauss quadrature*;

+ The authors recognized that their free surfaces were too flat far from 
the apex, due to the coarse mesh used, and can explain the difference observed 
in the number of capillary waves with `Basilisk` simulations;

+ *Keller \& Miksis* used a monophasic formulation, whereas a ratio of 
$10^3$ between liquid and gaz is employed here;

+ For $\theta_0 = 32.5$°, an important deviation is observed at the apex 
between our results and *Keller \& Miksis*, 
but we are positive that it might be a lack of convergence from the reference 
paper, as the gap between $\theta_0 = 27.5$° and $\theta_0 = 32.5$° is 
quite peculiar concerning the authors'curves. Also because all the other 
angles are well reproduced.
*/

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