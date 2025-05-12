/** 
# 2D surface tension driven recoil of a liquid wedge: *a convergence study*

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

All the modifications to the native libraries used here 
are well explained in [this test case](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/).
Do not refrain to take a look to better understand the "tricks" 
used in this code for managing non-closed, moving interfaces.</div>


## Introduction

We try to reproduce the numerical results of 
[Keller \& Miksis, (1983)](#keller1983) for a 2D liquid
wedge recoil driven by surface tension. 

This problem is **scale invariant**, as one can 
infer from the beautiful sketch below by 
[Peregrine *et al.*, (1983)](#peregrine1990):

![Wedge recoil drawn by [Peregrine *et al.*, (1983)](#peregrine1990)](img_keller/recoil_peregrine.png)

and by the resulting movie from the current `Basilisk C` code page:

![Wedge recoil with `Basilisk`](keller_fig2_conv/k+m_phys_space_N11_test.mp4)(width="600" height="600")


Indeed, *Keller \& Miksis* have shown that the physical variables 
at stake are the space and time coordinates $x, y$ and $t$, 
the density $\rho$ and the surface tension $\sigma$. 
Hence, for an inviscid, incompressible and irrotational flow, 
there exists a *velocity potential* $\varphi$ that can be 
non-dimensionalized by the triplet $\left\{\sigma, \rho, t \right\}$ 
by dimensional analysis (*Vaschy-Buckingham theorem*):

$$
\overline{\varphi}(\xi, \eta) 
= 
\left( \dfrac{\sigma^2 t}{\rho^2} \right)^{-1/3} 
\, \varphi(x,y,t)
$$

with:

$$
\xi = \left( \dfrac{\rho}{\sigma t^2} \right)^{1/3} \, x 
\quad , \quad
\eta = \left( \dfrac{\rho}{\sigma t^2} \right)^{1/3} \, y 
$$

where $\xi$ and $\eta$ are the **self-similar variables**.

**In particular:** lengths scale like $t^{2/3}$ and velocities like $t^{-1/3}$. 
When $t \to 0$, the velocity diverges: there is a *singularity*.

Here, we demonstrate this *scale invariance* by only 
showing the validation of the convergence study done 
for a wedge with a $45$° half-angle. 
For the collection of all the initial wedge angles 
tested by *Keller \& Miksis* 
($\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$°), 
please refer to [this page](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_all.c).


## Methodology

In order to exhibit the scale invariance, one can perform simulations in 
the *physical* domain, gets the interface profiles for different times, 
and rescales them by dividing them by $t^{2/3}$. 
If all rescaled profiles collapse onto a single one (*master curve*), 
then the problem is *steady* in this new rescaled frame called 
"the *self-similar space*": the problem is said to be *scale invariant*.

The simulation box of size $L_0$ has a resolution $\Delta$ related to the 
maximum level of discretization $N_{max}$: $\Delta = L_0/2^{N_{max}}$. 
The smallest structure resolved has a size of approximately 
$\ell_{min} \approx 10 \Delta$.

### Time Limits

Lengths scale as $t^{2/3}$, therefore wavelengths of surface waves generated 
by the capillary recoil evolve also like $t^{2/3}$. Now, the smallest 
structure numerically resolved wirtes:

$$
\ell_\sigma :=
\left(\dfrac{\sigma}{\rho} \right)^{1/3}{t_{min}}^{2/3} = 10 \, \Delta 
\equiv \ell_{min}
\Rightarrow 
t_{min} = \sqrt{\dfrac{\rho}{\sigma}} 
\, \left(10 \dfrac{L_0}{2^{N_{max}}}\right)^{3/2}
$$

where $t_{min}$ is the time from which observed spatial structures are 
large enough to be properly simulated.

As for the *final time* of the simulation, it is obtained when capillary waves 
have reached the opposite border of the apex and reflect. The *dispersion 
relation* is given by [Lamb, 1916, p.252, §267](#lamb1916): 

$$
\omega^2 = \dfrac{\sigma}{\rho} \, k^3
$$

with $\omega$ the circular frequency, and $k$ the wavenumber. 
Hence, the group velocity $v_g = \mathrm{d}\omega/\mathrm{d}k$ is:

$$
v_g = \dfrac{3}{2} \sqrt{\sigma/\rho} \, k^{1/2} 
$$

The wavenumber scales like 
$k \sim {\ell_\sigma}^{-1} \sim (10 \, \Delta)^{-1}$, so the *maximum timestep* 
is obtained for:

$$
t_{max} 
\sim 
\dfrac{L_0}{v_g}
\sim \dfrac{2}{3} \sqrt{\dfrac{\rho}{\sigma}} \, L_0 \, 
\sqrt{10 \, \Delta}
$$



### Non-dimensional variables

We non-dimensionalize the problem to work with quantities of the 
same order of magnitude. A natural choice is then:

$$
u = \sqrt{\dfrac{\sigma}{L_0 \,\rho} } \, \tilde{u} 
\quad ; \quad 
t = \sqrt{\dfrac{L_0^3 \, \rho}{\sigma} } \, \tilde{t} 
\quad ; \quad 
p = \dfrac{\sigma}{L_0} \tilde{p}
\quad ; \quad 
\kappa = \dfrac{1}{L_0} \widetilde{\kappa} 
$$
$$
x = L_0 \, \tilde{x}
\quad ; \quad 
y = L_0 \, \tilde{y}
\quad ; \quad 
\Delta = L_0 \, \widetilde{\Delta}
$$

Now, the non-dimensionalized time limits read as:
$$
\tilde{t}_{min} = \left(\dfrac{10}{2^{N_{max}}}\right)^{3/2} 
\quad ; \quad 
\tilde{t}_{max} = \dfrac{2}{3}\left(\dfrac{10}{2^{N_{max}}}\right)^{1/2} 
$$

Finally, the non-dimensionalization is done with the liquid density 
and surface tension. A ratio $\rho_{liq}/\rho_{gas} = 10^3$ is taken 
for modelling the water/air interface. 


### Simulations done to show self-similarity

An adaptive mesh is used to refine the wedge apex where most of the 
physics happens. The convergence study is done only for the initial 
angle $\theta_0 = 45$° in order to optimize the resolution $\widetilde{\delta}$ 
for the final time of the simulation $\tilde{t}_{max}$ assuring a good 
balance between a sharp interface and quick computations. 

Besides, very separated time scales are scanned to show the 
self-similarity behaviour. To minimize the time computation, the 
following strategy has been applied: a bunch of simulations is launched 
for various end times corresponding to the smallest structures 
correctly observed at a given maximal resolution, that is to say, for 
$\tilde{t}_{end} \equiv \tilde{t}_{min}$. 

The following set of parameters/simulations has been chosen for proving 
the scale invariance: 

  *  $N_{max} = 11 \Rightarrow \widetilde{\Delta} = 4.9 \times 10^{-4}$ 
      for $\tilde{t} \leqslant \tilde{t}_{max} \approx 5 \times 10^{-2}$;  
  *  $N_{max} = 13 \Rightarrow \widetilde{\Delta} = 1.2 \times 10^{-4}$ 
      for $\tilde{t} \leqslant 10^{-4}$;  
  *  $N_{max} = 16 \Rightarrow \widetilde{\Delta} = 1.5 \times 10^{-5}$ 
      for $\tilde{t} \leqslant 10^{-5}$;  
  *  $N_{max} = 18 \Rightarrow \widetilde{\Delta} = 3.8 \times 10^{-6}$ 
      for $\tilde{t} \leqslant 10^{-6}$. 





*/


/**
## Code 

### General Parameters
*/


#define LEVEL 7 
#define SIZE 1e0
#define THETA_0 45.*pi/180.
#define X_OFFSET (1./3.)*SIZE

#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase_keller_fig2.h"
#include "tension.h"
#include "f_BC_keller_fig2.h"

// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 

scalar f0[]; // auxiliary volume fraction field used for BC 
vector n_front[]; // normal vector of the interface
scalar alpha_front[]; // intercept linked to n_front[] during VoF-reconstruction 

scalar u_L2[]; // for movie

int k; // counter
int lvl_max[] = {18, 16, 13, 13, 11};
double t_end[] = {1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2};

/**
### Boundary Conditions

The following boundary conditions are applied, 
keeping in mind that the velocity is considered to be zero 
on the opposite side of the apex:

<img src="img_keller/keller_init_scheme.png" alt="Initial configuration with Boundary Conditions" height="350" class="center"/>

As soon as the simulation begins, the contact angle on the 
symmetry axis becomes $90$°, and therefore does not remain at 
the initial $\theta_0$ value. Nonetheless, the far-field angle 
$\beta_0$ does not change.

> **N.B.:** 

>  + all the "right" boundary conditions are designed for 
    $\theta_0 = 27.5$°, $32.5$°, $45$°;

>  + all the "top" boundary conditions are designed for 
    $\theta_0 = 65$°, $80$°;
*/


  /* Normal vectors of the interface and related intercepts BC */

// !!! Normal vectors are renormalized with the L1-norm !!!
// pay attention to the "foreach_dimension()" declaration for the normal components
n_front.x[right] = -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)));
n_front.y[right] = cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)));
alpha_front[right] = plane_alpha (f[ghost], (coord){
    -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
    cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
});

n_front.y[top] = -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)));
n_front.x[top] = cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)));
alpha_front[top] = plane_alpha (f[ghost], (coord){
   -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
    cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
});

  /* Volume Fractions BC */
f[right] = f0[1,0];
f[top] = f0[0,1];

  /* Contact Angles */
h.t[bottom] = contact_angle (pi/2.);
h.t[right] = contact_angle (pi/2. - THETA_0); 
h.t[top] = contact_angle (pi - THETA_0); 

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
    fraction(f0, - (y - x*tan(THETA_0)) ) ; 
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, lvl_max[k]).nf); 

    // Initialization of the ghot cells'auxiliary field: 
    foreach_boundary(right)
      f0[1,0] = f_BC_right(f0[], f0[0,1], f0[0,-1]);
    foreach_boundary(top)
      f0[0,1] = f_BC_top(f0[], f0[1,0], f0[-1,0]);

    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
  #else
    fraction(f0, - (y - x*tan(THETA_0)) ) ; 

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
    adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, lvl_max[k], LEVEL-1);
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


event end (t = t_end[k]){
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


event profiles (t = t_end[k]){               
  char filename[200] ;
  sprintf(filename, "keller_shape_th%g_t%g_s1.dat", THETA_0*180/pi, t) ;
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}


/**
Once that we have got the interface profiles, we divide all the lengths 
by $t^{2/3}$ to exhibit the self-invariance. In the figure below, all 
the curves for different times collapse onto a master curve: the problem has 
become *steady* in the *self-similar space*, compared to the non-steady 
situation in the *physical space* represented in the inset. 

Notice also the perfect match between our results with `Basilisk` and 
those of *Keller \& Miksis*.

~~~pythonplot Numerical convergence \& Scale Invariance Proof
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Use LaTeX fonts in the plot:
plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=20) # size for the title
plt.rc('axes', labelsize=17) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 


L_0 = 1.

times_str = ["1e-06", "0.001", "0.01"] 
times_float = np.array([1e-6, 0.001, 0.01])
times_scaled = np.power(times_float, -2./3.)

#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
keller_45deg = np.loadtxt ("../keller_data/raw_keller_fig2_45.dat", 
                             unpack=False)

#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 10
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_aspect('equal')
axins = inset_axes(ax, width=1.5, height=1.5, borderpad=1, loc=3, 
                   bbox_to_anchor=(.55, .0625, .3125, .3125),
                   bbox_transform=ax.transAxes)
axins.patch.set_color('#fffcf4')

ax.set_xlim(0,4)
ax.set_ylim(0,4)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')

# Labels/Lim physical space
axins.set_title(r'$\tilde{x}$', fontsize=18)
axins.set_ylabel(r'$\tilde{y}$')
axins.set_xlim([-0.001/L_0, .042/L_0])
axins.set_ylim([0, .042/L_0])

j = 0
for k, t in enumerate(times_str):
  facets = np.loadtxt ("keller_shape_th45_t{}_s1.dat".format(t))
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments):
    if l == 0:
      ax.plot(segment[:, 0]*times_scaled[k], segment[:, 1]*times_scaled[k], 
              lw=1.2, color = mapcolors[k+1+j], 
              label=r'$\tilde{t}$' + r'$\,\, = {}$'.format(t))
      axins.plot(segment[:, 0]/L_0, segment[:, 1]/L_0, lw=1.2, color = mapcolors[k+1+j])
    else:
      ax.plot(segment[:, 0]*times_scaled[k], segment[:, 1]*times_scaled[k], 
              lw=1.2, color = mapcolors[k+1+j]) 
      axins.plot(segment[:, 0]/L_0, segment[:, 1]/L_0, lw=1.2, color = mapcolors[k+1+j])
  j += 1

# Keller & Miksis results:
ax.plot(keller_45deg[:,0], keller_45deg[:,1], '--', lw=0.8, color="darkred",
        label=r'Keller & Miksis')

# Time arrow:
axins.annotate(r' ', xy=(0.026/L_0, 0.035/L_0), 
              xytext=(.001/L_0, .01/L_0),
              arrowprops=dict(arrowstyle="->", 
                            color = "grey"
                            # width=0.005, 
                            # head_width=6*0.005,
                            ),
              color = "grey",
              fontsize=12
              )
axins.text(.013/L_0, .028/L_0, r'$\tilde{t}$', fontsize=12,
         ha='center', va='center', color="grey")

ax.legend(frameon=False, loc=2)

plt.savefig('keller_fig2_selfim_conv.svg') 
~~~

**CONCLUSION:** to reconcile time execution and precision up to 
$\tilde{t}_{max}$, the maximum level of refinement 
$N_{max} = 11 \widetilde{\Delta}_{opt} = 4.9 \times 10^{-4}$ will be chosen 
for the [full study](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_all.c) 
of all wedge angles of *Keller \& Miksis*.

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

Unfortunately, due to [a specific bug](http://basilisk.fr/sandbox/bugs/successive_runs.c), 
there is a problem with the timestep `dt` computation between two consecutive 
runs, and the profiles for $\tilde{t} = 10^{-5}$ and $\tilde{t} = 10^{-4}$ 
are not resolved at all (so we do not show them here).

**However**, this problem does not exist by running independtly the simulations. 
For instance for $\tilde{t} = 10^{-5}$, we can replace the current code 
with the following lines:
```C
int lvl_max[] = {16};
double t_end[] = {1.e-5};

...

for (k = 0; k < 1; k++){ // main
  ...
}

Comment movie event or change the condition on counter `k`value to `k == 1`. 
```
and everything runs smoothly (*i.e.* the timestep is correctly initialized and 
computed until reaching $\tilde{t} = 10^{-5}$ when we do run a single 
simulation and not successive ones.)
</div>

*/


event vel_p_maps (t = t_end[k]) // (t = {1.e-5, 1.e-4, 1.e-3, 1.e-2})
{
  scalar x_interf[];
  scalar y_interf[];
  char fileup[200];
  char fileinterf[200];
  sprintf (fileup, "keller_u_p_th%g_t%g_s1.dat", THETA_0*180/pi, t);
  sprintf (fileinterf, "keller_interf_th%g_t%g_s1.dat", THETA_0*180/pi, t);
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
### Movie
*/

#include "view.h"
#include "cmaps.h"

event mov_p_u (t+=1e-4) {
  if (k > 3){
    char legend[100];
    sprintf(legend, "t = %0.2g", t);
    view (tx = -0.05, ty = -0.01, width = 800., height = 800., fov = 8.); // for t_end=1e-2
    box();
    squares ("p", spread=-1, linear=true, map=cool_warm);
    draw_vof ("f", lw = 1.5);
    mirror ({0.,1}) {
      draw_vof ("f", lw = 1.5);
      squares("u_L2", spread=-1, linear=true, map=viridis);
      vectors("u", scale=0.0001, lc = {1, 1, 1});
    }
    draw_string(legend, 1, size = 20., lw = 2.1);
    save ("k+m_phys_space_N11_test.mp4");
  }
}


/**
## Some Remarks 

Despite our efforts concerning the improvements of the curvature 
computations at contact angles, the simulation needs to be stopped 
as soon as capillary waves reach the opposite side of the wedge tip, 
since it still creates reflexive waves (they are not really going 
outside of the box in spite of the changes made in 
[`vof_keller_fig2.h`](http://basilisk.fr/sandbox/cailler/keller_miksis/vof_keller_fig2.h)).
Also, spurious capillary waves are still generated by the AMR at 
opposite contact angle.

Therefore, for the 2D case of *Keller \& Miksis*, all the custom 
modifications are not necessarily needed, and symmetry conditions 
("naive code") could be totally justified in this case.

**However:** for very cute contact angles at opposite wedge tip side ($< 20$°), 
our treatment of height functions limits significantly vorticity 
perturbations contrary to default boundary conditions. 
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

@article{peregrine1990,
  title={The bifurcation of liquid bridges},
  author={Peregrine, D Howell and Shoker, G and Symon, A},
  journal={Journal of Fluid Mechanics},
  volume={212},
  pages={25--39},
  year={1990},
  publisher={Cambridge University Press},
  URL = {https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/bifurcation-of-liquid-bridges/A2E5D4B09C6214E085A1AE49C589F5C3#},
}

@book{lamb1916,
  title={Hydrodynamics},
  author={Lamb, H},
  publisher={Cambridge University Press},
  URL = {https://archive.org/details/hydrodynamics02lambgoog},
}
~~~
*/