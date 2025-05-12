/** 
# Self-similar recoiling cones in the presence of a negative *dipolar flow* 

This page is mirroring [this one](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_dipolar_flow.c). 
As such, please **refer to the link provided** for explanations upon motivations, 
parameters, source code tips...  

The current page deals with the case of a **negative** far-field dipolar flow:

~~~pythonplot Parametric study for $\theta_0 = 120$° and negative dipolar flow
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

t = 0.02
t_str = '0.02'
th = 120
N = 7
mu0_list = range(0,-11,-1)

t_scaled = np.power(t,-2./3.)


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 25
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_aspect('equal')

# If we want a vertical axis of symmetry:
# ---------------------------------------

ax.spines['left'].set_color('red')    # change color
ax.spines['left'].set_linestyle((0, (3, 5, 1, 5))) # tuple (0, (3, 5, 1, 5)) for `dashdotted`` line
ax.set_xlim(0,8)
ax.set_ylim(0,10)
ax.set_xlabel(r'$\eta$')
ax.set_ylabel(r'$\xi$')


for k, mu0 in enumerate(mu0_list):
  if k == 0:
    path_shape = "../sierou_data/self_focusing_sing_shape_th{}_t{}_N9_m{}.dat".format(th,t_str,mu0)
  else:
    path_shape = "self_focusing_sing_shape_th{}_t{}_N{}_m{}.dat".format(th,t_str,N,mu0)
  facets = np.loadtxt (path_shape, unpack=False)
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments):
    if ((k > 0) and (k < len(mu0_list) - 1)):
      ax.plot(segment[:, 1]*t_scaled, -segment[:, 0]*t_scaled, 
              lw=1.2, color = mapcolors[2*k]) 
    else:
      ax.plot(segment[:, 1]*t_scaled, -segment[:, 0]*t_scaled, 
              lw=1.2, color = mapcolors[2*k])      


# Labels for μ0 = 0 and μ0 = -10:
plt.text(6.5, 3.3, r'$\widetilde{\mu}_0 = 0$', rotation=30, color=mapcolors[0], fontsize=12)
plt.text(5.8, 5.7, r'$\widetilde{\mu}_0 = -10$', rotation=25, color=mapcolors[2*len(mu0_list)], fontsize=12)

# Plot Sierou curves for θ_0 = 120°:
for k in range(1,6):
  sierou = np.loadtxt("../sierou_data/fig16b_sierou_curve{}.dat".format(k))
  if k == 1:
    ax.plot(sierou[:,0], sierou[:,1], 
            '--', lw = 0.8, color='darkred',
            label=r'Sierou & Lister')
  else:
    ax.plot(sierou[:,0], sierou[:,1], '--', lw = 0.8, color='darkred')

ax.legend(frameon=False, loc=1)

plt.savefig('fig16b_repro_th120_t0.02_N7_vert_axis_mu0_neg.svg') 
~~~

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Note to the reader**

The above plot was generated for a maximum level of refinement of 
$N_{max} = 7 \Rightarrow \widetilde{\Delta} = 7.8 \times 10^{-3}$, hence, a 
***very low resolution*** to run on the `Basilisk` server. For your own 
purpose, it is suggested to go as high as $N_{max} = 9$ for better accuracy 
and convergence of the apex shapes (see the 
[*Discussion* section](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_neg_dipolar_flow.c.#discussion)).
</div>


## Code

### General Parameters
*/

#define LEVEL 6

#define MAXLEVEL 7
#define SIZE 1 
#define T_END 0.02 

#define THETA_0 120.*pi/180.
#define BETA_0 pi - THETA_0
#define MU_0  -1

#define Y_OFFSET 0. 
#define X_OFFSET 0.8*SIZE

#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"

#include "f_BC_dipolar_flow.h" // BCs in ghost cells


// Vectors, fields...
scalar omega[];
scalar f0[];
vector h[];
scalar curv_viz[]; // to store curvature 
vector n_front[];
scalar alpha_front[];
scalar visc_dissip[];

int k; // counter
int mu_0[] = {0, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1}; // dipole potential strength

#include "u_BC_dipolar_flow.h" // far-field velocities for MU_0 ≠ 0
#include "tag.h" // to remove erroneous micro-bubbles


/** 
### Boundary Conditions
*/

// Normal vector in the ghost cells, to break the default symmetry condition:
// Pay attention to the "foreach_dimension()" declaration for the normal components
// (Initially located in the `vof.h` file, moved here for convenience.)
//                  !! L1 norm is used for normals !!

n_front.x[left] = - sin(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)));
n_front.y[left] = - cos(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)));
alpha_front[left] = plane_alpha (f[ghost], (coord){
      - sin(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0))), 
      - cos(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)))
      });

n_front.y[top] = - sin(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)));
n_front.x[top] = - cos(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)));
alpha_front[top] = plane_alpha (f[ghost], (coord){
      - sin(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0))), 
      - cos(BETA_0)/(fabs(-sin(BETA_0)) + fabs(-cos(BETA_0)))
      });


// Volume Fraction BC:
f[left] = f0[-1,0];
f[top] = f0[0,1];

// Contact Angles:
h.t[left] = contact_angle (pi/2. + BETA_0) ; 
h.t[top] = contact_angle (BETA_0) ;

// On bottom boundary (AXIS of SYMMETRY):
u.n[bottom] = dirichlet (0.);
u.t[bottom] = neumann (0.);


#if MU_0 != 0
// without far-field flow --> dirichlet(0) on u!!

  // On left boundary, outflow conditions are imposed for the external fluid.
  // WARNING! The interface can jump from top to left border!
  u.n[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uxA(x,y) );
  u.t[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uyA(x,y) );

  // On top boundary, exterior and interior fluids coexist for most of the
  // simulation:
  u.n[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uyA(x,y) );
  u.t[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uxA(x,y) );

  // On right boundary (will always be liquid phase):
  u.n[right] = dirichlet ( uxA(x,y) );
  u.t[right] = dirichlet ( uyA(x,y) );


  // Pressure BC:
  p[left] = (f[] == 0. ? dirichlet (0.) 
    : neumann ( - ( fm.n[0]/alpha.n[0] )
    *( a.n[0] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) );

  pf[left] = (f[] == 0. ? dirichlet (0.) 
    : neumann ( - ( fm.n[0]/alpha.n[0] )
    *( a.n[0] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) );


  p[top] = ( f[] == 1. ?
    neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uyA(x,y)*duyA_dy(x,y) - uxA(x,y)*duyA_dx(x,y)) )
    : dirichlet (0.) );

  pf[top] = ( f[] == 1. ?
    neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uyA(x,y)*duyA_dy(x,y) - uxA(x,y)*duyA_dx(x,y)) )
    : dirichlet (0.) );


  p[right] = neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) ;

  pf[right] = neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) ;


#else // MU_0 = 0
  // On left boundary (will always be gas phase):
  u.n[left] = dirichlet (0.);
  uf.n[left] = 0.;
  u.t[left] = dirichlet (0.);
  uf.t[left] = 0.;

  // On top boundary, exterior and interior fluids coexist:
  u.n[top] = dirichlet (0.);
  uf.n[top] = 0.;
  u.t[top] = dirichlet (0.);
  uf.n[top] = 0.;

  // On right boundary (will always be liquid phase):
  u.n[right] = dirichlet (0.);
  u.t[right] = dirichlet (0.);
#endif




/** 
### Generic Events
*/

int main() {

  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 0, mu2 = 0;

  for (k = 1; k < 11; k++){    
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
    fraction(f0, (y - (- (x) * tan(BETA_0) + Y_OFFSET) ) ) ; 
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, MAXLEVEL, LEVEL).nf); 
  #else //!TREE
    fraction(f0, (y - (- (x) * tan(BETA_0) + Y_OFFSET) ) ) ;
  #endif

  foreach_boundary(left)
    f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (f0[], f0[1,0], f0[-1,0]);

  #if TREE
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
  #endif 

  boundary({f0});

  foreach(){
    f[] = f0[];
    #if MU_0 != 0
      u.x[] = f[]*uxA(x,y) + (1. - f[])*uxO(x,y); 
      u.y[] = f[]*uyA(x,y) + (1. - f[])*uyO(x,y);
    #endif
  }

  #if MU_0 != 0
    boundary((scalar *){u});
  #endif

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
  // vorticity in AXI:
  foreach()
    omega[] = y*(u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);;
  boundary({omega});

  // curvature
  curvature (f, curv_viz, 1, add = false);
}


event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1e-3,1e-2,1e-2}, MAXLEVEL, LEVEL);
    // unrefine(f[] == 0);
  #endif 

  #if MU_0 != 0
    foreach_boundary(left)
      f0[-1,0] = f_BC_left (f[], f[0,1], f[0,-1]);
    foreach_boundary(top)
      f0[0,1] = f_BC_top (f[], f[1,0], f[-1,0]);
  #else // MU_0 = 0
    foreach_boundary(left)
      f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);
    foreach_boundary(top)
      f0[0,1] = f_BC_top (f0[], f0[1,0], f0[-1,0]);
  #endif

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

  // Try to remove all parasitic micro bubbles forming 
  remove_droplets (f, minsize = 1, threshold = 1.e-3);
  remove_droplets (f, minsize = 1, threshold = 1.e-3, bubbles=true);

  foreach (){
    if (f[] < 5e-4){
      f[] = 0.;
    }
    if ( (1. - f[]) < 5e-4 ){
      f[] = 1.;
    }
  }
  boundary({f});
}



event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // vider le buffer
}

event end (t = T_END){}


/** 
### Outputs
*/

event profiles (t = T_END){                   
  char filename[200] ;
  sprintf(filename, "self_focusing_sing_shape_th%g_t%g_N%d_m%d.dat", 
          THETA_0*180/pi, t, MAXLEVEL, mu_0[k]) ;
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}



/** 
## Discussion

A full study -- not suited for the `Basilisk` server -- has been conducted 
locally for two different angles ($\theta_0 = 120$° and $145$°) and a 
broad range of negative dipolar flow strengths 
($\widetilde{\mu}_0 \in [-10 \, ; 0]$), 
with optimized parameters for a better refinement, without sacrifying 
too much computational speed: 
$\tilde{t}_{max} = 2 \times 10^{-2}$, 
$N_{max} = 9 \,\,  \longleftrightarrow \,\, 
\widetilde{\Delta}_{opt} = 1.95 \times 10^{-3}$. 

The results are gathered on the below figures:

![`Basilisk` results for $\theta_0 = 120$° (*left*) and $\theta_0 = 145$° (*right*, 
not done by *S&L*)](img_sierou/sierou_mu-neg_v2.png){width="50%"}

+ Excellent agreement between the `Basilisk` results and 
[Sierou \& Lister](#sierou2004) for $\widetilde{\mu}_0 \in [-4 \, ; 0]$. 

+ **Extension** of the reference paper results for *higher* $\widetilde{\mu}_0$ 
*values* and other $\theta_0$ angles: we show for the first time capillary-inertial 
**self-similar jet structures** for high $\widetilde{\mu}_0$ values, therefore 
that a strong far-field dipolar flow can be responsible for a *curvature reversal!*

+ Under the influence of a *negative* dipolar flow field, the overall movement 
(*capillary flow*) pushing forward the obtuse liquid cone is proportionally 
**strengthened** by the intensity of the dipolar flow. 

+ This behaviour is characterized in the self-similar space by a shift of the 
interface profiles according to ascending $\xi$. 

+ Gaps between nadirs increase with the absolute value of the dipolar flow 
intensity. This is the opposite for *positive* dipolar flows.  




## References
~~~bib
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

