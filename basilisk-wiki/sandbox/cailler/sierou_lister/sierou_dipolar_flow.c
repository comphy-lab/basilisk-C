/** 
# Self-similar recoiling cones in the presence of a (positive) *dipolar flow* 

## Context of the study

When bubbles are bursting at the ocean free surface, they form during their 
collapse a *conical cavity* whose angle $\theta_0$ with the axis of symmetry is 
greater than $90$°, counting from the liquid phase:

![Collapsing cavity of a bursting bubble, 
[Poujol *et al.*, (2021)](#poujol2021)](img_sierou/sierou_poujol.png){width="40%"}

For instance, in the experiment of [Zeff *et al.*, (2000)](#zeff2000) the 
observed conical angle is $\theta_0^{Zeff} = 120$°, whereas in the numerical 
simulations of [Duchemin *et al.*, (2002)](#duchemin2002), an angle 
$\theta_0^{Duchemin} = 145$° has been computed:

![*Left:* [Zeff *et al.*, (2000)](#zeff2000). 
*Right* [Duchemin *et al.*, (2002)](#duchemin2002)](img_sierou/sierou_zeff_duchemin.png){width="40%"}

They exhibit a capillary-inertial **self-similar behaviour** following the 
scaling of [Keller \& Miksis, (1983)](#keller1983), where lengths scale as 
$r_\sigma = (\sigma t^2 / \rho)^{1/3}$; see also [this `Basilisk` page](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c).

According to [Sierou \& Lister, (2004)](#sierou2004), and under the assumption 
of *potential flows*, collapsing conical cavities correspond to the ***time 
reversal*** of the capillary-inertial recoil of conical interfaces having a 
$\theta_0$ angle greater than $90$°. 
Therefore, the authors have used this property to model *self-focusing 
singularities*: instead of trying to start from the burst of the bubble like 
in the [simulation of A. Berny](http://basilisk.fr/sandbox/aberny/bubble/bursting2D.c), 
which requires a complex initialization, they simply start from the supposed 
final state of a *cone*, and study its behaviour under the spontaneous capillary 
flow developed in 3D--AXI, but also in conjunction with a far-field *dipolar 
flow* explained in [this documentation](http://basilisk.fr/sandbox/cailler/sierou_lister/README).

**The goal of this `Basilisk` code** is to reproduce the results of the 
*Fig. 16b* of [Sierou \& Lister, (2004)](#sierou2004), when the far-field 
dipolar flow is *positive* and $\theta_0 = 120$°.
For *negative* far-field dipolar flow intensities, 
please see [this source code](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_neg_dipolar_flow.c).

## Code

### General Parameters & Needed Libraries

Physical quantities are non-dimensionalized thanks to the dimensional independent 
triplet quantities $\{\sigma, \rho, L_0 \}$, with $L_0$ the simulation size 
domain:

$$
\mathbf{u} = \sqrt{\dfrac{\sigma}{L_0 \,\rho} } \, \widetilde{\mathbf{u}} 
\quad ; \quad 
t = \sqrt{\dfrac{L_0^3 \, \rho}{\sigma} } \, \tilde{t} 
\quad ; \quad 
p = \dfrac{\sigma}{L_0} \tilde{p}
\quad ; \quad 
\kappa = \dfrac{1}{L_0} \widetilde{\kappa} 
$$
$$
r = L_0 \, \tilde{r}
\quad ; \quad 
z = L_0 \, \tilde{z} 
\quad ; \quad 
\Delta = L_0 \, \widetilde{\Delta} 
\quad ; \quad 
\mu_d = \sqrt{\sigma \, L_0 / \rho} \, \widetilde{\mu}_0
\quad ; \quad 
C = \widetilde{C}
$$

where $C$ is the volume fraction, and $\mu_d$ is the far-field *dipolar flow* 
detailed on [*this page*](http://basilisk.fr/sandbox/cailler/sierou_lister/u_BC_dipolar_flow.h). 
The non-dimensionalization is done with the 
**liquid** density and surface tension. Hence, the non-dimensionalized 
*Navier--Stokes* equations solved with `Basilisk` read as:

$$
\left\{\begin{array}{rcl} 
  \partial_{\tilde{t}} \mathbf{\widetilde{u}} 
  + \left( \mathbf{\widetilde{u}} \cdot \widetilde{\boldsymbol{\nabla}}\right) 
  \mathbf{\widetilde{u}}
  &=&
  - \widetilde{\boldsymbol{\nabla}} \tilde{p}
  + \textcolor{red}{\mathrm{Oh}} \, \widetilde{\boldsymbol{\nabla}}^2 \mathbf{\widetilde{u}} 
  + \widetilde{\kappa} \, \widetilde{\boldsymbol{\nabla}} \widetilde{C} 
  \quad \text{in liquid}\\
  \partial_{\tilde{t}} \mathbf{\widetilde{u}} 
  + \left( \mathbf{\widetilde{u}} \cdot \widetilde{\boldsymbol{\nabla}}\right) 
  \mathbf{\widetilde{u}}
  &=&
  \textcolor{Orchid}{\dfrac{\rho_{liq}}{\rho_{gaz}}} \left(
  - \widetilde{\boldsymbol{\nabla}} \tilde{p} 
  + \textcolor{orange}{\dfrac{\mu_{gaz}}{\mu_{liq}}} \,
  \textcolor{red}{\mathrm{Oh}} \, \, \widetilde{\boldsymbol{\nabla}}^2 \mathbf{\widetilde{u}} 
  + \widetilde{\kappa} \, \widetilde{\boldsymbol{\nabla}} \widetilde{C} 
  \right)
  \quad \text{in gas}
\end{array}
\right.
$$

where $\mathrm{Oh} = \mu_{liq}/\sqrt{\sigma \, \rho_{liq} \, L_0}$ is the 
*Ohnesorge number* for the liquid, that measures viscous effects compared to 
inertia and surface tension. 

The simulations are **inviscid**, therefore 
$\mu_{liq}^{num} = \mu_{gas}^{num} = 0 \Rightarrow \mathrm{Oh} = 0$, and done 
with a density ratio of $10^3$, hence 
$\rho_{liq}^{num} = 1$, $\rho_{gas}^{num} = 10^{-3}$, whereas $\sigma^{num} = 1$. 

Preliminary convergence studies (see p.181 of [Nicolas Cailler's thesis](https://theses.hal.science/tel-04965377v1)) 
have shown that, for a simulation box of size 
unity and a dipolar flow intensity in the range of 
$\widetilde{\mu}_0 \in [-10 \, ; 10]$, the interface shapes are well resolved 
and self-similarly converged when rescaled by $\tilde{t}^{2/3}$ for the 
following optimized parameters: $\tilde{t}_{max} = 2 \times 10^{-2}$, 
$N_{max} = 9 \,\,  \longleftrightarrow \,\, 
\widetilde{\Delta}_{opt} = 1.95 \times 10^{-3}$. 

**However**, the `Basilisk` server running time limitations constrain us to 
decrease the optimum level of refinement for completing the parametric study 
explained hereafter: 

$$
\tilde{t}_{max} = 2 \times 10^{-2}, \quad \quad 
N_{max} = 8 \,\,  \longleftrightarrow \,\, 
\widetilde{\Delta} = 3.9 \times 10^{-3}
$$

This is done to run the [figure output](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_dipolar_flow.c#outputs) 
as expected.
*/

#define LEVEL 6

#define MAXLEVEL 8
#define SIZE 1 
#define T_END 0.02 

#define THETA_0 120.*pi/180.
#define BETA_0 pi - THETA_0
#define MU_0  1

#define Y_OFFSET 0. 
#define X_OFFSET 0.6*SIZE


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

/** 
In this study, only the values $\widetilde{\mu}_0 \in [1 \, ; 4]$ are 
simulated for the far-field dipolar flow, because:

  + they correspond to the values numerically achieved by *Sierou \& Lister*;
  + to speed up the computation on the `Basilisk` server.

The special case of $\widetilde{\mu}_0 = 0$ was treated in local, in order to 
not mess up with the way Boundary Conditions are treated in this file 
(would need a separate run by changing the preprocessor definition `MU_0` to `0`).
*/

int k; // counter
int mu_0[] = {0, 1, 2, 3, 4}; // dipole potential strength

/** 
We then need to import the files implementing the far-field velocities 
developing under the presence of a far-field dipolar flow:
*/
#include "u_BC_dipolar_flow.h" // far-field velocities for MU_0 ≠ 0

/** 
Running the code, some micro-bubbles could emerge on the axis. We remove them, 
as they are not expected:
*/
#include "tag.h" // to remove erroneous micro-bubbles


/** 
### Boundary Conditions 


The general numerical configuration is schematized on the below figure:

![Boundary conditions definitions 
\& geometrical considerations (here for $\widetilde{\mu}_0 = -10$ though)](img_sierou/sierou_BCs.png){width="35%"}

The normal to the interface used for the management of ghost cells' 
volume fractions is given by the $L_1$--normalized vector:

$$
\mathbf{n} = 
[|\cos \theta_0| + |\sin \theta_0| ]^{-1}\, {}^{\mathrm{t}} 
(-\sin \theta_0 ,\, -\cos \theta_0)
$$
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


/** 
The inviscid liquid cone is initialized with an open angle $\theta_0 > 90$°, 
so that the liquid phase occupies the whole right side of the simulation box. 
Two configurations are then possible for the contact angle $\beta_0$ on the 
opposite side of the apex:

+ if the interface crosses the top border, then: $\beta_0 = \pi - \theta_0$;
+ if the interface crosses the left border, then: 
$\beta_0 = \dfrac{3\pi}{2} - \theta_0$.
*/


// Contact Angles:
h.t[left] = contact_angle (pi/2. + BETA_0) ; 
h.t[top] = contact_angle (BETA_0) ;

/** 
On the axis of symmetry $(\tilde{r} = 0$), *symmetry conditions* are logically 
imposed upon velocity:
*/

// On bottom boundary (AXIS of SYMMETRY):
u.n[bottom] = dirichlet (0.);
u.t[bottom] = neumann (0.);

/** 
Then, if a *dipolar flow* does exist:
*/

#if MU_0 != 0
// without far-field flow --> dirichlet(0) on u!!

/** 
Due to the ratio of density $\rho_{gas}^{num}/\rho_{liq}^{num} = 10^{-3}$, we 
consider a zero-pressure *Dirichlet* condition for the boundaries 
$\partial \Omega$ belonging to the gas phase:

$$
\forall \, (\tilde{r},\tilde{z}) \in \partial \Omega, \quad 
{C}(\tilde{r}, \tilde{z}) = 0 \,\, \Rightarrow \,\, 
\tilde{p}_\infty^+(\tilde{r},\tilde{z}) = 0
$$

which defines the *reference pressure* for the simulation. A *Neumann* condition 
for the velocity field is then associated for $C = 0$ to model an *outflow 
condition* in the gas phase: 

$$
\forall \, (\tilde{r},\tilde{z}) \in \partial \Omega, \quad 
{C}(\tilde{r}, \tilde{z}) = 0 \,\, \Rightarrow \,\, 
\boldsymbol{\nabla}_n \cdot 
\widetilde{\mathbf{u}}_\infty^+ (\tilde{r},\tilde{z}) = 0
$$

*/

  // On left boundary, outflow conditions are imposed for the external fluid.
  // WARNING! The interface can jump from top to left border!
  u.n[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uxA(x,y) );
  u.t[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uyA(x,y) );

  // Pressure BC:
  p[left] = (f[] == 0. ? dirichlet (0.) 
    : neumann ( - ( fm.n[0]/alpha.n[0] )
    *( a.n[0] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) );

  pf[left] = (f[] == 0. ? dirichlet (0.) 
    : neumann ( - ( fm.n[0]/alpha.n[0] )
    *( a.n[0] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) );

/** 
On the boundaries of the liquid phase $(C = 1)$, *Dirichlet* conditions are 
imposed on the velocity field according to the formulae defined in 
[`u_BC_dipolar_flow.h`](http://basilisk.fr/sandbox/cailler/sierou_lister/u_BC_dipolar_flow.h). 
The associated pressure is dealt with a *Neumann* condition to take into 
account the far-field velocity whose normal components are not zero:

$$
\forall \, (\tilde{r},\tilde{z}) \in \partial \Omega, \quad 
{C}(\tilde{r}, \tilde{z}) = 1 \,\, \Rightarrow \,\, 
\left\{\begin{array}{rcl}
  \partial_{\tilde{z}} \tilde{p}_\infty^-(\tilde{r},\tilde{z})  
  &=& \dfrac{1}{\rho_{liq}^{num}} 
  \left( \mathrm{f}_{\tilde{z}} 
  - \tilde{u}_{\tilde{z}} \partial_{\tilde{z}} \tilde{u}_{\tilde{z}} 
  - \tilde{u}_{\tilde{r}} \partial_{\tilde{r}} \tilde{u}_{\tilde{z}}\right)  \\
  \partial_{\tilde{r}} \tilde{p}_\infty^-(\tilde{r},\tilde{z})  
  &=& \dfrac{1}{\rho_{liq}^{num}} 
  \left( \mathrm{f}_{\tilde{r}} 
  - \tilde{u}_{\tilde{r}} \partial_{\tilde{r}} \tilde{u}_{\tilde{r}} 
  - \tilde{u}_{\tilde{z}} \partial_{\tilde{z}} \tilde{u}_{\tilde{r}}\right)
\end{array}
\right.
$$

where $\mathbf{f}$ refers to the volume forces. 
*/


  // On top boundary, exterior and interior fluids coexist for most of the
  // simulation:
  u.n[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uyA(x,y) );
  u.t[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uxA(x,y) );

  p[top] = ( f[] == 1. ?
    neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uyA(x,y)*duyA_dy(x,y) - uxA(x,y)*duyA_dx(x,y)) )
    : dirichlet (0.) );

  pf[top] = ( f[] == 1. ?
    neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uyA(x,y)*duyA_dy(x,y) - uxA(x,y)*duyA_dx(x,y)) )
    : dirichlet (0.) );


  // On right boundary (will always be liquid phase):
  u.n[right] = dirichlet ( uxA(x,y) );
  u.t[right] = dirichlet ( uyA(x,y) );
  
  p[right] = neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) ;

  pf[right] = neumann ( ( fm.n[ghost]/alpha.n[ghost] )
    *( a.n[ghost] - uxA(x,y)*duxA_dx(x,y) - uyA(x,y)*duxA_dy(x,y)) ) ;

/** 
If no *dipolar flow* is considered, imposing a zero-velocity on all borders 
(excluding the axis of symmetry) gives very close results compared to the 
slightly more complex considerations done 
[*there*](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_sierou_all.c#boundary-conditions) 
(direct comparisons with the results of *Sierou \& Lister* tend to confirm 
this simple choice, albeit not rigorous):
*/

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

  // loop for exploring different μ0 values:
  for (k = 1; k < 5; k++){    
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


/** 
The initial cone is the "ultimate" self-similar form: the initial condition is 
the *exact* field for the cone; in other words, the velocity field defined 
in [`u_BC_dipolar_flow.h`](http://basilisk.fr/sandbox/cailler/sierou_lister/u_BC_dipolar_flow.h) 
is valid in the whole domain:
*/

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

/** 
We update *height functions* to compute with precision the curvature of the 
interface. This update might be non necessary, but we ensure that way that 
our modifications in ghost cells' volume fractions have been well taken into 
account at the moment of computing the curvature at crossed boundaries.
*/
  // computation of curvature for visualization:
  boundary({f});

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});

  curvature (f, curv_viz, 1, add = false);
}

/** 
Vorticity and curvature are computed at each timestep:
*/
event vorti (i++){
  // cf. bug report: https://groups.google.com/g/basilisk-fr/c/ok-OhtzO1Pk
  // vorticity in AXI:
  foreach()
    omega[] = y*(u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);;
  boundary({omega});

  // curvature
  curvature (f, curv_viz, 1, add = false);
}

/** 
We use the adaptive mesh refinement of `Basilisk`, in order to capture the 
various shapes at the apex of the cone under the influence of the 
*dipolar flow strength*:
*/
event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1e-3,1e-2,1e-2}, MAXLEVEL, LEVEL);
    // unrefine(f[] == 0);
  #endif 


/** 
In the far-field, there exists a movement of the interface resulting from 
the balance between the *capillary flow* coming from the *Laplace* pressure 
gradient and the *dipolar flow* (`MU_0 != 0`). In that case, it is mandatory 
to update at each timestep the volume fractions in ghost cells belonging to 
the borders where the interface is moving/is susceptible to move. 

The update shall be *dynamical* to allow the translation of the interface 
along the borders, that is to say, feeding the auxiliary volume fraction field 
`f0[]` by applying the functions `f_BC_myborder()` to `f[]` 
(not to `f0[]` which would freeze any free displacement at contact angle 
$\beta_0$; we do not want a pinned point when `MU_0 != 0`). 

Otherwise, when `MU_0 = 0` (no dipolar flow), we apply a *static* update (no 
displacement expected at the contact angle $\beta_0$).
*/  
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

/** 
We try to remove all parasitic micro bubbles forming:
*/

  remove_droplets (f, minsize = 1, threshold = 1.e-3);
  remove_droplets (f, minsize = 1, threshold = 1.e-3, bubbles=true);

  // Trick to be sure that all parasitic bubbles will be removed:
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





/** 
### Outputs
*/

event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // vider le buffer
}

event end (t = T_END){}

/** 
At the end of each run, we store the shape of the interface into a dedicated 
file:
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
We are now able to plot all the `Basilisk` curves against *Sierou \& Lister* 
ones, for $\widetilde{\mu}_0 \in [0 \, ; 4]$ and providing from local 
computation the special case of no dipolar flow 
(*cf.* reasons given [*here*](file:///Users/nicolas/wiki/sandbox/cailler/sierou_lister/sierou_dipolar_flow.c.html#general-parameters-needed-libraries)):

~~~pythonplot Parametric study for $\theta_0 = 120$° 
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
N = 8
mu0_list = range(0,5,1)

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
ax.set_ylim(-4,6)
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


# Labels for μ0 = 0:
plt.text(6., 3.8, r'$\widetilde{\mu}_0 = 0$', rotation=30, color=mapcolors[0], fontsize=12)
plt.text(6., 2.2, r'$\widetilde{\mu}_0 = 4$', rotation=34, color=mapcolors[8], fontsize=12)

# Plot Sierou curves for θ_0 = 120°:
for k in range(5,10):
  sierou = np.loadtxt("../sierou_data/fig16b_sierou_curve{}.dat".format(k))
  if k == 5:
    ax.plot(sierou[:,0], sierou[:,1], 
            '--', lw = 0.8, color='darkred',
            label=r'Sierou & Lister')
  else:
    ax.plot(sierou[:,0], sierou[:,1], '--', lw = 0.8, color='darkred')

ax.legend(frameon=False, loc=1)

plt.savefig('fig16b_repro_th120_t0.02_N8_vert_axis_mu0_pos.svg') 
~~~



*/


/** 
## Discussion

A full study -- not suited for the `Basilisk` server -- has been conducted 
locally for two different angles ($\theta_0 = 120$° and $145$°) and a 
broad range of dipolar flow strengths ($\widetilde{\mu}_0 \in [0 \, ; 10]$), 
with optimized parameters for a better refinement, without sacrifying 
too much computational speed: 
$\tilde{t}_{max} = 2 \times 10^{-2}$, 
$N_{max} = 9 \,\,  \longleftrightarrow \,\, 
\widetilde{\Delta}_{opt} = 1.95 \times 10^{-3}$. 

The results are gathered on the below figures:

![`Basilisk` results for $\theta_0 = 120$° (*left*) and $\theta_0 = 145$° (*right*, 
not done by *S&L*)](img_sierou/sierou_mu-pos_v2.png){width="50%"}

+ For $\theta_0 = 120$°, $\widetilde{\mu}_0 \in [0 \, ; 4]$, the results in 
the scale invariant space of [Sierou \& Lister, (2004)](#sierou2004) are 
almost indistinguishable from our simulations with `Basilisk`. These results 
confirm that in the self-similar space, the interface position depends on the 
strength of the dipolar flow. 

+ Interface profiles widely differ given the far-field angle $\theta_0$, with 
notably visible capillary waves for $\theta_0 = 145$° (large open angles) 
depleting rapidly with increasing parameter $|\widetilde{\mu}_0|$. 

+ Combined with our results for $\theta_0 = 145$°, 
it shows that an infinity of self-similar solutions exists, indexed by 
$(\theta_0, \widetilde{\mu}_0)$.

+ Without dipolar flow $(\widetilde{\mu}_0 = 0)$, an obtuse cone of liquid
$(\theta_0 > 90$°$)$ "moves forward" the gas cavity, whereas an acute liquid 
cone $($[$\theta_0 < 90$°](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_sierou_all.c)$)$ 
*recoils*. Both phenomena are due to the overall movement of the *capillary flow*. 

+ A *dipolar flow* for $\widetilde{\mu}_0 > 0$ tends to recoil the cone and 
counteracts the *capillary flow* moving forward the (obtuse) cone of liquid.  

+ The effective recoil threshold (for obtuse liquid cones) is obtained when 
the dipolar-flow-generated velocity field overcomes the *Laplace* pressure 
gradient flow (*capillary flow*): on the $\theta_0 = 120$°--figure, this 
threshold value is obtained for $\widetilde{\mu}_0^* = 2$, since the abscissa 
of the apex is *negative*, which demonstrates the complete recoil of the 
interface compared to its initial state. 
However, for $\theta_0 = 145$° (above figure) the effective recoil of the 
obtuse liquid cone is only achieved when $\widetilde{\mu}_0^* = 7$. 
With larger open angles, greater dipolar flow strengths are needed to observe 
a reversal of the translation of the obtuse cone along the axis of symmetry.

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

@article{zeff2000,
  title={Singularity dynamics in curvature collapse and jet eruption on a fluid surface},
  author={Zeff, Benjamin W and Kleber, Benjamin and Fineberg, Jay and Lathrop, Daniel P},
  journal={Nature},
  volume={403},
  number={6768},
  pages={401--404},
  year={2000},
  url = {https://complex.umd.edu/papers/nature2000.pdf},
  doi = {doi: 10.1038/35000151},
  publisher={Nature Publishing Group UK London}
}

@article{duchemin2002,
  title={Jet formation in bubbles bursting at a free surface},
  author={Duchemin, Laurent and Popinet, St{\'e}phane and Josserand, Christophe and Zaleski, St{\'e}phane},
  journal={Physics of fluids},
  volume={14},
  number={9},
  pages={3000--3008},
  url = {https://laurentduchemin.gitlab.io/Papers/dpjz.pdf},
  year={2002},
  doi = {https://doi.org/10.1063/1.1494072},
  publisher={American Institute of Physics}
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

@article{poujol2021,
  title={Sound of effervescence},
  author={Poujol, Mathis and Wunenburger, R{\'e}gis and Ollivier, Fran{\c{c}}ois and Antkowiak, Arnaud and Pierre, Juliette},
  journal={Physical Review Fluids},
  volume={6},
  number={1},
  pages={013604},
  year={2021},
  url ={https://hal.science/hal-03121607/document}, 
  doi = {https://doi.org/10.1103/PhysRevFluids.6.013604},
  publisher={APS}, 
}


~~~
*/