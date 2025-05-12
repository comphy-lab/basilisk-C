/** 

# Self-Similar DNS convergence study for the *Keller \& Miksis (1983)* problem

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**DISCLAIMER** 

Due to inherent running time constraints of the server, the simulations for 
generating the graphs and movies were done at ***low** resolution* and a 
*shorter end time* than wished. 
As a direct consequence, spurious vorticity pockets along the interface are 
disturbing the interface due to this low-resolution, and on the long run 
$(\tau > 10)$, a divergence of the position of the interface is hence observed. 

**This is not the case when running with `MAXLEVEL = 9`, where the self-similar 
solver really behaves as expected, no matter the time of the simulation, 
and vorticity disturbances are not seen on movies.**
</div>


Reproduction in the self-similar domain of the [Keller \& Miksis, (1983)](#keller1983) 
wedge recoil driven by surface tension. 
Extension also made for the simple axisymmetric case of a recoiling 
cone in vacuum without any dipolar far-field flow $(\mu_0 = 0)$, 
as described in [Sierou \& Lister, (2004)](#sierou2004). 



![Wedge recoil with the self-similar solver developed](selfsim_keller_conv/selfsim_keller_N8_movie.mp4)(width="300" height="500")




## Motivations

The shapes of the interface shall match the Fig. 2 of 
[Keller \& Miksis, (1983)](#keller1983) when **directly simulated 
into the self-similar space**, 
for the following initial wedge angles:
$\theta_0 = 27.5$°, $32.5$°, $45$°, $65$° and $80$°. 

For the AXI case, the interface shall match the Fig. 2 of 
[Sierou \& Lister, (2004)](#sierou2004) when **directly simulated 
into the self-similar space**, 
for the following initial cone angles:
$\theta_0 \text{ [rad] } = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2.$

## Methodology (Recap)

Thanks to the [2D results of Keller \& Miksis in the physical space](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c), 
we know the self-similar space-extension of the profiles. 
Therefore, a $12 \times 12$ box is sufficient for both avoiding border effects 
and allowing the wedge tip to reach its steady $\xi$--position.

Viscosity is considered here (*cf. infra*).

For volume fractions boundary conditions, the applied method is the same as 
the one used in the physical space (`f0` auxiliary volume fraction field 
storing the initial state, and a better management of normals at borders). 
 
<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**For detailed explanations upon the self-similar solver: 
*[go to this page](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README)!***
</div>

In short, the self-similar solver is composed of the following additions:

  + in the file [`selfsim_centered_keller.h`](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_centered_keller.h):

    - new variables: `vector xi[]` | `face vectors uf_temp[], lambdaf[]`;
    - boundary conditions modified for the pressure for symmetry conditions;
    - boundary conditions for the vector `xi[]` and face vector `lambdaf[]` 
      (`uf_temp[]` is a *temporary* face vector for storing the $n+1/2$ state of the 
      predicted-projected face velocity `uf[]`, and as such, no need for 
      securing the BCs);
    - **init event:** position vector `xi[]` and advection face vector 
    $\mathbf{\overline{\Lambda}} = \mathbf{\overline{u}} - (2/3) \boldsymbol{\xi}$ 
    are defined;
    - **stability event:** modified with the new advection velocity $\mathbf{\overline{\Lambda}}$; 
    
    - **advection_term event:** 
    
      1. `uf[]` is predicted with an upwinded state of the BCG algorithm taking 
      into account the new advection velocity $\mathbf{\overline{\Lambda}}$ 
      and additional source terms;
      2. `uf[]` is then projected onto the divergence-free space;
      3. `lambdaf[]` is built with the predicted-projected `uf[]` face vector 
      and a linear interpolation to the faces of the position vector;
      4. `uf_temp[]` is built for storing the $n+1/2$ state of the 
        predicted-projected face velocity `uf[]`, as it will be used later 
        in the acceleration event for evaluating source terms;
      5. fluxes are then computed using the BCG algorithm (slightly modified 
        with source terms, *cf. infra*): material points are advected thanks to 
        the new advection velocity $\mathbf{\overline{\Lambda}}$;
        
    - **acceleration_term event:** 
    
      1. add `dt*((1. - 2.*Nd)/3.)*uf_temp.x[]` as a $n+1/2$ source term;
      2. `Nd` = number of dimensions ($2$ in 2D / $3$ in AXI);
      
    - **projection event:** 
    
      1. do not forget to add again in the `centered_gradient()` function 
          the additionnal source term `((1. - 2.*Nd)/3.)*uf_temp.x[]`;
      2. `lambdaf[]` face vector is re-built as `uf[]` is projected to 
      determine the $n+1$ state;

  + in the file [selfsim_bcg_keller.h](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_bcg_keller.h):
    - supplementary source terms have to be taken into account, as the original 
      BCG algorithm assumes divergence-free advection velocity 
      [and so $\boldsymbol{\nabla} \cdot \left(\mathbf{u} \otimes \mathbf{u} \right) 
      = \left(\mathbf{u} \cdot \boldsymbol{\nabla} \right) \mathbf{u}$, 
      which is not the case here, as we do have 
      $\boldsymbol{\nabla} \cdot \left(\mathbf{\overline{\Lambda}} \otimes 
      \mathbf{\overline{u}} \right) = \left(\mathbf{\overline{\Lambda}} \cdot 
      \boldsymbol{\nabla} \right) \mathbf{\overline{u}} + \mathbf{\overline{u}} 
      \left(\boldsymbol{\nabla} \cdot \mathbf{\overline{\Lambda}}\right) = 
      \left(\mathbf{\overline{\Lambda}} \cdot 
      \boldsymbol{\nabla} \right) \mathbf{\overline{u}} - (2 \, N_d / 3) 
      \mathbf{\overline{u}}$;
      therefore we do need to add the following line for the computation of 
      the re-predicted face velocity:
      
```C
      double f2 = f[i] + ( (src[] + src[-1]) + (2.*Nd/3.)*(f[] + f[-1]) )*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;
```

  + in the file [selfsim_two-phase-generic_keller.h](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_two-phase-generic_keller.h):
    - viscosity has to "vanish" exponentially in time with a modified 
      non-dimensional viscosity:
      
```C 
      muv.x[] = alphav.x[]*pow(1./(alphav.x[]*f.sigma), 2./3.)
                *fm.x[]*mu(ff)*exp(-t/3.);
```

  + in the file [selfsim_vof_keller.h](http://basilisk.fr/sandbox/cailler/self_sim_DNS/selfsim_vof_keller.h): 
  use $\mathbf{\overline{\Lambda}}$ now as an external face vector.
      


*/



/** 
## Code

### General Parameters 
*/


#define LEVEL 7 
#define MAXLEVEL 8
#define SIZE 12e0
#define T_END 50 
#define X_OFFSET (1./3.)*SIZE

// #include "axi.h"
#if AXI // Sierou & Lister (2004) value
  #define THETA_0 0.8 
#else // Keller & Miksis (1983) value
  #define THETA_0 45*pi/180.
#endif

#include "selfsim_centered_keller.h"
#include "contact.h"
#include "selfsim_two-phase_keller.h"
#include "tension.h"
#include "selfsim_f_BC_keller.h"



// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 

scalar f0[]; // auxiliary volume fraction field used for BC 
vector n_front[]; // normal vector of the interface
scalar alpha_front[]; // intercept linked to n_front[] during VoF-reconstruction 


/** 
### Boundary Conditions
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


#if AXI // capillary flow due to the curvature gradients in AXI:
  u.t[right] = dirichlet ( y/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
  u.n[right] = dirichlet ( x/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
  u.n[top] = dirichlet ( y/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
  u.t[top] = dirichlet ( x/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
#else // No flow on the right and top borders in 2D:
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
#endif


// Outflow on the left gas-border (to set a pressure reference):
// [ONLY in AMR, since we need low-resolution on this border to avoid backflow]
// #if TREE
//   p[left] = dirichlet(0.);
//   pf[left] = dirichlet(0.);
//   u.n[left] = neumann(0.);
// #endif


/** 
### Generic Events
*/

int main() {
  size(SIZE) ;
  init_grid(1 << LEVEL) ;
  origin(-X_OFFSET, 0) ;

  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 1e-3, mu2 = 1e-5; 

  run() ;
}


event init (t = 0){
  // Initial free-surface definition:
  #if TREE
    do{
    fraction(f0, - (y - x*tan(THETA_0)) ) ; 
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, MAXLEVEL).nf); 

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
  // cf. bug report: https://groups.google.com/g/basilisk-fr/c/ok-OhtzO1Pk
  // vorticity:
  #if AXI
    foreach()
      omega[] = y*(u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    boundary({omega});
  #else // 2D-case
    vorticity (u, omega);
  #endif

  curvature (f, curv_viz, 1, add = false);
}


event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1.e-3,1e-2,1e-2}, MAXLEVEL, 6);
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


event end (t = T_END){}


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


/** 
### Outputs

We consider the simulation converged if, independently of the number of 
elements in the mesh, the simulation taking place in the self-similar space 
has reached a **steady state** that is equivalent to the interface shapes 
obtained numerically by [Keller \& Miksis, (1983)](#keller1983). We choose 
an initial angle of $\theta_0 = 45$° because it is the most disadvantageous case 
since there is no preferred direction for computing *height functions* 
(the interface is not sufficiently horizontal nor vertical, hence both 
directions have to be explored to estimate the distance to the interface and 
consequently the curvature). 

For the sake of simplicity for running this code on the `Basilisk` server, 
we limit our study to $\tau_{end} = 50$ and $N_{max} = 8 \Rightarrow 
\Delta = 0.047$ for the maximum level of refinement of the grid, since 
preliminary studies have shown that convergence was reached as early as this 
maximum level of refinement. 

*/


event bottom_pos (t += 0.1){
  char filecurv[200] ;
  #if AXI
    sprintf(filecurv, "selfsim_sierou_xi-pos_N%d_th%g.dat", MAXLEVEL, THETA_0*180/pi) ;
  #else
    sprintf(filecurv, "selfsim_keller_xi-pos_N%d_th%g.dat", MAXLEVEL, THETA_0*180/pi) ;
  #endif
  static FILE * fcurv = fopen(filecurv, "a") ;

  double xpos = 0., curv = 0.;
  foreach_boundary (bottom){
    if (interfacial (point,f)){ // if we are at the interface bottom
      if ( (f[] > 0.1) && (f[] < 0.9) ){
        xpos = x;
        curv = curv_viz[]; //avg_neighbor (point, curv_viz, f);
      }
    }
  }
  fprintf(fcurv, "%g %g %g\n", t, xpos, curv);
  fflush(fcurv); // empty the buffer
}

/** 
~~~pythonplot Time convergence of the axis position
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
label_size = 15
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
t_N8, xi_N8, curv_N8 = np.loadtxt("selfsim_keller_xi-pos_N8_th45.dat", unpack=True) 
t_N8_filt = t_N8[xi_N8 > 0]
xi_N8_filt = xi_N8[xi_N8 > 0]


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 18
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_xlim(0,10)
ax.set_ylim(0,1)
ax.set_xlabel(r'$\tau$')
ax.set_ylabel(r'$\xi_{axis}$')

ax.hlines(0.7369596975356669, 0, 10, linestyles='dashed', linewidth=0.8,
          colors = "darkred", label=r'Keller & Miksis')
ax.hlines(0.715, 0, 10, linestyles='dotted', linewidth=0.8,
          colors = "mediumpurple", label=r'Basilisk in physical space')
ax.plot(t_N8_filt, xi_N8_filt, lw=1.2, label=r'$N_{max}=8$', 
        color=mapcolors[9])

ax.legend(frameon=False, loc=4)
plt.savefig('selfsim_keller_xi-pos_conv.svg') 
~~~ 

As shown on the above figure, there is a **quick convergence** for the shape 
of the interface, as soon as $\tau = 4$, by comparing the interface position 
on the axis with the results collected from [Keller \& Miksis, (1983)](#keller1983) 
and [our previous results in the physical space](http://basilisk.fr/sandbox/cailler/keller_miksis/keller_fig2_conv.c#outputs).


We then plot various shape interfaces up to $\tau = 50 \Leftrightarrow 
\mathrm{e}^{50} \approx 5 \times 10^{21}$ in order to show the **efficiency 
of the self-similar solver**: 
*/


event profiles (t = {1., 5., 10., 20., 50.}){                   
  char filename[200] ;
  #if AXI
    sprintf(filename, "selfsim_sierou_shape_th%g_t%g_N%d.dat", 
            THETA_0*180/pi, t, MAXLEVEL) ;
  #else
    sprintf(filename, "selfsim_keller_shape_th%g_t%g_N%d.dat", 
      THETA_0*180/pi, t, MAXLEVEL) ;
  #endif
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}


/** 
~~~pythonplot Convergence towards a scale invariant state with the self-similar DNS solver
times_str = ["1", "5", "10", "20", "50"] 
times_float = np.array([1, 5, 10, 20, 50])


#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
keller_45deg = np.loadtxt ("../../keller_miksis/keller_data/raw_keller_fig2_45.dat", 
                             unpack=False)

#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 12
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,4)
ax.set_ylim(0,4)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')

j = 0
for k, t in enumerate(times_str):
  facets = np.loadtxt ("selfsim_keller_shape_th45_t{}_N8.dat".format(t))
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments):
    if l == 0:
      ax.plot(segment[:, 0], segment[:, 1], 
              lw=1.2, color = mapcolors[k+1+j], label=r'$\tau = {}$'.format(t))
    else:
      ax.plot(segment[:, 0], segment[:, 1], 
              lw=1.2, color = mapcolors[k+1+j]) 
  j += 1

# Keller & Miksis results:
ax.plot(keller_45deg[:,0], keller_45deg[:,1], '--', lw=0.8, color="darkred",
        label=r'Keller & Miksis')


ax.legend(frameon=False, loc=2)

plt.savefig('selfsim_solver_keller_45deg_conv_N8.svg') 
~~~

Indeed, the interface is *steady* as early as $\tau \sim 4$ and remains in this 
steady state **until the end of the simulation**. To keep things in perspective, 
it means that during a total running time of $\tau = 50$, the wedge recoils 
over a distance $L \sim (\mathrm{e}^{50})^{2/3} \sim 3 \times 10^{14}$... To keep a 
resolution of $\Delta = 0.047$ used for the self-similar DNS, it would need 
in the physical space a maximum level of refinement of $N_{max}^{phys} = 53 \,\,
(3 \times 10^{14}/2^{53} \sim 0.033)$, so **53 levels of refinement**!! 
This would be absolutely not possible to reproduce in the physical space 
(*the design would not be very human*).

***We now understand the relevance of such a self-similar solver for scale 
invariant problems.***

*/

event vel_p_maps (t = {1., 5., 10., 20., 50.})
{
  scalar x_interf[];
  scalar y_interf[];
  char fileup[200];
  char fileinterf[200];
  #if AXI
    sprintf (fileup, "selfsim_sierou_u_p_th%g_t%g.dat", THETA_0*180/pi, t);
    sprintf (fileinterf, "selfsim_sierou_interf_th%g_t%g.dat", THETA_0*180/pi, t);
  #else
    sprintf (fileup, "selfsim_keller_u_p_th%g_t%g.dat", THETA_0*180/pi, t);
    sprintf (fileinterf, "selfsim_keller_interf_th%g_t%g.dat", THETA_0*180/pi, t);
  #endif
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

event viewing (t +=0.05) {
  view (width = 600, height = 1000, fov = 25, tx = -0.25, ty = 1e-4,
	bg = {.7, .7, .7});

  clear();
  draw_vof ("f", lw = 1.5);
  squares ("omega", linear = true, map = blue_white_red);
  box (notics = false);
  mirror ({0,1}) {
    draw_vof ("f", lw = 1.5);
    squares ("p", min = -1., max = 1., linear = false, map = cool_warm); 
    box (notics = true);
  }
  char legend[100];
  sprintf (legend, "ln t = %0.2g", t);
  draw_string (legend, 1, size = 20., lw = 1.7); // “0” bl, “1” tl, “2” tr and “3” br 

  #if AXI
    save ("selfsim_sierou_N8_movie.mp4");
  #else
    save ("selfsim_keller_N8_movie.mp4");
  #endif

}

/** 
### Role of the viscosity

As seen on the movie, at the beginning of the simulation in the scale invariant 
space $(\tau = 1)$ we can see a thick boundary layer detachment 
in the regions of strong curvature variations and diffuses in the domain. 
This can be understand by the fact that the viscous term is still relevant 
$(\textcolor{red}{\mathrm{e}^{-\tau/3} \sim 0})$. 

However, for the steady regime $\tau = 10 \sim 50$, these boundary layers 
are thinning and are closer to the interface: viscous diffusion still applies, 
but the **radial dilatation** exerts by the velocity advection term dominates 
this diffusion, along with the fact that 
$(\textcolor{red}{\mathrm{e}^{-\tau/3} \xrightarrow[\tau \to + \infty]{} 0})$. 



*/



/** 
## Technical Remarks and Empirical Observations

DO NOT PUT `p[border] = dirichlet (0.);` with the self-similar solver!!!!
Enormous vorticity layers are arriving from this border and destroy 
everything.
At the moment, default symmetry conditions on the gas border seem fine.

The following remarks are important for overcoming unexpected behaviours 
when using mesh refinement:

* In AMR, the following patches can be used:

  + if we still want to use the outflow condition `p[left] = dirichlet (0.);`, 
    then the (for instance) instruction `unrefine(x < -7./24.*SIZE);` will 
    limit the backflow, but this is highly not recommended;

  + there are spurious velocity trails on the axis that can be diminished 
    thanks to `unrefine( ((f[] == 1.) || (f[] == 0)) && (y < 10.*Delta) );`

* In AXI: for better accuracy of the shape interface position on the axis
    compared to the [Sierou & Lister (2004)](#sierou2004) results, 
    far-field flow is needed as boundary conditions:

    ```C
      u.t[right] = dirichlet ( y/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
      u.n[right] = dirichlet ( x/(pow(x*x + y*y, 3./2.) * (tan(THETA_0))) );
    ```
    
    These boundary conditions are expected, as curvature gradients that exist 
    for cones in AXI are creating capillary flows.
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