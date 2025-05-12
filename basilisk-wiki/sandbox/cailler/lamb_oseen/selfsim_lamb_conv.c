/**
# Convergence study of the *Lamb--Oseen* vortex simulated with the self-similar solver 

## Motivations

We want to recover the analytical self-similar solution of the [*Lamb--Oseen 
vortex*](https://en.wikipedia.org/wiki/Lamb%E2%80%93Oseen_vortex) 
when the initial velocity is already the exact solution, in order 
to verify the steady state achieved by the self-similar solver, and 
also the numerical order of convergence (*is it still* $2^{\mathrm{nd}}$ *order or* 
$1^{\mathrm{st}}$ *order accurracy?*). 
A recap of all the main analytical results is summed up in the ending 
[*Appendix*](#appendix-analytical-solutions-in-cartesian-coordinates).

If the self-similar solver works as intended, its numerical solution 
shall remain *steady* for the whole simulation since the initial guess is 
already the correct one. To ensure that everything works well, we first measure 
at a given time of the simulation the difference with the exact solution for 
different discretization levels $N_{max}$ of the grid (*spatial convergence*). 
Then, at a satisfying resolution, we measure the *error diffusion* over time, 
since we should stay on the initial condition solution. 


## Performing a self-similar DNS: *a recap*

Thanks to the [results obtained in the physical space](http://basilisk.fr/sandbox/cailler/lamb_oseen/lamb.c), 
we know the self-similar space-extension of the profiles for appreciable 
reproduction of the self-similar solution: $[-20 \, ; 20]$. 
Therefore, a box *five times larger* should be enough for avoiding border effects, 
that is to say, a square box $[-100 \, ; 100] \times [-100 \, ; 100]$ of size 
$L_0 = 200$ centered on the origin of the Cartesian coordinate system.

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**For detailed explanations upon the self-similar solver: 
*[go to this page](http://basilisk.fr/sandbox/cailler/self_sim_DNS/README)!***
</div>

In short, the self-similar solver for the *Lamb--Oseen* problem is composed of 
the following additions:

  + in the file [`selfsim_centered_lamb.h`](http://basilisk.fr/sandbox/cailler/lamb_oseen/selfsim_centered_lamb.h):

    - new variables: `vector xi[]` | `face vectors uf_temp[], lambdaf[]`;
    - boundary conditions modified for the pressure for symmetry conditions;
    - boundary conditions for the vector `xi[]` and face vector `lambdaf[]` 
      (`uf_temp[]` is a *temporary* face vector for storing the $n+1/2$ state of the 
      predicted-projected face velocity `uf[]`, and as such, no need for 
      securing the BCs);
    - **init event:** position vector `xi[]`and advection face vector 
    $\mathbf{\overline{\Lambda}} = \mathrm{Re} \, \mathbf{\overline{u}} - (1/2) \boldsymbol{\xi}$ 
    are defined;
    - **stability event:** modified with the new advection velocity $\mathbf{\overline{\Lambda}}$; 

    - **advection_term event:** 

      1. `uf[]` is predicted with an upwinded state of the BCG algorithm taking 
        into account the new advection velocity $\mathbf{\overline{\Lambda}}$ 
        and additional source terms;
      2. `uf[]` is then projected onto the divergence-free space;
      3. `lambdaf[]` is built with the predicted-projected `uf[]` face vector and 
        a linear interpolation to the faces of the position vector;
      4. `uf_temp[]` is built for storing the $n+1/2$ state of the 
        predicted-projected face velocity `uf[]`, as it will be used later 
        in the acceleration event for evaluating source terms;
      5. fluxes are then computed using the BCG algorithm (slightly modified 
        with source terms, *cf. infra*): material points are advected thanks to 
        the new advection velocity $\mathbf{\overline{\Lambda}}$;

    - **acceleration_term event:** add `-.5*dt*uf_temp.x[]` as a $n+1/2$ source term;

    - **projection event:**

      1. do not forget to add again in the `centered_gradient()` function 
          the additional source term `-.5*uf_temp.x[]`;
      2. `lambdaf[]` face vector is re-built as `uf[]` is projected to determine 
          the $n+1$ state;

  + in the file [`selfsim_bcg_lamb.h`](http://basilisk.fr/sandbox/cailler/lamb_oseen/selfsim_bcg_lamb.h):
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
      \boldsymbol{\nabla} \right) \mathbf{\overline{u}} - \mathbf{\overline{u}}$; 
      therefore we do need to add the following line for the computation of 
      the re-predicted face velocity:

```C
double f2 = f[i] + ( (src[] + src[-1]) + 1.*(f[] + f[-1]) )*dt/4. 
            + s*(1. - s*un)*g.x[i]*Delta/2.;
```
*/


/** 
## Code

### General Parameters 

The convergence study is done with a regular Cartesian mesh. The total running 
time of the simulation is $\tau_{end} = 10$, large enough for evaluating the 
self-similar solver behaviour (how evolves the difference with the analytical 
solution). 
*/


#include "grid/multigrid.h"
#include "selfsim_centered_lamb.h"

#define LEVEL 8 // minimum grid size for AMR
#define MAX_LEVEL 9 // maximum refinement level for AMR
#define SIZE 200
#define T_END 10

scalar uy_sol[]; // analytical self-sim. sol. for the η-component of the velocity 
scalar uy_diff[]; // difference between the numerical result and the analytical sol.
scalar omega[]; // vorticity 


/** 
Several discretization levels $N_{max}$ are tested:

$$
\forall \, N_{max} \in \{6 \, ; 10 \}
\,\, \Rightarrow \,\, 
\Delta = L_0 / 2^{N_{max}} \in \{0.195 \, ; 3.125 \}
$$

However, the case $N_{max} = 10$ is done remotely, as it would exceed the 
`Basilisk` server capacities. Related results can be found in 
[`lamb_data`](http://basilisk.fr/sandbox/cailler/lamb_oseen/lamb_data/). 
*/

int k; // counter
int lvl_grid[] = {6, 7, 8, 9};


/** 
The initial velocity field is the analytical self-similar solution of 
circulation unity:

$$
\overline{u}_\theta^\star(\xi,\eta) 
= \dfrac{1}{2 \pi \sqrt{\xi^2 + \eta^2}} 
\left(1 - \mathrm{e}^{- (\xi^2 + \eta^2) / 4} \right)  
$$
*/

double velocity_profile(double z){
  /* Exact analytical self-similar solution for η = 0 */
  return ((1. / (2. * pi * z)) * (1. - exp(-(z * z + 0 * 0)/4.)));
}


/**
### Boundary conditions 

To close the problem, the vortex circulation must be preserved when changing 
to the self-similar coordinates. The following constraint shall be therefore 
considered:

$$
\displaystyle \oint_{\mathcal{C}} 
\mathbf{u} \cdot \mathbf{t} \,\, \mathrm{d} \ell = \Gamma
\quad \Rightarrow \quad
\displaystyle \oint_{\mathcal{C}} 
\mathbf{\overline{u}} \cdot \mathbf{t} \,\, \mathrm{d} \xi_\ell = 1
$$

where $\xi_\ell$ is a contour element in the scale invariant space. 
Since numerical simulations are done in a square box of size $L_0$, we need to 
pay attention to boundary conditions as there is no rotation invariance of the 
domain about the origin. Consequently, *Dirichlet* conditions optimally 
complying the previous constraint have to be addressed: 
far from the vortex core, we can then suppose a *potential* flow, hence, 
*irrotational*. Under this assumption, there exists a constant $\alpha$ for 
which $\mathbf{u} = \alpha/r \, \mathbf{e}_\theta$, *i.e.*, in self-similar 
Cartesian variables:

$$
\mathbf{\overline{u}}
=
{}^{\raisebox{2pt}{\scriptsize $\mathrm{t}$}}\!
\begin{pmatrix}
  \overline{u} & \overline{v}
\end{pmatrix} 
=
\alpha \, 
{}^{\raisebox{5pt}{\scriptsize $\mathrm{t}$}}\!
\begin{pmatrix}
  - \dfrac{\eta}{\xi^2 + \eta^2} 
  &
  \dfrac{\xi}{\xi^2 + \eta^2}
\end{pmatrix}
$$

In that case, the circulation constraint reads as:

$$
\displaystyle \oint_{\mathcal{C}} 
\mathbf{\overline{u}} \cdot \mathbf{t} \,\, \mathrm{d} \xi_\ell = 1
\Rightarrow
4 \, \alpha \displaystyle \int_{-L_0/2}^{L_0/2} 
  \dfrac{L_0/2}{{\xi'}^2 + {L_0}^2/4} \,\, \mathrm{d} \xi' 
= 1 
\Rightarrow
\alpha = \dfrac{1}{2 \, \pi}
$$

so finally:

$$
  \overline{u} = - \dfrac{1}{2 \pi} \dfrac{\eta}{\xi^2 + \eta^2}
  \quad ; \quad
  \overline{v} = \dfrac{1}{2 \pi} \dfrac{\xi}{\xi^2 + \eta^2}
$$
*/

u.t[right] = dirichlet ((1./(2.*pi)) * (x)/(x*x + y*y)) ;
u.t[left] = dirichlet ((1./(2.*pi))  * (x)/(x*x + y*y)) ;
u.t[top] = dirichlet (- (1./(2.*pi)) * (y)/(x*x + y*y)) ;
u.t[bottom] = dirichlet (- (1./(2.*pi)) * (y)/(x*x + y*y)) ;

u.n[right] = dirichlet (- (1./(2.*pi)) * (y)/(x*x + y*y)) ;
u.n[left] = dirichlet (- (1./(2.*pi)) * (y)/(x*x + y*y)) ;
u.n[top] = dirichlet ((1./(2.*pi)) * (x)/(x*x + y*y)) ;
u.n[bottom] = dirichlet ((1./(2.*pi)) * (x)/(x*x + y*y)) ;

/** 
We also impose that the far-field pressure respects the potential solution:

$$
\overline{p}_\star(\xi, \eta)
= \mathrm{Re} \, \overline{p}(\xi, \eta)
= - \dfrac{\mathrm{Re}}{8 \pi^2 } \, \dfrac{1}{\xi^2 + \eta^2}
$$
*/

p[right] = dirichlet(- Rey*(1./(8.*pi*pi))/(x*x + y*y) );
p[left] = dirichlet(- Rey*(1./(8.*pi*pi))/(x*x + y*y) );
p[top] = dirichlet(- Rey*(1./(8.*pi*pi))/(x*x + y*y) );
p[bottom] = dirichlet(- Rey*(1./(8.*pi*pi))/(x*x + y*y) );

/** 
Here, the *Reynolds* number is arbitrary; for our needs, we will choose 
$\mathrm{Re} = 100$. 
*/



/** 
### Generic Events
*/

int main()
{

  Rey = 1.e2; // Reynolds number

  // viscosity
  const face vector mu_sim[] = {1., 1.};
  mu = mu_sim; 

  for (k = 0; k < 4; k++){    
    L0 = SIZE;
    #if !TREE
      init_grid(1 << lvl_grid[k]);
    #else
      init_grid(1 << LEVEL);
    #endif
    origin(-SIZE/2., -SIZE/2.);
    run() ;
  }
}


event init (t = 0){
  foreach() {
    double radius = sqrt(x*x + y*y) ;                                                                                                                                                                 
    u.x[] =  -1. * velocity_profile(radius) * y/radius ;
    u.y[] =  velocity_profile(radius) * x/radius ;
  }  
  #if TREE
    unrefine(sqrt(x*x + y*y) > 0.5*SIZE/2.);
    boundary({uf});
    boundary({lambdaf});
  #endif
  // Self-Similar solution for the η-component of the velocity field
  foreach()
    uy_sol[] = (1./(2.*pi)*(x/(x*x + y*y)))*(1. - exp(-(x*x + y*y)/4.)) ; 
}


event vorti (i++){
  vorticity (u, omega);
}

#if TREE
event adapt (i++){
  unrefine(sqrt(x*x + y*y) > 0.5*SIZE/2.);
  boundary({uf});
  boundary({lambdaf});
}
#endif

event end (t = T_END) {
  printf("\ni = %d t = %g\n", i, t);
}


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
}


/** 
### Outputs 

We store the vertical velocity component values for different times, according 
to the current level of discretization of the grid: 
*/



event profiles (t = {1, 2, 5, 10}){
  char filename[100];
  #if TREE
    sprintf(filename, "uy_exact_func_AMR_N%d_t%g.dat", LEVEL, t);
  #else
    sprintf(filename, "uy_exact_func_N%d_t%g.dat", lvl_grid[k], t);
  #endif
  FILE * fp = fopen(filename, "w");
  for (double x = -L0/2. ; x <= L0/2. ; x += 0.01)
      fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}


/** 
Then, we plot the results for the time and spatial convergence study (figure 
below):

~~~pythonplot Time and Spatial convergences
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

times_str = ["1", "2", "5", "10"] 
times_float = np.array([1, 2, 5, 10])

lvl_str = ["6", "7", "8", "9", "10"]
lvl_float = np.array([6, 7, 8, 9, 10])

Nmax = 9

path_uy_prof = "uy_exact_func_" 



#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 5
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]


xi_array = np.linspace(-20, 20, 1000)

def selfsim_sol_lamb(xi_array):
  a = 0.5/np.pi
  b = 1./xi_array
  c = 1. - np.exp(- np.square(xi_array)/4.)
  return (a*b*c)

f_y = selfsim_sol_lamb(xi_array)


# Velocity Profiles directly obtained in Self-Similar Space 

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

#ax.set_box_aspect(1) # allows a square box independent of data limits

#axins = zoomed_inset_axes(ax, 10, loc='upper right') # frame for zoom

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
#mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.5")   
#axins.set_facecolor('#fffcf4')

# try other functions for the Python Basilisk parser:
axins = inset_axes(ax, width=1, height=.7, borderpad=.1, loc=1, 
                   bbox_to_anchor=(.55, .0625, .45, .9125),
                   bbox_transform=ax.transAxes)
axins.patch.set_color('#fffcf4')


ax.set_xlim(-20,20)
ax.set_ylim(-0.06,0.06)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\overline{u}_\eta$')

axins.set_xlim([1.8, 2.7])
axins.set_ylim([0.049, 0.051])

for k, t in enumerate(times_str):
  uy_prof = np.loadtxt (path_uy_prof + "N{}_t{}.dat".format(Nmax, t))
  ax.plot(uy_prof[:, 0], uy_prof[:, 1], 
          lw=1.2, color = mapcolors[k], label=r'$\tau = {}$'.format(t))
  axins.plot(uy_prof[:, 0], uy_prof[:, 1], lw=1.2, color = mapcolors[k])

ax.plot(xi_array, f_y, '--', lw=0.8, color="darkred", label=r'Self-sim. sol.')
axins.plot(xi_array, f_y, '--', lw=0.8, color="darkred")
  
ax.legend(frameon=False, loc=2)
ax.text(10,-0.055, r'$N_{max} = 9$', fontsize = 15)



# Velocity Profiles for différent Δ obtained directly in Self-Similar Space 

ax2.set_xlim(-20,20)
ax2.set_ylim(-0.06,0.06)
ax2.set_xlabel(r'$\xi$')

t = 10

Nmax_str = r'$N_{max} \, \, $'

for k, Nmax in enumerate(lvl_str):
  if k == 4:
    uy_prof = np.loadtxt ("../lamb_data/" + path_uy_prof + "N{}_t{}.dat".format(Nmax, t))
  else:
    uy_prof = np.loadtxt (path_uy_prof + "N{}_t{}.dat".format(Nmax, t))
  ax2.plot(uy_prof[:, 0], uy_prof[:, 1], 
          lw=1.2, color = mapcolors[k], label= Nmax_str + r'$= {}$'.format(Nmax))

ax2.plot(xi_array, f_y, '--', lw=0.8, color="darkred", 
        label=r'Self-sim. sol.')
  
ax2.legend(frameon=False, loc=2)
ax2.text(13,-0.055, r'$\tau = 10$', fontsize = 15)

plt.savefig('lamb_selfsim_uy_exact_time+lvl_conv.svg') 
~~~

For a low number of discretized elements ($N_{max} < 8$, maximum $128$ elements), 
only the far-field vertical component $\overline{u}_\eta$ is correctly 
reproduced. However, profiles for levels $9$ and $10$ fit perfectly the 
analytical solution, showing the necessary resolution needed to capture all 
the physics inside the vortex core. 
The time convergence proves that the self-similar solver solution practically 
does not deviate from the steady state solution. A small deviation with the running 
time of the simulation is though observed when looking very closely (*cf. inset*), 
that we explain by *numerical diffusion*, due to successive projection solvings 
of a steady velocity field. A finer mesh tends to limit these slight discrepancies. 

Anyway, this excellent agreement with the theory **validates** the self-similar 
solver for single-phase fluids. 


It is also of utmost importance to ensure that the self-similar solver, due to 
its numerous changes in the numerical schemes, did not downgrade the time and 
spatial order of `Basilisk` (**second** order scheme). 
A classical test used in several codes on this website is to measure the 
[*uniform norm*](https://en.wikipedia.org/wiki/Uniform_norm) between the 
reference solution and the numerical one obtaineds with `Basilisk`, when 
variying the spatial resolution of the grid. 

Hence, we compute such a difference and collect some statistics upon it thanks 
to [utility functions](http://basilisk.fr/src/utils.h#simple-field-statistics) 
of `Basilisk`, for the different grid levels tested (reminder: $N_{max} = 10$ 
done locally):
*/

event stats (t = {1, 2, 5, 10}){
  foreach()
    uy_diff[] = uy_sol[] - u.y[];
  norm n = normf(uy_diff);
  char fileconv[200];
  #if TREE
    sprintf(fileconv, "diff_uy_conv_AMR_N%d.dat", LEVEL);
  #else
    sprintf(fileconv, "diff_uy_conv_N%d.dat", lvl_grid[k]);
  #endif
  static FILE * fpconv = fopen(fileconv, "a");
  fprintf (fpconv, "%g %i %g %g %g\n", t, lvl_grid[k], n.avg, n.rms, n.max);
} 



/** 
~~~pythonplot Showing 2nd-order spatial convergence
times_str = ["1", "2", "5", "10"] 
times_float = np.array([1, 2, 5, 10])
ntimes = len(times_float)

marker_str = ["*", "o", "^", "d"]

lvl_str = ["6", "7", "8", "9", "10"]
lvl_float = np.array([6, 7, 8, 9, 10])

path_diff_uy = "diff_uy_conv_" 

xi_array = np.linspace(1e1, 2e3)


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 4
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]


# Convergence Error between the Analytical Self-Similar solution of the 
# Lamb-Oseen vortex and the self-similar DNS starting from the analytical solution 

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e1,2e3)
ax.set_ylim(1e-5,1e-1)
ax.set_xlabel(r'$2^{N_{max}}$')
ax.set_ylabel(r'$\Vert \overline{u}_\eta^{\mathrm{sol}} - \overline{u}_\eta^{\mathrm{num}} \Vert_{\infty}$')


for i, N in enumerate(lvl_str):
  diff_uy = np.loadtxt (path_diff_uy + "N{}.dat".format(6))
  for k, t in enumerate(times_str):
    if i == 0:
      ax.scatter(np.power(2.,diff_uy[k + ntimes*i, 1]), 
              diff_uy[k + ntimes*i, 4], 
              s=20, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[k], 
              label=r'$\tau = {}$'.format(t)
              )
    elif i == len(lvl_str) - 1:
      diff_uy = np.loadtxt ("../lamb_data/" + path_diff_uy + "N{}.dat".format(10))
      ax.scatter(np.power(2.,diff_uy[k, 1]), 
              diff_uy[k, 4], 
              s=20, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[k], 
              )        
    else:
      ax.scatter(np.power(2.,diff_uy[k + ntimes*i, 1]), 
              diff_uy[k + ntimes*i, 4], 
              s=20, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[k], 
              )      
ax.plot(xi_array, 1*np.power(xi_array, -1),
        '--', lw=0.8, color="darkred", label=r'$1/N_{max}$')
ax.plot(xi_array, 120*np.power(xi_array, -2),
        '--', lw=0.8, color="darkblue", label=r'$1/N_{max}^2$')
  
ax.legend(frameon=False, loc=3)

plt.savefig('lamb_selfsim_conv_exact_norminf.svg')
~~~

The above figure shows that two different regimes exist, with a first-order error 
at low number of elements $(2^{N_{max}} = 64$ and $128)$, and a **second-order** 
error for finer meshes $(2^{N_{max}} \geqslant 512)$. Therefore, the second-order 
spatial convergence is maintained! *Good job, self-similar solver!* 

Notice also the small error-increase in time due to the *numerical diffusion* 
discussed earlier.
*/












/**
## Additional Remarks

  + For overcoming unexpected behaviours when using mesh refinement, 
the following patch can be used: unrefine all borders to avoid propagation 
of undesirable boundary layers from corners.

  + For convergence studies, ALWAYS use the multigrid solver, as the quadtree 
path creates (for unknown reasons) a deviation of the solution around $\tau = 5$ 
with the exact analytical solution, due to a convergence of vorticity layers 
coming rapidly from the borders (can be seen when enabling the cells in BView).


# APPENDIX: analytical solutions in Cartesian coordinates

It has been seen for the 
[*Lamb--Oseen* problem simulated in the physical space](http://basilisk.fr/sandbox/cailler/lamb_oseen/lamb.c#towards-self-similarity) 
that the exact solution of the azimuth velocity is:

$$
u_\theta(r,t) 
= \dfrac{\Gamma}{2 \pi r} \left(1 - \mathrm{e}^{-r^2 / 4\nu t} \right )  
\quad (*)
$$

with $\Gamma$ the flow circulation. The Cartesian components related to the 
azimuth velocity are simply:

$$
u := u_x = - u_\theta \sin(\theta) = - u_\theta \dfrac{y}{\sqrt{x^2 + y^2}} 
\quad ; \quad 
v := u_y = u_\theta \cos(\theta) = u_\theta \dfrac{x}{\sqrt{x^2 + y^2}}  
\quad (**)
$$

Now, searching for scale invariant solutions to the *Navier--Stokes* 
equations describing the vortex (along with the initial and boundary conditions) 
one can exhibit the following self-similar variables:

$$
\xi(x,t) = \dfrac{x}{\sqrt{\nu t}} 
\quad ; \quad 
\eta(y,t) = \dfrac{y}{\sqrt{\nu t}} 
\quad (***)
$$

Using $(*)$, $(**)$, and $(***)$, we can finally derive the Cartesian 
self-similar analytical expressions of the velocity and pressure fields:

$$
\left\{\begin{array}{rcl}
u(x,y,t) & = & \dfrac{\Gamma}{\sqrt{\nu t}} \overline{u} \left[ \xi, \eta 
\right] = \dfrac{\Gamma}{\sqrt{\nu t}}
\left[
-\dfrac{1}{2 \pi} \dfrac{\eta}{\xi^2 + \eta^2}  
\left(1 - \mathrm{e}^{-\dfrac{\xi^2 + \eta^2}{4}} \right)
\right] \\ 
& & \\
v(x,y,t) & = & \dfrac{\Gamma}{\sqrt{\nu t}} \overline{v} \left[ \xi, \eta 
\right] = \dfrac{\Gamma}{\sqrt{\nu t}}
\left[
\dfrac{1}{2 \pi} \dfrac{\xi}{\xi^2 + \eta^2}  
\left(1 - \mathrm{e}^{-\dfrac{\xi^2 + \eta^2}{4}} \right)
\right] \\ 
& & \\
\dfrac{p(x,y,t)}{\rho} & = & \dfrac{\Gamma^2}{\nu t} \overline{p} \left[ \xi, 
\eta\right] \\
& & \\
& = & \dfrac{\Gamma ^2}{16 \pi^2 \text{$\nu $t}}
\left[
\text{Ei}\left(-\dfrac{\xi^2 + \eta^2}{4}\right)
- \text{Ei}\left(-\dfrac{\xi^2 + \eta^2}{2}\right)
- \dfrac{2}{\xi^2 + \eta^2}
  \left(1-\mathrm{e}^{-\dfrac{\xi^2 + \eta^2}{4}}\right)^2
\right]
\end{array}\right.  
$$

The analytical expression for the pressure field was derived thanks to 
the [German Wiki page](https://de.wikipedia.org/wiki/Hamel-Oseenscher-Wirbel#Druck) 
on the *Lamb--Oseen* vortex, where $\mathrm{Ei}$ is the *exponential integral* 
function.

*/