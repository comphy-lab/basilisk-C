/**
# The *Lamb--Oseen* self-similar solver acts like a *blackhole* 

![The "*blackhole*"](selfsim_lamb_blackhole/movie_selfsim_blackhole_N8.mp4)(width="500" height="500")


> On the movie is well illustrated the crucial role of the additional non-linear 
advecting term 
$\textcolor{Orchid}{\left[
  \boldsymbol{\nabla} \cdot 
    \left(\mathbf{\overline{\Lambda}} \otimes \overline{\mathbf{u}} \right) 
  \right]^{n+1/2}}$ in the numerical scheme of the self-similar centered *Navier--Stokes* solver. 
It acts literally as the "right speed" to move with the physical phenomenon 
in order to observe continuously the same state, *i.e.* to be *steady*. 

> The absence of any *diffusion* in the solid core of the vortex once the 
scale invariant solution is obtained results from the source term 
$\textcolor{green}{-\dfrac{1}{2}  \, \mathbf{\overline{u}}^{n+1/2}}$: 
it hinders the vorticity layer spread in the fluid domain.


## Motivations

This is a numerical experiment extending drastically the results of the 
[first test of convergence](http://basilisk.fr/sandbox/cailler/lamb_oseen/selfsim_lamb_conv.c) 
of the self-similar solver dedicated to the *Lamb--Oseen* vortex. 
Starting from a **random** initial velocity field associated to a 
[*white noise*](http://basilisk.fr/src/common.h), the question at stake is: 
do we inexorably converge towards a *steady state*? If so, is this the exact 
same self-similar solution than the one analytically predicted? 
We are hence interested in the long-time behaviour of the self-similar solver. 

Aside from the initial condition, the problem's nature is the same as before, 
therefore, we do not change any other parameter (grid size, boundary conditions, 
*Reynolds* number...), so **please** refer to the 
[previous test case](http://basilisk.fr/sandbox/cailler/lamb_oseen/selfsim_lamb_conv.c) 
for any details upon the choices made and other technical/theoretical aspects.
*/


/** 
## Code

### General Parameters 

The total running time of the simulation is $\tau = 100$ to let enough time 
for the system to evolve. 
*/


#include "grid/multigrid.h"
#include "selfsim_centered_lamb.h"
#include <time.h> // for setting seeds (initializers in pseudorandom number generators)

#define LEVEL 8 // minimum grid size for AMR
#define MAX_LEVEL 9 // maximum refinement level for AMR
#define SIZE 200
#define T_END 100

scalar uy_sol[];
scalar uy_diff[];
scalar omega[];

/** 
Several discretization levels $N_{max}$ are tested:

$$
\forall \, N_{max} \in \{6 \, ; 10 \}
\,\, \Rightarrow \,\, 
\Delta = L_0 / 2^{N_{max}} \in \{0.195 \, ; 3.125 \}
$$

However, the cases $N_{max} = 9$ and $10$ are done remotely, as it would exceed the 
`Basilisk` server capacities due to the total running time of the simulation needed. 
Related results can be found in 
[`lamb_data`](http://basilisk.fr/sandbox/cailler/lamb_oseen/lamb_data/). 
*/

int k; // counter
int lvl_grid[] = {6, 7, 8};

/**
### Boundary conditions

The far-field velocity respects the potential solution:
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
So does the far-field pressure:

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
### Generic Events
*/

int main()
{
  Rey = 1.e2; // Reynolds number

  // viscosity
  const face vector mu_sim[] = {1., 1.};
  mu = mu_sim; 

  for (k = 0; k < 3; k++){    
    /** 
    We first initialize a 
    [*random number generator*](https://www.quora.com/Why-do-we-use-srand-time-NULL) 
    in order to be sure that the initial chaotic state of the velocity field 
    is not the same for each run, thus, showing the **true universality** of 
    the simulation's final result:
    */
    srand((unsigned) time(NULL));
    /** 
    
    */
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
  /** 
  We use a *white noise* as an initial guess for the velocity field, 
  that is to say, the most random field possible to be the furthest from 
  the analytical solution:
  */
  foreach() {
    u.x[] = noise() ;
    u.y[] = noise() ;
  }
  /** 
   If using Adaptive Mesh Refinement, we 
   unrefine all borders to avoid propagation of undesirable 
   boundary layers of vorticity coming from corners:
   */
 #if TREE
    unrefine(x < - .8*SIZE/2.);
    unrefine(y < - .8*SIZE/2.);
    unrefine(x > .8*SIZE/2.);
    unrefine(y > .8*SIZE/2.);

    boundary({uf});
    boundary({lambdaf});
  #endif
  /** 
  We store the self-similar analytical solution for the $\eta$--component 
  of the velocity field:
  */
  foreach()
    uy_sol[] = (1./(2.*pi)*(x/(x*x + y*y)))*(1. - exp(-(x*x + y*y)/4.)) ; 
}

/** 

*/

event vorti (i++){
  vorticity (u, omega);
}


#if TREE
event adapt (i++){
  adapt_wavelet ((scalar *){u}, (double []){1e-3, 1e-3}, MAX_LEVEL, LEVEL-1);
  unrefine(x < - .8*SIZE/2.);
  unrefine(y < - .8*SIZE/2.);
  unrefine(x > .8*SIZE/2.);
  unrefine(y > .8*SIZE/2.);
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

We store the vertical velocity component values for different times, here, for 
the maximum level of discretization of the grid running on the `Basilisk` server, 
so $N_{max} = 8$: 
*/

event profiles (t = {1, 2, 5, 10, 20, 50, 100}){
  if (lvl_grid[k] == 8){
    char filename[100];
    #if TREE
      sprintf(filename, "uy_noise_func_AMR_N%d_t%g.dat", LEVEL, t);
    #else
      sprintf(filename, "uy_noise_func_N%d_t%g.dat", lvl_grid[k], t);
    #endif
    FILE * fp = fopen(filename, "w");
    for (double x = -L0/2. ; x <= L0/2. ; x += 0.01)
        fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
    fclose (fp);
  }
}


/**
~~~pythonplot Time convergence for a *white noise* as initial condition of the velocity 
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

times_str = ["1", "2", "5", "10", "20", "50", "100"] 
times_float = np.array([1, 2, 5, 10, 20, 50, 100])

Nmax = 8

path_uy_prof = "uy_noise_func_" 



#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 7
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


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

ax.set_xlim(-20,20)
ax.set_ylim(-0.06,0.06)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\overline{u}_\eta$')

ax2.set_xlim(-20,20)
ax2.set_ylim(-0.06,0.06)
ax2.set_xlabel(r'$\xi$')

for k, t in enumerate(times_str):
  uy_prof = np.loadtxt (path_uy_prof + "N{}_t{}.dat".format(Nmax, t))
  if k < 3: # Transitory state in Self-Similar Space 
    ax.plot(uy_prof[:, 0], uy_prof[:, 1], 
            lw=1.2, color = mapcolors[k], label=r'$\tau = {}$'.format(t))
  else: # Towards the steady state in Self-Similar Space 
    ax2.plot(uy_prof[:, 0], uy_prof[:, 1], 
            lw=1.2, color = mapcolors[k], label=r'$\tau = {}$'.format(t)) 

ax.plot(xi_array, f_y, '--', lw=0.8, color="darkred")
ax2.plot(xi_array, f_y, '--', lw=0.8, color="darkred", 
        label=r'Self-sim. sol.')
  
ax.legend(frameon=False, loc=2)
ax.text(10,-0.055, r'$N_{max} = 8$', fontsize = 15)

ax2.legend(frameon=False, loc=2)
ax2.text(10,-0.055, r'$N_{max} = 8$', fontsize = 15)

plt.savefig('lamb_selfsim_uy_noise_time_conv.svg') 
~~~ 

On the above figure, the velocity profile $\overline{u}_\eta$ evolves from 
a *random* distribution to the steady self-similar solution around $\tau \sim 50$ 
for this numerical experiment and for $N_{max} = 8$ (a shorter time of convergence 
can be reached for a finner mesh grid and/or depending on the intial seed 
impacting the unsteady regime). 

*/






/** 
We register some valuable informations upon various fields of interest:
*/

event vel_p_maps (t = {1, 2, 5, 10, 20, 50, 100})
{
  if (lvl_grid[k] == 8){
    char fileup[200];
    sprintf (fileup, "lamb_noise_u_p_N%d_t%g.dat", lvl_grid[k], t);
    FILE * fup = fopen (fileup, "w");
    foreach ()
      fprintf(fup, "%g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[], omega[]);
    fclose(fup);
  }
}

/** 
Finally, we verify the convergence order of the numerical scheme used in the 
self-similar solver, to ensure that we are recovering the conclusions of the 
[first test case](http://basilisk.fr/sandbox/cailler/lamb_oseen/selfsim_lamb_conv.c). 
We measure for various numbers of grid elements (between $64$ and $1024$) the 
Euclidian and Uniform norms of the differences between our numerical results 
and the analytical solution:
*/

event stats (t = T_END){
  foreach()
    uy_diff[] = uy_sol[] - u.y[];
  norm n = normf(uy_diff);
  char fileconv[200];
  #if TREE
    sprintf(fileconv, "diff_uy_noise_conv_AMR_N%d.dat", LEVEL);
  #else
    sprintf(fileconv, "diff_uy_noise_conv_N%d.dat", lvl_grid[k]);
  #endif
  static FILE * fpconv = fopen(fileconv, "a");
  fprintf (fpconv, "%g %i %g %g %g\n", t, lvl_grid[k], n.avg, n.rms, n.max);
} 

/**
~~~pythonplot Showing 2nd-order spatial convergence
times_str = ["100"] 
times_float = np.array([100])

marker_str = ["d"]

lvl_str = ["6", "7", "8", "9", "10"]
lvl_float = np.array([6, 7, 8, 9, 10])

path_diff_uy = "diff_uy_noise_conv_" 
path_diff_uy_data = "../lamb_data/diff_uy_noise_conv_" 

xi_array = np.linspace(1e1, 2e3) 


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 7.5))

ax.set_xscale('log')
ax.set_yscale('log')

ax2.set_xscale('log')
ax2.set_yscale('log')

ax.set_xlim(1e1,2e3)
ax.set_ylim(1e-6,1e0)
ax.set_xlabel(r'$2^{N_{max}}$')
ax.set_ylabel(r'$\Vert \overline{u}_\eta^{\mathrm{sol}} - \overline{u}_\eta^{\mathrm{num}} \Vert_{2}$')

ax2.set_xlim(1e1,2e3)
ax2.set_ylim(1e-5,1e1)
ax2.set_xlabel(r'$2^{N_{max}}$')
ax2.set_ylabel(r'$\Vert \overline{u}_\eta^{\mathrm{sol}} - \overline{u}_\eta^{\mathrm{num}} \Vert_{\infty}$')


for i, N in enumerate(lvl_str):
  diff_uy = np.loadtxt (path_diff_uy + "N{}.dat".format(6))
  for k, t in enumerate(times_str):
    if i == 0:
      ax.scatter(np.power(2.,diff_uy[i, 1]), diff_uy[i, 3], 
              s=30, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[5], 
              label=r'$\tau = {}$'.format(t)
              )
      ax2.scatter(np.power(2.,diff_uy[i, 1]), diff_uy[i, 4], 
              s=30, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[5], 
              label=r'$\tau = {}$'.format(t)
              )
    elif i > len(lvl_str) - 3:  
      diff_uy = np.loadtxt (path_diff_uy_data + "N9-10.dat")  
      ax.scatter(np.power(2.,diff_uy[i-3, 1]), diff_uy[i-3, 3], 
        s=30, # marker size 
        marker = marker_str[k],
        facecolors = 'none',
        edgecolors = mapcolors[5], 
        )    
      ax2.scatter(np.power(2.,diff_uy[i-3, 1]), diff_uy[i-3, 4], 
        s=30, # marker size 
        marker = marker_str[k],
        facecolors = 'none',
        edgecolors = mapcolors[5], 
        )    
    else:
      ax.scatter(np.power(2.,diff_uy[i, 1]), diff_uy[i, 3], 
              s=30, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[5], 
              )      
      ax2.scatter(np.power(2.,diff_uy[i, 1]), diff_uy[i, 4], 
              s=30, # marker size 
              marker = marker_str[k],
              facecolors = 'none',
              edgecolors = mapcolors[5], 
              )   
ax.plot(xi_array, 1.5*np.power(xi_array, -1),
        '--', lw=0.8, color="darkred", label=r'$1/N_{max}$')
ax.plot(xi_array, 3*np.power(xi_array, -2),
        '--', lw=0.8, color="darkblue", label=r'$1/N_{max}^2$')
  
ax2.plot(xi_array, 1.5*np.power(xi_array, -1),
        '--', lw=0.8, color="darkred", label=r'$1/N_{max}$')
ax2.plot(xi_array, 170*np.power(xi_array, -2),
        '--', lw=0.8, color="darkblue", label=r'$1/N_{max}^2$')


ax.legend(frameon=False, loc=3)
ax2.legend(frameon=False, loc=3)

plt.savefig('lamb_selfsim_conv_noise_norm2+inf.svg')
~~~

From the above results, it is very clear now that the **self-similar solver 
verifies a second-order scheme**, since for long simulations $(\tau = 100)$ and 
$N_{max} \geqslant 7$, all errors for both norms collapse on the line 
$y = 1/N_{max}^2$ in *log--log* scale. 
*/


/** 
### Movie

We make a movie of this beautiful numerical experiment acting like a blackhole:
*/

#include "view.h"

event viewing (t += 0.1) {
  if (lvl_grid[k] == 8){
  view (width = 480, height = 480, fov = 18, tx = 1e-4, ty = 1e-4,
	bg = {.7, .7, .7});

  clear();
  stats s = statsf (omega);
  double max_range = max(s.max, -s.min);
  squares ("omega", min = -max_range, max = max_range, 
    linear = true, map = blue_white_red); 
  box (notics = true);
  char legend[100];
  sprintf (legend, "ln t = %0.2g", t);
  draw_string (legend, 1, size = 20., lw = 1.7); // “0” bl, “1” tl, “2” tr and “3” br 
  save ("movie_selfsim_blackhole_N8.mp4");
  }
}


