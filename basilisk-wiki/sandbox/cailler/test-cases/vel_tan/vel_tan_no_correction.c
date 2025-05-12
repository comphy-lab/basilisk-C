/** 
# The problem of advecting a 2D-liquid interface in a tangent velocity field

## Brief Report

### Description

Simple test of a straight gas-liquid interface crossing left and 
bottom borders in a tangent flow of intensity $U_0$, 
with surface tension taken into account.
No viscosity considered.
Code done in a "naive" way, *i.e.* only with the default solvers 
/ without any custom modifications, to demonstrate its inoperability.

### Expectation

In a tangent flow, the initial straight interface shall remain the same, 
especially at the contact angles.

### Results with Basilisk

![Vorticity field and gas-liquid interface in Multigrid](vel_tan_no_correction/movie_MG_v2.mp4)(width="600" height="600")


The interface is unexpectedly translating along the contact borders 
instead of staying steady, and is deforming at contact angles.

Of all the tests done, the following criteria were NOT the cause of the 
problem:

* velocity direction is correctly set and remains the same at all 
  borders during the DNS;
* AMR --> same pb in MG, but AMR generates additional capillary deformations;
* density ratio $\rho_1/\rho_2$ --> same pb with $\rho_1 = \rho_2$;
* Boundary Conditions:
    + with/without inflow & outflow conditions;
    + with/without pressure reference at top gas border;
* velocity intensity $U_0$ --> the larger the value, the more the interface 
  is deformed, but deformation + translation still appear when $U_0$ 
  is set to unity;

### Conclusion

As soon as initialization and step $i = 1$, bottom curvature 
at contact angle is $\mathcal{O}(1)$ whereas the rest of the interface 
is at $10^{-7}$.
Curvature computed by height functions is therefore not well 
calculated and propagates errors along the simulation, creating vorticity 
pockets moving the interface.
It means that volume fractions in ghost cells are not evaluated correctly.
Ways of improvement shall investigate their update for overcoming 
this bad behavior.
*/


/** 
## Code

### General Parameters
*/

#define LEVEL 7
#define MAXLEVEL 9
#define SIZE 1.e2
#define T_END 100. 
#define THETA_0 120.*pi/180.
#define BETA_0 pi - THETA_0
#define Y_OFFSET 0. 
#define X_OFFSET 0.55*SIZE 
#define U0 1e0

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"

// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 




/** 
###  Boundary Conditions

The simulation takes place in a box of arbitrary dimensions 
$[0 \, ; 100] \times [0 \, ; 100]$, where the liquid wedge 
in the bottom-left corner is described by the slope: 
$$
y = (x_{off} - x) \, \tan \beta_0
$$
where $x_{off}$ is an offset on the bottom axis 
for which $y = 0$.

Outflow conditions are considered on the top side, 
inflow on the bottom, and other boundaries are defined 
by the analytical tangent flow expression of intensity $U_0$.

![Initialization of the interface](img_vel_tan/vel_tan_no_correc_init.svg)

~~~pythonplot Initialization of the interface
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

# https://stackoverflow.com/a/11158224
# importing a module from a parent directory
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
import angle_helper_func as ang

# Use LaTeX fonts in the plot:
plt.rc('text', usetex=True)
plt.rc('font', family='serif') 
plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
plt.rc('text.latex', preamble=r'\usepackage{mlmodern}')
mpl.rcParams['figure.constrained_layout.use'] = True
label_size = 15
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 



#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
path_shape = "vel_tan_no_correc_shape_"
path_vp = "vel_tan_no_correc_u_p_" # x, y, u.x, u.y, p, omega

facets_init = np.loadtxt (path_shape + "t0.dat", unpack=False)
x, y, ux, uy, p, omega =  np.loadtxt(path_vp + "t0.dat", unpack=True)

# Get segments of the output_facets() function of Basilisk C:
N_seg_init = int (0.5*facets_init.shape[0])
segments_init = np.split (facets_init, indices_or_sections=N_seg_init)



#                            DATA INTERPOLATION
# ------------------------------------------------------------------------------
# interpolate 1D data to a 2D grid (for streamlines plot)
# https://stackoverflow.com/questions/35877478/matplotlib-using-1-d-arrays-in-streamplot

x2 = np.linspace(x.min(), x.max(), 1000)
y2 = np.linspace(y.min(), y.max(), 1000)
X, Y = np.meshgrid(x2, y2)
OMEGA = interpolate.griddata((x, y), omega, (X, Y), method='cubic')
UX = interpolate.griddata((x, y), ux, (X, Y), method='cubic')
UY = interpolate.griddata((x, y), uy, (X, Y), method='cubic')

#                                   PLOTS
# ------------------------------------------------------------------------------
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,100)
ax.set_ylim(0,100)
# ax2 = ax.twinx()
# ax2.set_ylim(0,100)
ax.set_title(r'$\partial_n u_n = 0, \,\, p = 0$')
ax.set_xlabel(r'$u_n = U_0 \sin(\beta_0), \,\, u_t = - U_0 \cos(\beta_0),$' 
              + r'$\,\, \partial_n p = 0$')
ax.set_ylabel(r'$u_n = - U_0 \cos(\beta_0), \,\, u_t = U_0 \sin(\beta_0)$')
secax_y = ax.secondary_yaxis('right')
secax_y.set_ylabel(r'$u_n = - U_0 \cos(\beta_0), \,\, u_t = U_0 \sin(\beta_0)$')


# Draw the actual figure:
# -----------------------
for segment in segments_init:
    ax.fill_between(segment[:, 0], segment[:, 1], color = "#e9effb")
    ax.plot(segment[:, 0], segment[:, 1], '-', lw=1.5, color = "#5383dc") 

skip = (slice(None, None, 50), slice(None, None, 50))
ax.quiver(X[skip], Y[skip], UX[skip], UY[skip], color='black', width=0.002)

# Draw initial wedge angle:
# -------------------------
center = (55., 0.)
p1 = [(0., 95.2628), (55., 0.)]
p2 = [(0., 0.), (55., 0.)]

angle = ang.AngleAnnotation(center, p1[0], p2[0], ax=ax, 
                      size=75, text=r"$\beta_0$", color = "#af3127", 
                      text_kw=dict(fontsize=12, color="#af3127"))

plt.savefig('vel_tan_no_correc_init.svg', 
  bbox_inches='tight',
  pad_inches = 0.1) 
~~~
*/



// Velocity BC
u.n[left] = dirichlet ( - U0 * cos(BETA_0) ) ;
u.t[left] = dirichlet ( U0 * sin(BETA_0) ) ;
u.n[right] = dirichlet ( - U0 * cos(BETA_0) ) ;
u.t[right] = dirichlet ( U0 * sin(BETA_0) ) ;

// Face Velocity BC
uf.n[left] =  - U0 * cos(BETA_0)  ;
uf.t[left] =  U0 * sin(BETA_0)  ;

// Outflow Condition
u.n[top] = neumann(0.);
p[top] = dirichlet (0.);
pf[top] = dirichlet (0.);

// Inflow Condition
u.t[bottom] = dirichlet ( - U0 * cos(BETA_0) ) ;
u.n[bottom] = dirichlet ( U0 * sin(BETA_0) ) ;
uf.t[bottom] =  - U0 * cos(BETA_0)  ;
uf.n[bottom] =  U0 * sin(BETA_0)  ;
p[bottom]    = neumann(0.);
pf[bottom]   = neumann(0.);

// Contact Angles to maintain 
h.t[left] = contact_angle (THETA_0 - pi/2.) ; 
h.t[bottom] = contact_angle (BETA_0) ;

/** 
###  Generic Events
*/

int main() {
  size(SIZE) ;
  init_grid(1 << LEVEL) ;
  origin(0, 0) ;

  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 0, mu2 = 0; // disable viscosity

  run() ;
}

event init (t = 0){
  // Initial free-surface definition:

  #if TREE
    do{
    fraction(f, ((- (x - X_OFFSET) * tan(BETA_0) + Y_OFFSET) - y) ) ; 
    } while (adapt_wavelet ({f}, (double[]){1e-3}, MAXLEVEL).nf); 
  #else
    fraction(f, ((- (x - X_OFFSET) * tan(BETA_0) + Y_OFFSET) - y) ) ;
  #endif

  foreach(){
    u.x[] = - U0 * cos(BETA_0) ; 
    u.y[] = U0 * sin(BETA_0);
  }

  // curvature computation:
  curvature (f, curv_viz, 1, add = false);
}

event vorti (i++){
  vorticity (u, omega);
  curvature (f, curv_viz, 1, add = false);
}

event adapt (i++) {
  #if TREE
    adapt_wavelet ({f,u}, (double[]){1e-3,1e-3,1e-3}, MAXLEVEL, LEVEL-1);
  #endif

  // ensuring that height functions are up-to-date:
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
###  Outputs
*/

event profiles (t = {0, 50, 100}){                 
  char filename[200] ;
  sprintf(filename, "vel_tan_no_correc_shape_t%g.dat", t) ;
  FILE * fp = fopen(filename, "w") ;
  output_facets (f, fp) ;
  fclose (fp) ;
}


event vel_p_maps (t = {0, 50, 100})
{
  scalar x_interf[];
  scalar y_interf[];
  char fileup[200];
  char fileinterf[200];
  sprintf (fileup, "vel_tan_no_correc_u_p_t%g.dat", t);
  sprintf (fileinterf, "vel_tan_no_correc_interf_t%g.dat", t);
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
Results on a multigrid simulation show important 
vorticity layers generated from contact angles, 
compared to the expected steady position of the interface 
(*dashed line*):

![Vorticity at $t = 100$](img_vel_tan/vel_tan_no_correc_t100.svg)

~~~pythonplot Vorticity at $t = 100$
t = 100

#                             DATA EXTRACTION
# ------------------------------------------------------------------------------
path_interf = "vel_tan_no_correc_interf_" # x_interf, y_interf

facets = np.loadtxt (path_shape + "t{}.dat".format(t), unpack=False)
xi, yi = np.loadtxt(path_interf + "t{}.dat".format(t), unpack = True)
x, y, ux, uy, p, omega =  np.loadtxt(path_vp + "t{}.dat".format(t), unpack=True)

# Get segments of the output_facets() function of Basilisk C:
N_seg = int (0.5*facets.shape[0])
segments = np.split (facets, indices_or_sections=N_seg)

#                            DATA INTERPOLATION
# ------------------------------------------------------------------------------
# interpolate 1D data to a 2D grid (for streamlines plot)
# https://stackoverflow.com/questions/35877478/matplotlib-using-1-d-arrays-in-streamplot

x2 = np.linspace(x.min(), x.max(), 1000)
y2 = np.linspace(y.min(), y.max(), 1000)
X, Y = np.meshgrid(x2, y2)
OMEGA = interpolate.griddata((x, y), omega, (X, Y), method='cubic')
UX = interpolate.griddata((x, y), ux, (X, Y), method='cubic')
UY = interpolate.griddata((x, y), uy, (X, Y), method='cubic')

#                                   PLOTS
# ------------------------------------------------------------------------------
fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

plt.imshow(OMEGA, extent=[0, 100, 0, 100], origin='lower',
          norm=mcolors.CenteredNorm(halfrange=0.1), cmap="seismic", alpha=1) 
cbar = plt.colorbar(aspect=30, pad=0.02)
cbar.ax.set_title(r'$\omega$', fontsize=18)

skip = (slice(None, None, 50), slice(None, None, 50))
ax.quiver(X[skip], Y[skip], UX[skip], UY[skip], color='black', width=0.002)

for segment in segments:
    ax.plot(segment[:, 0], segment[:, 1], '-k', lw=1.2) 
for segment in segments_init[::2]:
    ax.plot(segment[:, 0], segment[:, 1], '-k', lw=1.2) 

plt.savefig('vel_tan_no_correc_t{}.svg'.format(t), 
  bbox_inches='tight',
  pad_inches = 0.1) 
~~~

### Curvature evolution 

This unexpected behaviour finds its origin in a 
**wrong computation** of the *curvature* at contact angles, 
when surface tension is taken into account. 
Indeed, when we collect curvature values on the bottom axis, 
it rapidly converges towards a $5.10^{-1}$ error:
*/




event bottom_curv (t += 1){
  char filecurv[200] ;
  #if TREE
    sprintf(filecurv, "vel_tan_no_correc_AMR_curv.dat") ;
  #else
    sprintf(filecurv, "vel_tan_no_correc_curv.dat") ;
  #endif
  static FILE * fcurv = fopen(filecurv, "a") ;

  double xpos = 0., curv = 0.;
  foreach_boundary (bottom){
    if (interfacial (point,f)){ // if we are at the interface bottom
      if ( (f[] > 0.1) && (f[] < 0.9) ){
        xpos = x;
        curv = curv_viz[]; 
      }
    }
  }
  fprintf(fcurv, "%g %g %g\n", t, xpos, curv);
  fflush(fcurv); // empty the buffer
}


/**
![Curvature value at the bottom axis](img_vel_tan/vel_tan_curv_no_correc.svg)

~~~pythonplot Curvature value at the bottom axis
t_noc, x_noc, curv_noc = np.loadtxt("vel_tan_no_correc_curv.dat", unpack=True)
t_noc_filt = t_noc[np.abs(curv_noc) > 0]
curv_noc_filt = curv_noc[np.abs(curv_noc) > 0]

Ncolors = 20
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

fig, ax = plt.subplots()
ax.set_yscale('log')

ax.set_xlim(0,100)
ax.set_ylim(1e-3,1)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$|\kappa_{bottom}|$')

ax.plot(t_noc_filt, np.abs(curv_noc_filt), lw=1.2, label=r'No corrections', 
        color=mapcolors[4])

ax.legend(frameon=False, loc=4)
plt.savefig('vel_tan_curv_no_correc.svg', 
            bbox_inches='tight', pad_inches = 0.1) 
~~~
*/


#if TREE
event print_mesh(t=T_END){
  FILE*fpmsh=fopen("vel_tan_no_correc_AMR_mesh","w");
  foreach(){
    if (f[]){
      fprintf(fpmsh,"%g %g\n", x-Delta/2., y-Delta/2.);
      fprintf(fpmsh,"%g %g\n", x-Delta/2., y+Delta/2.);
      fprintf(fpmsh,"%g %g\n", x+Delta/2., y+Delta/2.);
      fprintf(fpmsh,"%g %g\n", x+Delta/2., y-Delta/2.);
      fprintf(fpmsh,"%g %g\n\n", x-Delta/2., y-Delta/2.);
    }
  }
  fclose(fpmsh);
}
#endif


/**
### AMR distortions

If Adaptive Mesh Refinement (AMR) is enabled, then we 
do have additionally strong capillary fluctuations that 
are generated during the simulation and coming from 
contact angles:

![AMR distortions](img_vel_tan/vel_tan_no_correc_t100_closeup_AMR.svg)

~~~pythonplot Adaptive Mesh Refinement Plot (Zoom)
path_shape_AMR = "vel_tan_no_correc_AMR_shape_"
path_interf_AMR = "vel_tan_no_correc_AMR_interf_" # x_interf, y_interf
path_vp_AMR = "vel_tan_no_correc_AMR_u_p_" # x, y, u.x, u.y, p, omega

facets_AMR = np.loadtxt (path_shape_AMR + "t{}.dat".format(t), unpack=False)
facets_init_AMR = np.loadtxt (path_shape_AMR + "t0.dat", unpack=False)
facets_msh = np.loadtxt ("vel_tan_no_correc_AMR_mesh.dat", unpack=False)

# Get segments of the output_facets() function of Basilisk C:
N_seg_AMR = int (0.5*facets_AMR.shape[0])
segments_AMR = np.split (facets_AMR, indices_or_sections=N_seg_AMR)
N_seg_init_AMR = int (0.5*facets_init_AMR.shape[0])
segments_init_AMR = np.split (facets_init_AMR, indices_or_sections=N_seg_init_AMR)
N_seg_msh = int (0.2*facets_msh.shape[0])
segments_msh = np.split (facets_msh, indices_or_sections=N_seg_msh)


fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.set_xlim(0,20)
ax.set_ylim(70,100)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

for segment in segments_msh:
    ax.plot(segment[:, 0], segment[:, 1], lw=.5, color="darkgray") 
for segment in segments_AMR:
    ax.plot(segment[:, 0], segment[:, 1], '-k', lw=1.2) 
for segment in segments_init[::2]:
    ax.plot(segment[:, 0], segment[:, 1], '-k', lw=1.2) 

plt.savefig('vel_tan_no_correc_t{}_closeup_AMR.svg'.format(t), 
  bbox_inches='tight',
  pad_inches = 0.1) 
~~~
*/

/** 
## Movie
*/

#include "view.h"

event movie (t += 1; t <= T_END) {
  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  clear();
  view (tx = -0.5, ty = -0.5);
  box();
  draw_vof ("f", lw = 1.5);
  squares ("omega", linear = true, 
           max = 0.05, min = -0.05,
           map = blue_white_red);
  draw_string(legend, 1, size = 20., lw = 2.1);
  save ("movie_MG_v2.mp4");
}