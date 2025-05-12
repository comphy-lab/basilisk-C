/** 
# General Overview 

## Description

It has been shown in this ["naive" test case](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_no_correction.c) 
that a straight gas-liquid interface placed in a tangent flow 
and crossing left and bottom borders of the simulation box 
does not behave as intended when surface tension is taken into account 
and by using `Basilisk` ready-to-use functions. 

Indeed, the interface described by the slope: 
$$
y = (x_{off} - x) \, \tan \beta_0
$$
where $x_{off}$ is an offset on the bottom axis, 
shall remain steady in such a flow of intensity $U_0$, 
described in the figure below:

![Initialization of the interface](img_vel_tan/vel_tan_no_correc_init.svg)

Therefore, some modifications have been made to the initial code given 
[here](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_no_correction.c) 
to tackle *curvature computation* problems arising from contact angles 
subjected to inflow/outflow:

  + automatic update of ghost cells volumic fractions via 
    [f_BC_vel_tan.h](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/f_BC_vel_tan.h); 
  
  + introduction of an *auxiliary volume fraction field* `f0[]` 
    applied as a *Dirichlet* condition to the advected volume 
    fraction field `f[]`, in order to store the steady state 
    expected in ghost cells;

  + correction of the normal at the interface at boundaries, 
    as symmetry conditions are imposed in the default `vof.h` file: 
    some tweaks in `vof.h` were needed to pass the proper BCs, 
    see [vof_vel_tan.h](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vof_vel_tan.h);

  + systematic reevaluation of `f[]` boundaries and height functions after 
    each `foreach()` procedure they are used in (it is extremely important for 
    AMR use to do so, otherwise the aforementioned corrections are not applied 
    properly). 

## Results

![Vorticity field and gas-liquid interface with AMR](vel_tan_corrected/movie_AMR_v2.mp4)(width="600" height="600")

The colormap scale for the vorticity used here is 50 times lower than 
the one used in the *naive* code, otherwise one could think that nothing 
happens during the simulation 
(and this is precisely our goal: *nothing should happen*.) 

Contrary to the [*naive* test case](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_no_correction.c), 
Adaptive Mesh Refinement (AMR) does not produce capillary perturbations of 
the interface which remains straight as expected.  

# Code 

## General Parameters
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


// #include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase_vel_tan.h"
#include "tension.h"
#include "f_BC_vel_tan.h"

// Vectors, fields...
scalar omega[]; // vorticity field
vector h[]; // height functions vector
scalar curv_viz[]; // curvature 

scalar f0[]; // auxiliary volumic fraction field used for BC 
vector n_front[]; // normal vector of the interface
scalar alpha_front[]; // intercept linked to n_front[] during VoF-reconstruction 


/** 
## Boundary Conditions

At first, we modify the default BCs for the normals at the interface 
crossing borders by readily calculating the proper ones, since the 
interface is an oblique straight line:
*/

  /* Normal vectors of the interface and related intercepts BC */

// !!! Normal vectors are renormalized with the L1-norm !!!
// pay attention to the "foreach_dimension()" declaration for the normal components
n_front.y[bottom] = sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)));
n_front.x[bottom] = cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)));
alpha_front[bottom] = plane_alpha (f[ghost], (coord){
    sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
    cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
});


n_front.x[left] = sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)));
n_front.y[left] = cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)));
alpha_front[left] = plane_alpha (f[ghost], (coord){
    sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
    cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
});

/** 
These definitions are directly put in this user file for a better 
visibility, but are actually managed in 
[vof_vel_tan.h](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vof_vel_tan.h). 

Then, we apply the auxiliary volume fraction field 
to the advected one in ghost cells:

*/


  /* Volumic Fractions BC */
f[left] = f0[-1,0];
f[bottom] = f0[0,-1];

/** 
The rest of BCs remains the same as in the *naive* test case 
(see above figure): 
*/

  /* Contact Angles */
h.t[left] = contact_angle (THETA_0 - pi/2.) ; 
h.t[bottom] = contact_angle (BETA_0) ;


  /* Velocity BC */

u.n[left] = dirichlet ( - U0 * cos(BETA_0) ) ;
u.t[left] = dirichlet ( U0 * sin(BETA_0) ) ;
u.n[right] = dirichlet ( - U0 * cos(BETA_0) ) ;
u.t[right] = dirichlet ( U0 * sin(BETA_0) ) ;

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



/** 
## Generic Events
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
    fraction(f0, ((- (x - X_OFFSET) * tan(BETA_0) + Y_OFFSET) - y) ) ; 
    } while (adapt_wavelet ({f0}, (double[]){1e-3}, MAXLEVEL).nf); 

    // Initialization of the ghot cells'auxiliary field: 
    foreach_boundary(bottom)
      f0[0,-1] = f_BC_bottom (f0[], f0[1,0], f0[-1,0]);
    foreach_boundary(left)
      f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);

    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
  #else
    fraction(f0, ((- (x - X_OFFSET) * tan(BETA_0) + Y_OFFSET) - y) ) ;

    foreach_boundary(bottom)
      f0[0,-1] = f_BC_bottom (f0[], f0[1,0], f0[-1,0]);
    foreach_boundary(left)
      f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);
  #endif
  // ---------------------------

  foreach(){
    u.x[] = - U0 * cos(BETA_0) ; 
    u.y[] = U0 * sin(BETA_0);
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
    adapt_wavelet ({f,u}, (double[]){1e-3,1e-3,1e-3}, MAXLEVEL, LEVEL-1);
  #endif

/** 
At each step, we apply to ghost cells the same auxiliary volume fraction field 
as the initialized one, to keep the interface steady at borders:
*/

  foreach_boundary(bottom)
    f0[0,-1] = f_BC_bottom (f0[], f0[1,0], f0[-1,0]); 
  foreach_boundary(left)
    f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);

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

/** 
## Outputs
*/

event end (t = T_END){}


event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}


event profiles (t = {0, 50, 100}){                 
  char filename[200] ;
  sprintf(filename, "vel_tan_correc_shape_t%g.dat", t) ;
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
  sprintf (fileup, "vel_tan_correc_u_p_t%g.dat", t);
  sprintf (fileinterf, "vel_tan_correc_interf_t%g.dat", t);
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
Finally, we compute the curvature value of the interface 
at the bottom boundary:
*/


event bottom_curv (t += 1.){
  char filecurv[200] ;
  sprintf(filecurv, "vel_tan_correc_curv.dat") ;
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
  fflush(fcurv); // empty buffer
}

/**
~~~pythonplot Curvature value at the bottom axis
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=15) # size for the axes labels
label_size = 15
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 


t_c, x_c, curv_c = np.loadtxt("vel_tan_correc_curv.dat", unpack=True)
t_c_filt = t_c[np.abs(curv_c) > 0]
curv_c_filt = curv_c[np.abs(curv_c) > 0]

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

ax.plot(t_c_filt, np.abs(curv_c_filt), lw=1.2, label=r'With corrections [AMR]', 
        color=mapcolors[12])

ax.legend(frameon=False, loc=4)
plt.savefig('vel_tan_correc_curv.svg') 
~~~

As we can see, the curvature computation has been nicely improved, 
with a reduction by a factor of **25** compared to the 
[*naive* test case](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_no_correction.c#curvature-evolution).
That is to say a value of $2.10^{-2}$ instead of $5.10^{-1}$.
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
           max = 0.001, min = -0.001,
           map = blue_white_red);
  draw_string(legend, 1, size = 20., lw = 2.1);
  save ("movie_AMR_v2.mp4");
}