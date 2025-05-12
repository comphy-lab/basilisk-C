/** 
# Convergence test for inviscid collapsing cavities 

We want to verify that for inviscid collapsing cavities (see movie below), 
the following two properties are met:

  1. a capillary inertial scale invariant regime on approximately $2$ decades;
  2. the kinematic invariance of the inviscid flow. 

In order to do so, data are collected upon the axial velocity $\tilde{u}_{axis}$ 
at the apex of the cavity, maximum pressure $\tilde{p}_{max}$ in the domain, 
and axial curvature $\tilde{\kappa}$, and we plot them against the time to the 
singularity $\tilde{t}_0$. 

![Collapsing cavity | *top:* pressure | *bottom:* velocity field | $N_{max} = 8$](collapsing_data/invisc_data/invisc-conv_s10000_N8_m50_tinv1000.mp4)(width="40%")

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Important note**

All the notations employed, justification of the non-dimensionalization used and 
other general considerations are detailed in the [`README`](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README) 
file of this sandbox's directory.
</div>

## Numerical protocol 

![How to prepare a collapsing interface numerically](img_collapsing/collapsing_time_reversal_v2.png){width="30%"}

Our goal is simulating the convergence of a collapsing cavity towards a 
*singular state* corresponding to a perfect cone. To prepare a numerical 
configuration in which the interface is collapsing, we can exploit the conclusions 
of this [documentation's paragraph](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README#how-can-we-numerically-prepare-a-singular-state): 

  1. starting from an obtuse liquid cone of angle $\theta_0 > 90$° initially at 
  rest and placed in a [dipolar flow](http://basilisk.fr/sandbox/cailler/sierou_lister/README#dipolar-distribution-and-physical-meaning) 
  of strength $\widetilde{\mu}_0 > 0$ large enough to counteract the surface 
  tension driven recoil of the interface, the latter will *move backward* 
  compared to its initial state (*cf.* left side of the above figure);
  2. at a given time $\tilde{t} = \tilde{t}_{inv}$ is applied a 
  **kinematic reversal** $\mathbf{\tilde{u}} \to - \mathbf{\tilde{u}}$, 
  $\widetilde{\mu}_0 \to - \widetilde{\mu}_0$ (*right figure*), so that the 
  backward movement of the cone becomes a *finite-time singular collapse*, 
  with $\tilde{t}_0$ the time of curvature reversal (singularity), according to 
  [Sierou \& Lister, (2004)](#sierou2004);
  3. for an inviscid flow, we should verify the *kinematic invariance* of this 
  numerical experiment, *i.e.* at $\tilde{t} = 2 \, \tilde{t}_{inv} \approx \tilde{t}_0$ 
  there is the **cone reformation** (*see also the movie*).


## Code

### General Parameters

As demonstrated in the 
[`README`](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README#inviscid-non-dimensionalization) 
documentation, the inviscid problem doesn't have any characteristic length scale, 
so we can choose an arbitrary box size for our simulation. Let's take $L_0 = 10^4$:
*/

#define SIZE 1.e4

/** 
Preliminary studies have shown that a time reversal of $\tilde{t}_{inv} = 1000$ 
presents a good compromise between minimizing boundary effects and the number of 
decades revealing the capillary-inertial regime of the collapse:   
*/

#define T_INV 1.e3

/** 
We expect the complete reformation of the cone to be very close to twice the 
time reversal:
*/

#define T_END 2.2*T_INV

/** 
We conduct a spatial convergence study for maximum levels of the grid resolution 
belonging to the interval 
$N_{max} \in [7 \, ; 10] \Rightarrow \widetilde{\Delta}_{max} \in [9.8 \, ; 78.1]$. 
Due to the `Basilisk` server limitations for running codes in the *sandbox*, we 
limit the present code to the case where $N_{max} = 7$, and provide all the other 
needed data in the
[`invisc_data`](http://basilisk.fr/sandbox/cailler/collapsing_cavities/collapsing_data/invisc_data) 
directory. **With a lower minimum level, the simulation cannot converge.**
*/

#define LEVEL 7
#define MAXLEVEL 7 

/** 
An open angle of $\theta_0 = 120$°$> 90$° is chosen for the obtuse liquid cone, 
which coincides with the experimental value of [Zeff *et al.*, (2000)](#zeff2000). 
We also define a $\beta_0$ angle corresponding to the opposite angle of the apex. 
*/

#define THETA_0 120.*pi/180.
#define BETA_0 pi - THETA_0

/** 
Since we are also interested in the formation of liquid jets, and following the 
conclusions of the page
[`sierou_neg_dipolar_flow.c`](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_neg_dipolar_flow.c#discussion), 
strong values for the dipolar flow field are needed. We consequently choose 
$|\widetilde{\mu}_0| = 50$:
*/

#define MU_0 (t < T_INV ? 50.e0 : 50.e0) // dipole potential strength

/** 
An $x$--offset is necessary in order to avoid as much as possible boundary effects 
when the interface is recoiling/moving forward: 
*/

#define Y_OFFSET 0. 
#define X_OFFSET 0.4*SIZE 


/** 
At $\tilde{t}_{inv} = 1000$, the [mesh radius of the refined disk](#mesh-strategy) 
refining the growing cone apex should have approximately a ratio of $2/5$ with 
the size box, hence $\rho_\Delta = L_0 \times 2/5 = 4000$. 
In that case, the initial radius of the refined mesh can be inferred: 
$\rho_\Delta^0 = \rho_\Delta / {\tilde{t}_{inv}}^{2/3} = 40$, accordingly with 
the capillary-inertial dynamics.
Indeed, the disk radius evolves according to the capillary-inertial 
self-similar dynamics, growing before reaching the time of reversal, then 
decreasing when reforming the cone.
*/

#define R_0 40.
#define RADIUS_REFINE ( t < 1. ? R_0 : t < T_INV ? R_0*pow(t, 2./3.) : \
                        t < 2.*T_INV - 1. ? R_0*pow(2*T_INV - t, 2./3.) : \
                        t < 2.01*T_INV ? R_0 : R_0*pow(t - 2*T_INV, 2./3.) )


/**
If we want to make a movie, change the value below to $1$. 
To speed up the run on the `Basilisk` server, we will import a pre-registered 
file though. 
*/

#define MOVIE 0

/** 
We import all the needed libraries, especially those updating volume fractions 
in ghost cells to improve the treatment of moving contact angles for non-closed, 
moving interfaces, but also functions implementing the analytical expressions of 
the dipolar far-field flow:
*/                        

#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase_collapse.h" // change of direction for the interf. normal
#include "tension.h"
#include "f_BC_collapse.h" // updating vol. frac. ghost cells
#include "u_BC_collapse.h" // func. for dipolar flow
#include "tag.h" // to remove erroneous micro-bubbles

// Vectors, fields...
scalar f0[];
vector h[];
vector n_front[];
scalar alpha_front[];

scalar omega[];
scalar curv_viz[]; // to store curvature 
scalar u_rs[]; // L2-norm of velocity `u`


/** 
### Boundary Conditions 

**The geometric considerations, initial and boundary conditions are the same than 
those presented and *fully documented* in** 
[**`sierou_dipolar_flow.c`**](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_dipolar_flow.c#boundary-conditions), 
except for the box size obviously. Therefore, we instantiate the following domain: 

![Initial setup of the simulation](img_collapsing/collapsing_BC.png){width="30%"}

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

// Invariance of time in the equation:
event reversal_BC_velocity (i++) {
  if (t < T_INV){ // recoiling phase
    // On left boundary, symmetry conditions are imposed for the external fluid.
    // WARNING! The interface can jump from top to left border!
    u.n[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uxA(x,y) );
    u.t[left] = (1. - f[])*neumann(0.) + f[]*dirichlet ( uyA(x,y) );

    // On top boundary, ext. & int. fluids coexist for most of the simulation:
    u.n[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uyA(x,y) );
    u.t[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet ( uxA(x,y) );

    // On right boundary (will always be liquid phase):
    u.n[right] = dirichlet ( uxA(x,y) );
    u.t[right] = dirichlet ( uyA(x,y) );
  }
  else { // collapsing phase
    // On left boundary:
    u.n[left] = (1. - f[])*neumann(0.) + f[]*dirichlet (- uxA(x,y) );
    u.t[left] = (1. - f[])*neumann(0.) + f[]*dirichlet (- uyA(x,y) );

    // On top boundary:
    u.n[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet (- uyA(x,y) );
    u.t[top] = (1. - f[ghost])*neumann(0.) + f[ghost]*dirichlet (- uxA(x,y) );

    // On right boundary:
    u.n[right] = dirichlet (- uxA(x,y) );
    u.t[right] = dirichlet (- uyA(x,y) );
  }
}

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
  mu1 = 0, mu2 = 0;

  run() ;
}


event init (t = 0){
  fraction(f0, (y - (- (x) * tan(BETA_0) + Y_OFFSET) ) ) ; 

  foreach_boundary(left)
    f0[-1,0] = f_BC_left (f0[], f0[0,1], f0[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (f0[], f0[1,0], f0[-1,0]);

  foreach(){
    f[] = f0[];
    u.x[] =  f[]*uxA(x,y) + (1. - f[])*uxO(x,y); 
    u.y[] =  f[]*uyA(x,y) + (1. - f[])*uyO(x,y);
  }

  // computation of curvature for visualization:
  boundary({f});

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});

  curvature (f, curv_viz, 1, add = false);
}

event f_BC_update (i++) {
  foreach_boundary(left)
    f0[-1,0] = f_BC_left (f[], f[0,1], f[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (f[], f[1,0], f[-1,0]);

  boundary({f});

  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
}

event vorti (i++){
  // cf. bug report: https://groups.google.com/g/basilisk-fr/c/ok-OhtzO1Pk
  // vorticity:
  foreach()
    omega[] = (u.y[1,0]-u.y[-1,0]+u.x[0,-1]-u.x[0,1]) / (2.*Delta + SEPS);

  // L2-norm of `u`:
  foreach()
    u_rs[] = sqrt( sq(u.x[]) + sq(u.y[]) );

  // curvature:
  curvature (f, curv_viz, 1, add=false);
}



/** 
### Mesh Strategy

The most difficult part of this simulation is to collect clean data for fields 
quite unstable like the pressure or the curvature. Therefore, a refinement on the 
apex of the liquid phase is needed, because it is where all the physics is 
concentrated. Unfortunately, the preparation of a collapsing situation is 
directly linked with the curvature radius of the apex, which is maximum at 
$\tilde{t}_{inv}$. 
Furthermore, the use of a classical adaptive mesh refinement algorithm 
doesn't work well to capture the fields of interest, and may create "dendrites" 
of refinement in regions of poor physical interest, slowing down massively the 
run.

Hence, to avoid meshing too much regions of the whole simulation box 
that would slow down the computations, we need to find a way to balance grid 
resolution at the region of interest and computation speed. The strategy adopted 
is to develop a *manual mesh* that:

  1. locates the axial position of the interface (*apex*);
  2. refine the cone apex in the liquid phase within a certain range.

Additionally, for inviscid simulations, we need to **unrefine the gas phase** 
because:

  * we need to ensure the *kinematic invariance* by removing as much as 
  possible diffusion/viscosity effects;
  * it allows to match a "one-fluid" model (void/water).

Finally, the intersection point of the interface crossing the border opposite 
to the cone apex needs special treatment too: without a sufficient level of 
refinement in the vicinity of this intersection point, the *Poisson* solver diverges, 
so the location of this moving point is also needed to refine its region within a 
given radius.  
*/

event adapt (i++) {
  scalar pos[];
  scalar pos_top[];
  scalar pos_left[];
  coord D = {1., 0.};
  coord E = {0., 1.};
  double xorigin = -1e100;
  double xorigin_top = -1e100;
  double yorigin_top = 0.;
  double xtemp = 0.;
  double xtemp_top = 0.;
  position (f, pos, D, add=false);
  position (f, pos_top, D, add=false);
  position (f, pos_left, E,  add=false);
  foreach() {
    if ((y <= Delta) && interfacial(point, f)){
      xtemp = pos[];
      if (xtemp > xorigin){
        xorigin = xtemp;
      }
    }
    if ((y >= SIZE - Delta) && interfacial(point, f)){
      xtemp_top = pos_top[];
      yorigin_top = y;
      if (xtemp_top > xorigin_top){
        xorigin_top = xtemp_top;
      }
    }
    if ((x <= - X_OFFSET + Delta) && interfacial(point, f)){
      xorigin_top = x;
      yorigin_top = pos_left[];   
    }
  }

  /** For having the disk refined mesh ONLY in the fluid at the cone apex 
  and the opposite point of intersection with the simulation box: */
  refine ( (sqrt(sq (x - xorigin) + sq (y)) < RADIUS_REFINE) && level < MAXLEVEL);
  refine ( (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) < 0.5*RADIUS_REFINE) 
    && level < 7);
  
  /** We unrefine the gas phase, with the exception of the intersection point:*/
  unrefine ((f[] < 0.1) && (level > LEVEL - 1 ) 
    && (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) > 0.5*RADIUS_REFINE) );

  /** 
  We unrefine in the liquid phase all the regions not within the apex and 
  intersection's ranges:
  */ 
  unrefine ( (f[] > 0.1) && (sqrt(sq (x - xorigin) + sq (y)) > RADIUS_REFINE) 
    && (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) > 0.5*RADIUS_REFINE) 
    && (level > LEVEL - 1 ));



/** 
In the far-field, there exists a movement of the interface resulting from 
the balance between the *capillary flow* coming from the *Laplace* pressure 
gradient and the *dipolar flow*. In that case, it is mandatory 
to update at each timestep the volume fractions in ghost cells belonging to 
the borders where the interface is moving/is susceptible to move. 

The update shall be *dynamical* to allow the translation of the interface 
along the borders, that is to say, feeding the auxiliary volume fraction field 
`f0[]` by applying the functions `f_BC_myborder()` to `f[]` 
(not to `f0[]` which would freeze any free displacement at contact angle 
$\beta_0$; we do not want a pinned point). 
*/ 

  foreach_boundary(left)
    f0[-1,0] = f_BC_left (f[], f[0,1], f[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (f[], f[1,0], f[-1,0]);

  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0}); // for boundary conditions on levels

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

  foreach_boundary (bottom) {
    if (f[]<1 && f[1]>=1 && f[-1]>= 1) f[]=1;
    if (f[]>0 && f[1]<=0 && f[-1]<= 0) f[]=0;
  }
  foreach_boundary (top) {
    if (f[]<1 && f[1]>=1 && f[-1]>= 1) f[]=1;
    if (f[]>0 && f[1]<=0 && f[-1]<= 0) f[]=0;
  }
  foreach_boundary (left) {
    if (f[]<1 && f[0,1]>=1 && f[0,-1]>= 1) f[]=1;
    if (f[]>0 && f[0,1]<=0 && f[0,-1]<= 0) f[]=0;
  }
}

/** 
### Time Reversal
*/

event reversal (t=T_INV) {
  foreach_face(){
    uf.x[] = 0.;
  }
  foreach(){
    u.x[] *= -1.;
    u.y[] *= -1.;
  }
} 


/** 
### Outputs
*/

event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}

event end (t = T_END){}

/** 
We need some useful functions to store our data: 
*/

// Functions adapted for parallel computation:
double maximum_func (scalar a) {
  double maxi = - 1e100;
  foreach (reduction(max:maxi))
    if ((x < 0.38*SIZE) && (y < 0.4*SIZE)) // to avoid spurious effects on the top contact angle
      if (fabs(a[]) > maxi){
        maxi = fabs(a[]);
      }
  return maxi;
}

double minimum (scalar a) {
  double mini = 1e100;
  foreach (reduction(min:mini))
    if (y < 0.4*SIZE) // to avoid spurious effects on the top contact angle
      if (fabs(a[]) < mini)
        mini = fabs(a[]);
  return mini;
}

// cf. http://basilisk.fr/sandbox/ecipriano/src/common-evaporation.h
double avg_neighbor (Point point, scalar Y, scalar f) {
  /* Compute the average value of a scalar field Y 
  in a 3x3 stencil around the current cell */
  double fYnei = 0., fnei = 0.;
  foreach_neighbor(1) { // the `1` option is for enabling the 3x3 stencil
    if (Y[] < 1.e3){ // to avoid NaN values
      fYnei += f[]*Y[];
      fnei  += f[];
    }
  }
  return fYnei / fnei;
}

/** 
We store the maximum pressure at each timestep, along with the average curvature 
value on the axis, and the axial velocity: 
*/

event datafiles (i++){
  char filename[200] ;
  sprintf(filename, "invisc-conv_lvl%d_m%g_pressure_s%g_tinv%g.dat", 
    MAXLEVEL, MU_0, SIZE, T_INV) ;
  static FILE * fp = fopen(filename, "a") ;
  double pmax = maximum_func (p);
  double pmin = minimum (p);
  fprintf(fp,"%g %.6lf %.6lf\n", t, fabs(pmax - pmin), pmax);

  char filename2[200] ;
  sprintf(filename2, "invisc-conv_lvl%d_m%g_x-ux-curv_s%g_tinv%g.dat", 
    MAXLEVEL, MU_0, SIZE, T_INV) ;
  static FILE * fp2 = fopen(filename2, "a") ;
  double xpos = 0., curv = 0., u_axis = 0.;
  foreach_boundary (bottom, reduction(+:xpos) reduction(+:curv) reduction(+:u_axis)){
    if (interfacial (point,f)){ // if we are at the cone tip
      if ( (f[] > 0.1) && (f[] < 0.9) ){
        xpos += x;
        curv += avg_neighbor (point, curv_viz, f);
        u_axis += u.x[];
      }
    }
  }
  fprintf(fp2, "%g %g %g %g\n", t, xpos, u_axis, curv);

  fflush(fp); // empty buffer
  fflush(fp2); // empty buffer
}


/**
~~~pythonplot Spatial resolution
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Use LaTeX fonts in the plot:
plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=20) # size for the title
plt.rc('axes', labelsize=14) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 


lvl_lst = [7, 8, 9, 10]
resol = ['78.13', '39.06', '19.53', '9.77']
t_inv = 1000
markers = ["o", "v", "s", "d"]

t_0_lst = [2010, 1955, 1960, 1970]
#t_0_lst = [2040, 2010, 1985, 1987]
Ncolors = 5
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# U_AXIS
# ------
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e1, 1e3)
ax.set_ylim(1e-1, 1e1)
ax.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax.set_ylabel(r'$|\tilde{u}_{axis}|$')


for k, lvl in enumerate(lvl_lst):
  mu_0 = 50
  L_0 = 10000
  tinv = 1000
  if lvl == 7:
    path = 'invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)
  else:
    path = '../collapsing_data/invisc_data/invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)
  t, xpos, u_axis, curv = np.loadtxt(path + 'tinv{}.dat'.format(t_inv), unpack=True)
  # Let's filter the data where the curvature is not zero:
  t_filt = t[curv != 0]
  t_filt2 = t_filt[t_filt > t_inv]
  u_axis_filt = u_axis[curv != 0]
  u_axis_filt2 = u_axis_filt[t_filt > t_inv]
  ax.scatter(t_0_lst[k] - t_filt2, np.abs(u_axis_filt2), 
             marker = markers[k],
             facecolors='none', 
             edgecolors=mapcolors[k],
             label=r'$N_{max} = \,\,$' + r'${}, \,\,$'.format(lvl) + r'$\widetilde{\Delta}_{max} = \,\,$' + r'${}$'.format(resol[k])
             )

t_array = np.linspace(10, 5000, 200)
u_axis_scaling = 5.98762*np.power(t_array, -1./3.)
ax.plot(t_array, u_axis_scaling, '-', lw=0.8, color='darkred',
        label=r'$\tilde{t}^{-1/3}$' )

ax.legend(frameon=False, loc=3)


# PRESSURE
# --------
ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_xlim(1e1, 1e3)
ax2.set_ylim(1e-2, 1e0)
ax2.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax2.set_ylabel(r'$\tilde{p}_{max}$')

t_0_lst = [2040, 2010, 1985, 1988]


for k, lvl in enumerate(lvl_lst):
  mu_0 = 50
  L_0 = 10000
  if lvl == 7:
    path = 'invisc-conv_lvl{}_m{}_pressure_s{}_'.format(lvl, mu_0, L_0)
  else:
    path = '../collapsing_data/invisc_data/invisc-conv_lvl{}_m{}_pressure_s{}_'.format(lvl, mu_0, L_0)
  t, pdiff, pmax = np.loadtxt(path + 'tinv{}.dat'.format(t_inv), unpack=True)
  # Let's filter the data where the curvature is not zero:
  t_filt = t[t > t_inv+60]
  pmax_filt = pmax[t > t_inv+60]

  ax2.scatter(t_0_lst[k] - t_filt, pmax_filt, 
             marker = markers[k],
             facecolors='none', 
             edgecolors=mapcolors[k],
             label=r'$N_{max} = \,\,$' + r'${}, \,\,$'.format(lvl) + r'$\widetilde{\Delta}_{max} = \,\,$' + r'${}$'.format(resol[k])
             )


t_array = np.linspace(10, 5000, 200)
pmax_scaling = 5.2*np.power(t_array, -2./3.)
ax2.plot(t_array, pmax_scaling, '-', lw=0.8, color='darkred',
        label=r'$\tilde{t}^{-2/3}$' )

ax2.legend(frameon=False, loc=3)
plt.savefig('prelim_study_uaxis+pressure_Nconv.svg')
~~~
*/

/**
~~~pythonplot Curvature convergence
fig, ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e1, 1e3)
ax.set_ylim(1e-4, 1e-1)
ax.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax.set_ylabel(r'$|\tilde{\kappa}|$')

#t_0_lst = [2010, 1955, 1960, 1970]
t_0_lst = [1950, 1950, 1950, 1952]

for k, lvl in enumerate(lvl_lst):
  mu_0 = 50
  L_0 = 10000
  if lvl == 7:
    path = 'invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)
  else:
    path = '../collapsing_data/invisc_data/invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)
  t, xpos, u_axis, curv = np.loadtxt(path + 'tinv{}.dat'.format(t_inv), unpack=True)
  # Let's filter the data where the curvature is not zero:
  t_filt = t[curv != 0]
  t_filt2 = t_filt[t_filt > t_inv]
  curv_filt = curv[curv != 0]
  curv_filt2 = curv_filt[t_filt > t_inv]
  ax.scatter(t_0_lst[k] - t_filt2, np.abs(curv_filt2), 
             marker = markers[k],
             facecolors='none', 
             edgecolors=mapcolors[k],
             label=r'$N_{max} = \,\,$' + r'${}, \,\,$'.format(lvl) + r'$\widetilde{\Delta}_{max} = \,\,$' + r'${}$'.format(resol[k])
             )

t_array = np.linspace(10, 5000, 200)
curv_scaling = 0.229801 *np.power(t_array, -2./3.)
ax.plot(t_array, curv_scaling, '-', lw=0.8, color='darkred',
        label=r'$\tilde{t}^{-2/3}$' )


# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.55, box.height])

# Put a legend to the right of the current axis
ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig('prelim_study_curv_Nconv.svg')


~~~

We can notice that for almost $2$ temporal decades, the axial velocity and the 
maximum pressure in the domain follow the theorized capillary-inertial power laws, 
but only starting from $N_{max} = 10$. Indeed, the curvature criterion for 
numerical convergence:

$$
|\widetilde{\kappa}| < \dfrac{1}{2 \, \widetilde{\Delta}_{max}} 
$$

is not verified for levels of resolution lower than $N_{max} = 9$ and times 
close to the singularity. *Numerical viscosity* occurs at these times and the 
simulation is not anymore valid.  

The curvature is subject to strong oscillations and fluctuations due to the 
poor number of measured points on the interface (3 points) to average the axial 
curvature (*is there a better `Basilisk`--way of coding to compute this quantity 
rather than using the `avg_neighbor` function?*). Also, the need for unrefining 
the gas phase in the inviscid simulation creates a sharp transition in the 
resolution at the apex, and could explain the multiple curvature sign changes 
observed during the run.
*/


/**
We also check the **kinematic invariance** for the finer grid resolution, 
by comparing the preliminary phase of recoiling $(\tilde{t} \leqslant \tilde{t}_{inv})$ 
and the collasping phase $(\tilde{t}_{inv} < \tilde{t} \leqslant \tilde{t}_0 \approx 2 \, \tilde{t}_{inv})$. 
The following plots are thus tracking the axial velocity and maximum pressure 
for the recoiling and collapsing phases:

~~~pythonplot Kinematic Invariance [N~max~ = 10]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6))

lvl = 10
mu_0 = 50
L_0 = 10000

# U_AXIS
# ------

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e1, 1e3)
ax.set_ylim(1e-1, 1e1)
ax.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax.set_ylabel(r'$|\tilde{u}_{axis}|$')

if lvl == 7:
  path = 'invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)
else:
  path = '../collapsing_data/invisc_data/invisc-conv_lvl{}_m{}_x-ux-curv_s{}_'.format(lvl, mu_0, L_0)

t_inv_lst = [1000, 1000]
t_0_lst = [1970, 1970]
colors = ['darkred', 'darkblue']
markers = ['s', 'o']

for k, t_inv in enumerate(t_inv_lst):
  t, xpos, u_axis, curv = np.loadtxt(path + 'tinv{}.dat'.format(t_inv), unpack=True)
  # Let's filter the data where the curvature is not zero:
  if k == 0:
    t_filt = t[curv != 0]
    t_filt2 = t_filt[t_filt < t_inv_lst[k]]
    u_axis_filt = u_axis[curv != 0]
    u_axis_filt2 = u_axis_filt[t_filt < t_inv_lst[k]]
    ax.scatter(t_filt2, np.abs(u_axis_filt2), 
              marker = markers[k],
              facecolors='none', 
              edgecolors=colors[k],
              label=r'$\tilde{t} < \tilde{t}_{inv} = \,\,$' + r'${}$'.format(t_inv)
              )
  else:
    t_filt = t[curv != 0]
    t_filt2 = t_filt[t_filt > t_inv_lst[k]]
    u_axis_filt = u_axis[curv != 0]
    u_axis_filt2 = u_axis_filt[t_filt > t_inv_lst[k]]
    ax.scatter(t_0_lst[k] - t_filt2, np.abs(u_axis_filt2), 
              marker = markers[k],
              facecolors='none', 
              edgecolors=colors[k], 
              label=r'$\tilde{t} > \tilde{t}_{inv} = \,\,$' + r'${}$'.format(t_inv)
              )
    
t_array = np.linspace(10, 1000, 200)
u_axis_scaling = 6.2*np.power(t_array, -1./3.)
ax.plot(t_array, u_axis_scaling, '-', lw=0.8, color='black',
        label=r'$\tilde{t}^{-1/3}$' )

ax.legend(frameon=False)


# PRESSURE
# --------
ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_xlim(1e1, 1e3)
ax2.set_ylim(1e-2, 1e0)
ax2.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax2.set_ylabel(r'$\tilde{p}_{max}$')

if lvl == 7:
  path = 'invisc-conv_lvl{}_m{}_pressure_s{}_'.format(lvl, mu_0, L_0)
else:
  path = '../collapsing_data/invisc_data/invisc-conv_lvl{}_m{}_pressure_s{}_'.format(lvl, mu_0, L_0)


t_inv_lst = [1000, 1000]
t_0_lst = [1990, 1990]
colors = ['darkred', 'darkblue']
markers = ['s', 'o']

for k, t_inv in enumerate(t_inv_lst):
  t, pdiff, pmax = np.loadtxt(path + 'tinv{}.dat'.format(t_inv), unpack=True)
  # Let's filter the data where the curvature is not zero:
  if k == 0:
    t_filt = t[t < t_inv_lst[k]]
    pmax_filt = pmax[t < t_inv_lst[k]]
    ax2.scatter(t_filt, pmax_filt, 
              marker = markers[k],
              facecolors='none', 
              edgecolors=colors[k],
              label=r'$\tilde{t} < \tilde{t}_{inv} = \,\,$' + r'${}$'.format(t_inv)
              )
  else:
    t_filt = t[t > t_inv_lst[k]]
    pmax_filt = pmax[t > t_inv_lst[k]]
    ax2.scatter(t_0_lst[k] - t_filt, pmax_filt, 
              marker = markers[k],
              facecolors='none', 
              edgecolors=colors[k], 
              label=r'$\tilde{t} > \tilde{t}_{inv} = \,\,$' + r'${}$'.format(t_inv)
              )
    
t_array = np.linspace(10, 5000, 200)
pmax_scaling = 5.4*np.power(t_array, -2./3.)
ax2.plot(t_array, pmax_scaling, '-', lw=0.8, color='black',
        label=r'$\tilde{t}^{-2/3}$' )

ax2.legend(frameon=False)
plt.savefig('prelim_study_kin_rev_tinv1000_uaxis+pmax.svg')
~~~

For most part of the simulation, both quantities have a recoiling and collapsing 
phases that are overlapped, except during the temporal range where numerical 
viscosity occurs (at the very beginning and late stages of the simulation), 
that is to say, when the simulation is under resolved. It explains why the reformed 
cone is not *exactly* the same as the perfect initial one (the aggregated errors 
from numerical viscosity generate small deformations of the interface, that 
prevent it to revert to its initial solution). 

However, when the spatial resolution is ensured for the numerical convergence, 
the *kinematic invariance* is verified.
*/










/** 
### Movie
*/

#if MOVIE
#include "view.h"
#include "cmaps.h"

event mov_p_u (t+=10.) {
  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  view (tx = -0.05, ty = -0.01, width = 1000., height = 600., fov = 8.);
  box();
  squares ("p", spread=-1, linear=true, map=cool_warm); // map=cool_warm
  draw_vof ("f", lw = 1.5);
  mirror ({0.,1}) {
    draw_vof ("f", lw = 1.5);
    squares("u_rs", spread=-1, linear=true, map=viridis);
    cells ();
  }
  draw_string(legend, 1, size = 20., lw = 2.1);
  save ("invisc-conv_s10000_N8_m50_tinv1000.mp4");
}
#endif


/** 
## References 

~~~bib
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
