/** 
# Crossing a finite-time singularity through viscous collapse

This file is an extension of the 
[inviscid collapse of a cavity](http://basilisk.fr/sandbox/cailler/collapsing_cavities/invisc_collapse_conv.c) 
when taking into account **viscosity**. We hence want to simulate the collapse 
of a conical cavity (such as a bursted bubble) to investigate the impact of 
viscous effects close to the finite-time singularity. Our goal is to 
understand and show if viscosity is a regularizing mechanism able to explain 
the "*ultraviolet cutoff*" mentioned in [Zeff *et al.*, (2000)](#zeff2000). 

![**Viscous** collapse (go to 1') | *top:* pressure | *bottom:* velocity field | $N_{max} = 9$](collapsing_data/visc_data/visc-collapse_s230_N7-9_tinv15.mp4)(width="40%")

![**Inviscid** collapse | *top:* pressure | *bottom:* velocity field | $N_{max} = 9$](collapsing_data/visc_data/invisc-collapse_s230_N8-9_tinv15.mp4)(width="40%")

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Nota Bene**

To understand most of the code below, it is recommended to have read the following 
*sandbox*'s pages, as a lot of notations employed, physical explanations and 
numerical strategies for setting parameters are already defined there:

  1. First, take a look to the [`README`](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README) 
  file of this directory, as it details most of the notations used and the main 
  logic behind this code;
  2. **Boundary conditions** and **initialization** of this file are well documented 
  in a far [simpler numerical experiment](http://basilisk.fr/sandbox/cailler/sierou_lister/sierou_dipolar_flow.c#boundary-conditions) 
  concerning the impact of [*dipolar flows*](http://basilisk.fr/sandbox/cailler/sierou_lister/README#dipolar-distribution-and-physical-meaning);
  3. As said earlier, the current code is very similar to the 
  [*inviscid* case](http://basilisk.fr/sandbox/cailler/collapsing_cavities/invisc_collapse_conv.c) 
  where, for example, the [*mesh strategy*](http://basilisk.fr/sandbox/cailler/collapsing_cavities/invisc_collapse_conv.c#mesh-strategy) 
  is explained.
</div>

## Code 

### General Parameters

We want to describe the collapsing phase of a conical cavity for times 
very close to the finite-time singularity, *i.e.*, when $t - t_0 \sim t_\mu$. 
We use the [capillary-inertial non-dimensionalization](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README#inviscid-non-dimensionalization) 
by introducing the viscous length scale $L \equiv \ell_\mu$, along with the 
following numerical parameters for air/water:

$$
\rho_l^{num} = 1 
\,\, ; \,\,
\rho_g^{num} = 10^{-3}
\,\, ; \,\,
\mu_l^{num} = 1 
\,\, ; \,\,
\mu_g^{num} = 10^{-2}
\,\, ; \,\,
\sigma^{num} = 1
$$

so that:

$$
\dfrac{\rho_g}{\rho_l} = 10^{-3} \; ; \; \dfrac{\mu_g}{\mu_l} = 10^{-2}
$$

and:

$$
\ell_\mu = 1 \quad ; \quad t_\mu = 1
$$

Thus, in order to capture a transient between the capillary-inertial collapse 
and the viscous regime (if it exists), the simulation domain must have a size 
$L_0 \gg \ell_\mu$, *i.e.* $L_0 \sim 10^2 \ell_\mu  \sim 10^2$, so that the 
curvature radius of the interface at the kinematic reversal (when $t \equiv t_{inv}$) 
is far more larger than the viscous length. 
In the same time, we want to capture viscous effects, so we need a grid resolution 
at least ten times smaller than the viscous scale: $\Delta < 10 \, \ell_\mu$.

In practice, we are using the following parameters:

  1. for $\boxed{t_{inv} = 15 \, t_\mu}$ and $|\widetilde{\mu_0}| = 50$, the 
  curvature radius is roughly $\rho_\kappa \approx 50 \, \ell_\mu \gg \ell_\mu = 1$; 
  2. hence, a box size $\boxed{L_0 = 230 \, \ell_\mu}$ (corresponding to five times 
  the maximum curvature radius of the interface during the simulation) is chosen 
  to offer a good compromise between the spatial resolution needs 
  (therefore the CPU running time) and the boundary effects^[N.B.: these values 
  can be easily assessed by running the code at very low resolution and seeing 
  with *BView* where the interface stops on the axis of symmetry at $t \equiv t_{inv}$.];
  3. to ensure capturing all potential viscous effects, we decide to have around 
  twenty points below the viscous length scale. With a maximum level of refinement 
  $\boxed{N_{max} = 12}$, the maximum resolution is 
  $\Delta_{max} = L_0/2^{N_{max}} = 0,056 \, \ell_\mu$, so 
  $\lfloor \ell_\mu / \Delta_{max} \rfloor = 18$ which will be enough. 

As set out in the [`README`](http://basilisk.fr/sandbox/cailler/collapsing_cavities/README) 
file, the same parameters are used for inviscid simulations to compare the 
results (with the obvious exception of dynamic viscosities set to zero then).
*/


#define SIZE 2.3e2 
#define T_INV 1.5e1 
#define T_END 2.*T_INV+5.

#define THETA_0 120.*pi/180.
#define BETA_0 pi - THETA_0
#define MU_0 (t < T_INV ? 50.e0 : 50.e0) // dipole potential strength
#define MU_LIQ 1. // liquid dynamic viscosity
#define MU_GAS 1.e-2 // gas dynamic viscosity
#define Y_OFFSET 0. 
#define X_OFFSET 0.4*SIZE 

/** The minimum resolution level is set to $7$, otherwise lower values would 
result in a divergence of the *Poisson* solver (we need enough resolution on 
boundaries). */
 
#define LEVEL 7
 
/** As for the maximum resolution, the recoiling phase starts with a medium 
resolution value $(N_{max} = 9)$ to speed up this necessary but quite annoying 
preparatory phase for initiating the collapse, without having an impact on the 
whole simulation. Then, up to a certain time limit really close to the expected 
time of singularity, we increase by one level the grid resolution $(N_{max} = 10)$, 
in order to have a smooth connection with the last stage of the collapse, where 
viscous effects are expected to occur and the level of refinement switches to 
$N_{max} = 12$ .

~~~C
#define MAXLEVEL (t < T_INV ? 9 : t < 2*T_INV-5. ? 10 : 12)
~~~

However, to run the code on the `Basilisk` server, a very low-resolution 
configuration is used:
*/

#define MAXLEVEL 7
 
/** 
The disk radius evolves according to the capillary-inertial 
self-similar dynamics, growing before reaching the time of reversal, then 
decreasing when reforming the cone. The values set are optimized with respect 
to the time of singularity and reduction of computation time.
*/
#define R_0 15.
#define RADIUS_REFINE ( t < 1 ? R_0 : t < T_INV ? R_0*pow(t, 2./3.) : \
                        t < 2.*T_INV - 1. ? R_0*pow(2*T_INV - t, 2./3.) : \
                        t < 2.*T_INV + 2. ? R_0 : R_0*pow(t - 2*T_INV, 2./3.) )

/** 
Then we choose or not to store all the physical quantities of interest, make a 
movie, or if we need to refine (viscous simulation) or not (invisicd simulaiton) 
the gas phase:
*/
#define STOCK_DATA 1
#define MOVIE 0
#define GAS_REFINEMENT 1 // 1/ for a viscous simulation | 0/ for an inviscid one

/** 
We import all the needed libraries, especially those updating volume fractions 
in ghost cells to improve the treatment of moving contact angles for non-closed, 
moving interfaces, but also functions implementing the analytical expressions of 
the dipolar far-field flow:
*/    
#include "axi.h"
#include "navier-stokes/centered.h"

#if GAS_REFINEMENT
  #include "navier-stokes/double-projection.h" // improve pressure field stability
#endif

#include "contact.h"
#include "two-phase_collapse.h" // change of direction for the interf. normal
#include "tension.h"
#include "f_BC_collapse.h"
#include "u_BC_collapse.h"
#include "tag.h" // to remove erroneous micro-bubbles
#include "my-curvature_2D.h" // computing 2D curvature (to compare it with 3D-AXI)
#include "okubo-weiss.h" 

/** 
**A huge difference with the inviscid case is to rely on the 
[*double approximate projection method*](http://basilisk.fr/src/navier-stokes/double-projection.h).** 
Indeed, previous tests have shown that without a better estimate of the pressure 
field, *cusps* happened when approaching the conical singularity, slowing down 
massively the run while it was safe to say it was not a physical phenomenon but 
a purely numerical one.
*/ 


// Vectors, fields...
scalar f0[];
vector h[];
vector n_front[];
scalar alpha_front[];

scalar omega[];
scalar curv_viz[]; // to store curvature 
scalar curv_2D[];
scalar visc_dissip[];
scalar lambda2[]; // Okubo-Weiss criterion
scalar lambda2_rs[]; // rescaling over the range [-1;1]
scalar x_interf[];
scalar y_interf[];
scalar u_rs[]; // L2-norm of velocity `u`
scalar duz_dz[];
scalar p_interf[];
scalar x_p_interf[];
scalar y_p_interf[];


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

  #if GAS_REFINEMENT
    mu1 = MU_LIQ, mu2 = MU_GAS ; // enable viscosity
  #else
    mu1 = 0, mu2 = 0;
  #endif

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
  curvature_2D (f, curv_2D, 1, add=false);


  #if GAS_REFINEMENT
    // Okubo-Weiss criterion:
    compute_lambda2_v2 (lambda2, u);
    rescaled_lambda2 (lambda2_rs, lambda2);

    // viscous dissipation in cylindrical coordinates:
    foreach(){
      visc_dissip[] = (
        2.*(
          sq(u.x[1] - u.x[-1]) + sq(u.x[0,1] - u.x[0,-1])
          + sq(u.y[1] - u.y[-1]) + sq(u.y[0,1] - u.y[0,-1])
          + 2.*(u.y[1] - u.y[-1])*(u.x[0,1] - u.x[0,-1])
          )/sq(2.*Delta)
        + 2.*sq(u.y[]/y)
        )*(f[]*MU_LIQ + (1. - f[])*MU_GAS);
    }
  #endif
}


event get_p_interf(i++){
  coord D = {1., 0.};
  coord E = {0., 1.};
  position (f, x_p_interf, D, add=false);
  position (f, y_p_interf, E, add=false);
  foreach() {
    double xp = x_p_interf[] - 5.*Delta*n_front.x[];
    double yp = y_p_interf[] - 5.*Delta*n_front.y[]; 
    p_interf[] = interpolate (p, xp, yp);
  }
}


/** 
### Mesh Strategy

The strategy used for meshing the domain through time is explained 
[*here*](http://basilisk.fr/sandbox/cailler/collapsing_cavities/invisc_collapse_conv.c#mesh-strategy). 
See [*this movie*](http://basilisk.fr/sandbox/cailler/collapsing_cavities/invisc_collapse_conv.c#convergence-test-for-inviscid-collapsing-cavities) 
for having a visual idea of this strategy.
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


  /** For having the disk refined mesh in both phases (liq. + gas) at the cone apex 
    and the opposite point of intersection with the simulation box: */
  refine ( (sqrt(sq (x - xorigin) + sq (y)) < RADIUS_REFINE) && level < MAXLEVEL);
  refine ( (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) < 0.5*RADIUS_REFINE) 
    && level < 7);

  /** For the viscous simulation, the gas phase has to be refined, because the 
  generated vorticity behind the interface in the gas phase is physical and 
  therefore can have a role to play in the last stages of the collapse:
  */
  #if GAS_REFINEMENT
    /** We unrefine in both phases all the regions not within the apex and 
    intersection's ranges: */ 
    unrefine ( (sqrt(sq (x - xorigin) + sq (y)) > RADIUS_REFINE) 
      && (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) > 0.5*RADIUS_REFINE) 
      && (level > LEVEL - 1 ));

  /** For inviscid simulations, on the contrary, we need to **unrefine the gas phase** 
because:

  * we need to ensure the *kinematic invariance* by removing as much as 
  possible diffusion / numerical viscosity effects that are responsible for 
  significant deformations of the interface close to the singularity;
  * it allows to match a "one-fluid" model (void/water). */
  #else  
    /** We unrefine the whole gas phase, with the exception of the intersection point:*/
    unrefine ((f[] < 0.1) && (level > LEVEL - 1 ) 
      && (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) > 0.5*RADIUS_REFINE) );

    /** We unrefine in the liquid phase all the regions not within the apex and 
    intersection's ranges: */ 
    unrefine ( (f[] > 0.1) && (sqrt(sq (x - xorigin) + sq (y)) > RADIUS_REFINE) 
      && (sqrt(sq (x - xorigin_top) + sq (y - yorigin_top)) > 0.5*RADIUS_REFINE) 
      && (level > LEVEL - 1 ));
  #endif


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
<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**Nota Bene**

The following plots were all generated on a personal laptop for:


  #define MAXLEVEL (t < T_INV ? 9 : t < 2*T_INV-5. ? 10 : 12)


All the data used can be found in 
[*this directory*](http://basilisk.fr/sandbox/cailler/collapsing_cavities/collapsing_data/visc_data/).

</div>

*/


#if STOCK_DATA
  /** 
  We compare the evolutions for the curvature, the axial velocity and the maximum 
  pressure for two kinds of simulations: the viscous one and the inviscid one. 
  */

  event datafiles (i++){
    char filename[200] ;
    sprintf(filename, "visc-collapse_lvl%d_m%g_pressure_s%g_tinv%g.dat", 
      LEVEL, MU_0, SIZE, T_INV) ;
    static FILE * fp = fopen(filename, "a") ;
    double pmax = maximum_func (p);
    double pmin = minimum (p);
    fprintf(fp,"%g %.6lf %.6lf\n", t, fabs(pmax - pmin), pmax);


    char filename2[200] ;
    sprintf(filename2, "visc-collapse_lvl%d_m%g_x-ux-curv_s%g_tinv%g.dat", 
      LEVEL, MU_0, SIZE, T_INV) ;
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
~~~pythonplot Evolution of the axial curvature
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('legend', fontsize=12)
plt.rc('axes', titlesize=15) # size for the title
plt.rc('axes', labelsize=13) # size for the axes labels
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

N_min = 7
mu0 = 50 
L_0 = 230
t_inv = 15

path_invisc = "../collapsing_data/visc_data/invisc_comparison/invisc-collapse_lvl7_m50_x-ux-curv_s230_tinv15.dat"
path_visc = "../collapsing_data/visc_data/visc/visc-collapse_lvl7_m50_x-ux-curv_s230_tinv15.dat"

path_list = [path_invisc, path_visc]
labels = ["Inviscid", "Viscous"]

time_array = np.linspace(0.001, 10, 200)

#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 12
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]


# Curvature VS Time [Linear Plot] 

fig, ax = plt.subplots()

axins = inset_axes(ax, width=1.5, height=1.5, borderpad=1, loc=3, 
                   bbox_to_anchor=(.1, .45, .3125, .3125),
                   bbox_transform=ax.transAxes)
axins.patch.set_color('#fffcf4')


ax.set_xlim(27.5, 31.5)
ax.set_ylim(-5, 10)
ax.set_xlabel(r'$\tilde{t}$')
ax.set_ylabel(r'$\tilde{\kappa}$')

axins.set_xlim(25, 35)
axins.set_ylim(-1, 1)
axins.set_xlabel(r'$\tilde{t}$')
axins.set_ylabel(r'$\tilde{\kappa}$')


for k, path in enumerate(path_list):
  t, xpos, u_axis, curv = np.loadtxt(path, unpack=True)
  # Let's filter the data where the curvature is not zero:
  t_filt = t[curv != 0]
  curv_filt = curv[curv != 0]

  if k == 0:
    ax.plot(t_filt, curv_filt, lw=1.2, color = 'darkgray', 
            label = r'{}'.format(labels[k]) )
    axins.plot(t_filt[t_filt < 31.5], curv_filt[t_filt < 31.5], lw=1.2, color = 'darkgray')
  else:
    ax.plot(t_filt, curv_filt, lw=1.2, color = 'indianred', 
            label = r'{}'.format(labels[k]) )
    axins.plot(t_filt, curv_filt, lw=1.2, color = 'indianred')


threshold = 1./(L_0/np.power(2, 12))*0.5

ax.hlines(0., 27.5, 31.5, colors = 'black', linestyle = 'dashed', linewidth = 0.8)
ax.hlines(threshold, 27.5, 31.5, colors = 'purple', linestyle = 'solid', linewidth = 0.8)
axins.hlines(0., 25, 35, colors = 'black', linestyle = 'dashed', linewidth = 0.8)

ax.text(30.8, 9.2, r'$\left(2 \, \widetilde{\Delta}_{max} \right)^{-1}$', 
         color ='purple')

axins.vlines(29.82295, -1, 1, colors = 'purple', linestyle = 'dashed', linewidth = 0.6)
axins.text(29., 0.15, r'$\tilde{t}_0$', fontsize=15,
         ha='center', va='center', color="purple")

ax.legend(frameon=False, loc=3)

plt.savefig('curvature_vs_time_linear_scale_s230.svg') 
~~~

On the above figure we plot the axial curvature for the viscous simulation 
(in *red*) and the inviscid one (in *grey*). The conclusion is quite 
straightforward: both curvatures collapse onto the same curve initially, 
but rapidly deviate from each other when approaching the singularity time 
$\tilde{t}_0 \simeq 2 \, \tilde{t}_{inv} \approx 30$:

  * in the inviscid simulation, there is a curvature reversal due to the limit 
  of the spatial grid resolution, as the *divergence* of $\tilde{\kappa}$ 
  proves it, exceeding the threshold value $1/(2 \widetilde{\Delta}_{max})$: 
  beyond this limit, the simulation is no more valid as it is then under-resolved;
  * however, for the viscous simulation, a **regularizing mechanism** is triggered 
  at the curvature reversal: as it can be seen in the *inset*, 
  $\widetilde{\kappa} \in [-1 \,; 1]$ *before* and *after* the time $\tilde{t}_0$, 
  which is therefore de-singularized; to put it in other words, 
  *the singularity horizon has been crossed physically!*


To have a better judgement of the role played by viscosity close to the 
curvature reversal, we have plotted in log-log scales the time evolutions of 
the axial velocity and the maximum pressure (*see below*) in the vicinity of 
$\tilde{t}_0$. 

~~~pythonplot Axial Velocity VS Time [CAVITY COLLAPSE]
t_0_invisc = 29.75975 
t_0_visc = 29.82295 

t_0_list = [t_0_invisc, t_0_visc]

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xscale('log')
ax.set_yscale('log')

axins = inset_axes(ax, width=1.5, height=1.5, borderpad=1, loc=3, 
                   bbox_to_anchor=(.6, .55, .3125, .3125),
                   bbox_transform=ax.transAxes)
axins.patch.set_color('#fffcf4')


ax.set_xlim(1e-3, 1e1)
ax.set_ylim(1e0, 1e3)
ax.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax.set_ylabel(r'$\left|\tilde{u}_{axis} \right|$')

axins.set_xlim(25, 35)
axins.set_ylim(0, 15)
axins.set_xlabel(r'$\tilde{t}$')
axins.set_ylabel(r'$\left|\tilde{u}_{axis} \right|$')

for k, path in enumerate(path_list):
  t, xpos, u_axis, curv = np.loadtxt(path, unpack=True)
  # Let's filter the data where the curvature is not zero:
  t_filt = t[curv != 0]
  u_axis_filt = u_axis[curv != 0]

  if k == 0:
    ax.plot(t_0_list[k] - t_filt, np.abs(u_axis_filt), 
            lw=1.2, color = 'darkgray', 
            label = r'{}'.format(labels[k]) )
    axins.plot(t_filt[t_filt < 31.5], np.abs(u_axis_filt)[t_filt < 31.5], 
               lw=1.2, color = 'darkgray')
  else:
    ax.plot(t_0_list[k] - t_filt, np.abs(u_axis_filt), 
            lw=1.2, color = 'indianred', 
            label = r'{}'.format(labels[k]) )   
    axins.plot(t_filt, np.abs(u_axis_filt), lw=1.2, color = 'indianred') 


ax.plot(time_array, 5.8*np.power(time_array, -1./3.), 
        '--k', lw=0.8, label=r'$\tilde{t}^{-1/3}$')

axins.vlines(29.82295, 0, 15, colors = 'purple', linestyle = 'dashed', linewidth = 0.8)
axins.text(30.5, 7.5, r'$\tilde{t}_0$', fontsize=15,
         ha='center', va='center', color="purple")

ax.legend(frameon=False, loc=2)

plt.savefig('u-axis_vs_time_log_scale_s230.svg') 
~~~

The **axial velocities** of both simulations coincide on almost a decade, far 
from the singularity, and are sharing the same capillary-inertial self-similar 
regime. However, at one unity of viscous time from the curvature reversal, the 
curves are dividing into two routes:

  * for the viscous simulation, at $\tilde{t}_0 - \tilde{t} \simeq 1$ ($= t_{\mu}$), 
  we are reaching a *plateau* at around 
  $|\widetilde{u}_{axis}| \approx 8,5 \Rightarrow |u_{axis}| = 8,5 \, V_\mu = 612 \text{ m/s}$ 
  in water;
  * for the inviscid simulation, the self-similar regime continues on *two more 
  time decades*, before being under-resolved. 

This comparison allows us to clearly identify the time for which the **self-similarity 
is left** when viscous effects are enabled at time $t_\mu$, whereas without these 
effects the velocity diverges according to the capillary-inertial power law, 
up to the moment the simulation becomes under-resolved.


~~~pythonplot Pressure Max VS Time [CAVITY COLLAPSE]
path_invisc = "../collapsing_data/visc_data/invisc_comparison/invisc-collapse_lvl7_m50_pressure_s230_tinv15.dat"
path_visc = "../collapsing_data/visc_data/visc/visc-collapse_lvl7_m50_pressure_s230_tinv15.dat"
path_list = [path_invisc, path_visc]

t_0_invisc = 29.75975
t_0_visc = 29.82295  
t_0_list = [t_0_invisc, t_0_visc]

fig, ax = plt.subplots()
ax.set_aspect('equal')
axins = inset_axes(ax, width=1.5, height=1.5, borderpad=1, loc=3, 
                   bbox_to_anchor=(.1, .1, .3125, .3125),
                   bbox_transform=ax.transAxes)
axins.patch.set_color('#fffcf4')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 1e1)
ax.set_ylim(1e-1, 1e2)
ax.set_xlabel(r'$\tilde{t}_0 - \tilde{t}$')
ax.set_ylabel(r'$\tilde{p}_{max}$')

axins.set_xlim(25, 35)
axins.set_ylim(0, 12)
axins.set_xlabel(r'$\tilde{t}$')
axins.set_ylabel(r'$\tilde{p}_{max}$')


for k, path in enumerate(path_list):
  t, pdiff, pmax = np.loadtxt(path, unpack=True)

  if k == 0:
    ax.plot(t_0_list[k] - t, pmax, 
            lw=1.2, color = 'darkgray', 
            label = r'{}'.format(labels[k]) )
    axins.plot(t[t < 30], pmax[t< 30], 
            lw=1.2, color = 'darkgray')
  else:
    ax.plot(t_0_list[k] - t, pmax, 
        lw=1.2, color = 'indianred', 
        label = r'{}'.format(labels[k]) )
    axins.plot(t, pmax, lw=1.2, color = 'indianred') 

axins.vlines(29.82295, 0, 12, colors = 'purple', linestyle = 'dashed', linewidth = 0.8)
axins.text(30.5, 5, r'$\tilde{t}_0$', fontsize=15,
         ha='center', va='center', color="purple")

ax.plot(time_array, 2.3*np.power(time_array, -2./3.), 
        '--k', lw=0.8, label=r'$\tilde{t}^{-2/3}$')

ax.legend(frameon=False)

plt.savefig('pressure_vs_time_log_scale_s230.svg') 
~~~

Similar remarks can be said for the **maximum pressure**, since the self-similar 
regime (though only well captured in the inviscid simulation for 
$(\tilde{t}_0 - \tilde{t}) \in [10^{-2} \, ; 10^{-1}]$ due to boundary effects) 
is left much sooner in the viscous simulation. As before, the *inset* reveals 
the divergence of the maximum pressure for the inviscid case, whereas a 
de-singularization is observed when viscosity is enabled at 
$\tilde{t} - \tilde{t}_0 \lesssim t_{\mu}$. 

In the case of water in air, at the time of curvature reversal, 
$\tilde{p}_{max} \approx 7 \Rightarrow p_{max} = 362 \text{ bar}$. 
*/



/** 
We then store the interface shapes, with a smaller timestep really close to the 
singularity to focus on the potential significant morphological changes: 
*/

  event profiles (t = T_INV; t <= T_END; t += 1. ){                 
    char filename_v1[200] ;
    sprintf(filename_v1, "visc-collapse-shape_s%g_lvl%d_t%g_th%g_m%g_tinv%g.dat", 
            SIZE, MAXLEVEL, t, THETA_0*180./pi, MU_0, T_INV) ;
    FILE * fp_v1 = fopen(filename_v1, "w") ;
    output_facets (f, fp_v1) ;
    fclose (fp_v1) ;
  }

  event profiles_v2 (t = 2*T_INV-.9; t <= 2.3*T_INV; t += .1 ){                 
    char filename_v2[200] ;
    sprintf(filename_v2, "visc-collapse-shape_s%g_lvl%d_t%g_th%g_m%g_tinv%g.dat", 
            SIZE, MAXLEVEL, t, THETA_0*180./pi, MU_0, T_INV) ;
    FILE * fp_v2 = fopen(filename_v2, "w") ;
    output_facets (f, fp_v2) ;
    fclose (fp_v2) ;
  }


/** 
We now change manually the $\tilde{t}_0$ value to highlight the self-similar 
capillary inertial regime of the collapse for the inviscid simulation, as all 
the profiles are collapsing onto a master curve for different timesteps 
(*see left figure below*). On the contrary, a **clear departure** of the 
self-similar regime is observed for the viscous simulation when the simulation 
reaches one unity of viscous time before the presumed time of singularity 
(*right figure below*). 
Viscous effects are thus clearly responsible for this behaviour.

~~~pythonplot Pre-Singularity Capillary-Inertial Rescaling
N_min = 7
mu0 = 50 
L_0 = 230
t_inv = 15
lvl = 12
th = 120

t_0 = 30.1 

time_1 = np.arange(15, 30, 1)
time_str_1 = np.array(["{}".format(t) for t in time_1])
time_shift_1 =  t_0 - time_1
time_scaled_1 = np.power(time_shift_1, -2./3.)

time_2 = np.arange(29.1, 29.8, 0.1)
time_str_2 = np.array(["{:.1f}".format(t) for t in time_2])
time_shift_2 =  t_0 - time_2 
time_scaled_2 = np.power(time_shift_2, -2./3.)


#                                   PLOTS
# ------------------------------------------------------------------------------
Ncolors = 25
colormap = plt.cm.viridis# LinearSegmentedColormap
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
ax.set_aspect('equal')

ax.set_xlim(-50,20)
ax.set_ylim(0,70)

ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')

for k, t in enumerate(time_str_1):
  path_shape = "../collapsing_data/visc_data/invisc_comparison/shape/invisc-zenon-shape_s{}_".format(L_0) 
  if k < 10:
    facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(10, t, th, mu0, t_inv), unpack=False)
  else:
    facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(12, t, th, mu0, t_inv), unpack=False)
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments):
    if (k == 0) and (l == 0):
      ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = mapcolors[k], 
              label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1))) 
    elif (k == len(time_str_1) - 2):
      if (l == 0):
        ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
                '-', lw=1.2, color = 'darkred', 
                label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1)))
      else:
        ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = 'darkred') 
    elif (k == len(time_str_1) - 1):
      if (l == 0):
        ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
                '-', lw=1.2, color = 'purple', 
                label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1)))
      else:
        ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = 'purple') 
    else:
      ax.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
        '-', lw=1.2, color = mapcolors[k]) 

for k, t in enumerate(time_str_2):
  path_shape = "../collapsing_data/visc_data/invisc_comparison/shape/invisc-zenon-shape_s{}_".format(L_0) 
  facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(12, t, th, mu0, t_inv), unpack=False)
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments): 
    if (k == len(time_str_2) - 1) and (l == 0):
      ax.plot(segment[:, 0]*time_scaled_2[k], segment[:, 1]*time_scaled_2[k],
              '-', lw=1.2, color = mapcolors[len(time_str_1) + k], 
              label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1))) 
    else:
      ax.plot(segment[:, 0]*time_scaled_2[k], segment[:, 1]*time_scaled_2[k],
        '-', lw=1.2, color = mapcolors[len(time_str_1) + k]) 

ax.legend(frameon=False, loc=1)


# ------------------
# Viscous Simulation
# ------------------

time_1 = np.arange(15, 30, 1)
time_str_1 = np.array(["{:.1f}".format(t) for t in time_1])
time_shift_1 =  t_0 - time_1
time_scaled_1 = np.power(time_shift_1, -2./3.)

ax2.set_aspect('equal')

ax2.set_xlim(-50,20)
ax2.set_ylim(0,70)

ax2.set_xlabel(r'$\xi$')
ax2.set_ylabel(r'$\eta$')

for k, t in enumerate(time_str_1):
  path_shape = "../collapsing_data/visc_data/visc/shape/visc-zenon_gas_refined-shape_s{}_".format(L_0) 
  if k < 10:
    facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(10, t, th, mu0, t_inv), unpack=False)
  else:
    facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(12, t, th, mu0, t_inv), unpack=False)
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments):
    if (k == 0) and (l == 0):
      ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = mapcolors[k], 
              label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1))) 
    elif (k == len(time_str_1) - 2):
      if (l == 0):
        ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
                '-', lw=1.2, color = 'darkred', 
                label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1)))
      else:
        ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = 'darkred') 
    elif (k == len(time_str_1) - 1):
      if (l == 0):
        ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
                '-', lw=1.2, color = 'purple', 
                label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1)))
      else:
        ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
              '-', lw=1.2, color = 'purple') 
    else:
      ax2.plot(segment[:, 0]*time_scaled_1[k], segment[:, 1]*time_scaled_1[k],
        '-', lw=1.2, color = mapcolors[k]) 


for k, t in enumerate(time_str_2):
  path_shape = "../collapsing_data/visc_data/visc/shape/visc-zenon_gas_refined-shape_s{}_".format(L_0) 
  facets= np.loadtxt (path_shape + "lvl{}_t{}_th{}_m{}_tinv{}.dat".format(12, t, th, mu0, t_inv), unpack=False)
  N_seg = int (0.5*facets.shape[0])
  segments = np.split (facets, indices_or_sections=N_seg)
  for l, segment in enumerate(segments): 
    if (k == len(time_str_2) - 1) and (l == 0):
      ax2.plot(segment[:, 0]*time_scaled_2[k], segment[:, 1]*time_scaled_2[k],
              '-', lw=1.2, color = mapcolors[len(time_str_1) + k], 
              label = r'$\tilde{t}_0 - \tilde{t} = \,\,$' + r'${}$'.format(np.round(t_0 - float(t), 1))) 
    else:
      ax2.plot(segment[:, 0]*time_scaled_2[k], segment[:, 1]*time_scaled_2[k],
        '-', lw=1.2, color = mapcolors[len(time_str_1) + k]) 



# Time arrow:
ax2.annotate(r' ', xy=(-1.9, 29), 
            xytext=(-12, 20),
            arrowprops=dict(arrowstyle="->", 
                          color = "grey"
                          # width=0.005, 
                          # head_width=6*0.005,
                          ),
            color = "grey",
            fontsize=12
            )
ax2.text(.5, 27.5, r'$\tilde{t}$', fontsize=15,
        ha='center', va='center', color="grey")

ax2.legend(frameon=False, loc=1)


plt.savefig('selfsim_shape_pre_sing_invisc+visc_s230.svg') 
~~~
*/







  event vel_p_maps (t = 2*T_INV-.9; t <= 2.3*T_INV; t += .1 )
  {
    char fileup[200];
    char fileinterf[200];
    sprintf (fileup, "vel_p_dissip_maps_t%g_visc_collapse.dat", t);
    sprintf (fileinterf, "interf_maps_t%g_visc_collapse.dat", t);
    FILE * fup = fopen (fileup, "w");
    FILE * finterf = fopen (fileinterf, "w");
    position (f, x_interf, {1,0}, add=false);
    position (f, y_interf, {0,1}, add=false);
    foreach (){
      fprintf(fup, "%g %g %g %g %g %g\n",
                    x, y, u.x[], u.y[], p[], visc_dissip[]);
      if ( interfacial (point, f) ) {
        if (f[] > 0.1 && f[] < 0.9){
          fprintf(finterf, "%g %g\n", x_interf[], y_interf[]);
        }
      }
    } 
    fclose(fup);
    fclose(finterf);
  }

  event vel_p_maps_v2 (t = T_INV; t <= T_END; t += 1 )
  {
    char fileup_v2[200];
    char fileinterf_v2[200];
    sprintf (fileup_v2, "vel_p_dissip_maps_t%g_visc_collapse.dat", t);
    sprintf (fileinterf_v2, "interf_maps_t%g_visc_collapse.dat", t);
    FILE * fup_v2 = fopen (fileup_v2, "w");
    FILE * finterf_v2 = fopen (fileinterf_v2, "w");
    position (f, x_interf, {1,0}, add=false);
    position (f, y_interf, {0,1}, add=false);
    foreach (){
      fprintf(fup_v2, "%g %g %g %g %g %g\n",
                    x, y, u.x[], u.y[], p[], visc_dissip[]);
      if ( interfacial (point, f) ) {
        if (f[] > 0.1 && f[] < 0.9){
          fprintf(finterf_v2, "%g %g\n", x_interf[], y_interf[]);
        }
      }
    } 
    fclose(fup_v2);
    fclose(finterf_v2);
  }



  #if GAS_REFINEMENT
    event logfile (i++) {
      coord ubar;
      foreach_dimension() {
        stats s = statsf(u.x);
        ubar.x = s.sum/s.volume;
      }

      double ke = 0., vd = 0., vol = 0.;
      foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
        vol += dv();
        foreach_dimension() {
          // mean fluctuating kinetic energy
          ke += dv()*sq(u.x[] - ubar.x);
          // viscous dissipation
          vd += dv()*(2.*(
            sq(u.x[1] - u.x[-1]) + sq(u.x[0,1] - u.x[0,-1])
            + sq(u.y[1] - u.y[-1]) + sq(u.y[0,1] - u.y[0,-1])
            + 2.*(u.y[1] - u.y[-1])*(u.x[0,1] - u.x[0,-1])
            )/sq(2.*Delta)
          + 2.*sq(u.y[]/y)
          )*(f[]*MU_LIQ + (1. - f[])*MU_GAS);;
        }
      }
      ke /= 2.*vol;
      vd *= MU_LIQ/vol;
      static FILE * fd = fopen("stats_visc_collapse.dat","a");
      if (i == 0) {
        fprintf (fd, "t dt dissip. energy \n");
      }
      fprintf (fd, "%g %g %g %g\n", t, dt, vd, ke);

      fflush(fd); // empty buffer
    }
  #endif

  /** The following event is a clone of one of Aliénor Rivière's functions in
  [`stagnation.h`](http://basilisk.fr/sandbox/ariviere/stagnation.h) and adapted 
  for storing values in parallel computations. 
  It is used for post-treating these data directly with *Python*:
  */

  event getInterface(i++){
    char fixname[100];//x
    char fiyname[100];//y
    char fikname[100];//curvature
    char fik2Dname[100];//2D curvature
    char fixpname[100];//x pressure
    char fiypname[100];//y pressure
    char fipname[100];//pressure
    char fizname[100];//z-position
    char figname[100];//long. vel. gradient
    sprintf(fixname, "fix_%d.dat", pid() );
    sprintf(fiyname, "fiy_%d.dat", pid() );
    sprintf(fikname, "fik_%d.dat", pid() );
    sprintf(fik2Dname, "fik2D_%d.dat", pid() );
    sprintf(fixpname, "fixp_%d.dat", pid() );
    sprintf(fiypname, "fiyp_%d.dat", pid() );
    sprintf(fipname, "fip_%d.dat", pid() );
    sprintf(fizname, "fiz_%d.dat", pid() );
    sprintf(figname, "fig_%d.dat", pid() );
    FILE * fix = fopen(fixname, "a");
    FILE * fiy = fopen(fiyname, "a");
    FILE * fik = fopen(fikname, "a");
    FILE * fik2D = fopen(fik2Dname, "a");
    FILE * fixp = fopen(fixpname, "a");
    FILE * fiyp = fopen(fiypname, "a");
    FILE * fip = fopen(fipname, "a");
    FILE * fiz = fopen(fizname, "a");
    FILE * fig = fopen(figname, "a");
    fprintf (fix, "%g ", t);
    fprintf (fiy, "%g ", t);
    fprintf (fik, "%g ", t);
    fprintf (fik2D, "%g ", t);
    fprintf (fixp, "%g ", t);
    fprintf (fiyp, "%g ", t);
    fprintf (fip, "%g ", t);
    fprintf (fiz, "%g ", t);
    fprintf (fig, "%g ", t);

    scalar curve[]; //curvature
    scalar curve2D[]; //2D-curvature
    curvature (f, curve);
    boundary ((scalar*){curve});
    curvature_2D (f, curve2D);
    boundary ((scalar*){curve2D});

    //position with the height function
    scalar xpos[];
    scalar ypos[];

    position(f, xpos, {1,0});
    position(f, ypos, {0,1});

    foreach(){
      if (interfacial(point,f)){
        if (f[] > 0.1 && f[] < 0.9){
          double xp = xpos[] - 3.*Delta*n_front.x[];
          double yp = ypos[] - 3.*Delta*n_front.y[];
          fprintf (fix, "%g ", xpos[]);
          fprintf (fiy, "%g ", ypos[]);
          fprintf (fik, "%g ", curve[]);
          fprintf (fik2D, "%g ", curve2D[]);
          fprintf (fixp, "%g ", xp);
          fprintf (fiyp, "%g ", yp);
          fprintf (fip, "%g ", interpolate(p, xp, yp));
        }
      }
    }
    // Computation of the longitudinal velocity gradient:
    foreach() {
      if (y <= Delta) {
        fprintf (fiz, "%g ", x);
        duz_dz[] = ( u.x[1] - u.x[-1] )/( 2.*Delta );
        fprintf (fig, "%g ", duz_dz[]);
      }
    }
    fprintf (fix, "\n");
    fprintf (fiy, "\n");
    fprintf (fik, "\n");
    fprintf (fik2D, "\n");
    fprintf (fixp, "\n");
    fprintf (fiyp, "\n");
    fprintf (fip, "\n");
    fprintf (fiz, "\n");
    fprintf (fig, "\n");
    fclose(fix);
    fclose(fiy);
    fclose(fik);
    fclose(fik2D);
    fclose(fixp);
    fclose(fiyp);
    fclose(fip);
    fclose(fiz);
    fclose(fig);
    return 0;
  }
#endif



/** 
### Movie
*/

#if MOVIE
  #include "view.h"
  #include "cmaps.h"

  event mov_p_u (t+=.01) {
    char legend[100];
    sprintf(legend, "t = %0.2g", t);
    view (tx = -0.05, ty = -0.01, width = 1000., height = 600., fov = 8.);
    box();
    squares ("p", spread=-1, linear=true, map=cool_warm); 
    draw_vof ("f", lw = 1.5);
    mirror ({0.,1}) {
      draw_vof ("f", lw = 1.5);
      squares("u_rs", spread=-1, linear=true, map=viridis);
      cells ();
    }
    draw_string(legend, 1, size = 20., lw = 2.1);
    save ("visc-collapse_s230_N7-9_tinv15.mp4");
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

