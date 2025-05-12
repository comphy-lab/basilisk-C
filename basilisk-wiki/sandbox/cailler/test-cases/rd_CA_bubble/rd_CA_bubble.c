/** 
# Axisymmetric Bubble Growth and Kinematic Reversal with Varying Contact Angles 

## Introduction

A bubble in expansion (due to a source placed at the origin) is initialized with 
analytical contact angles at left and bottom border (therefore no longer at the 
axis of symmetry due to the Y-AXIS offset we are imposing).
As it grows, the bubble will have varying contact angles, and will even 
reach the right border.
Hence, we need to tell to `Basilisk` how the volume fractions in the ghost cells 
are updated, in order to be consistent with the advection of the cavity 
interface, by applying the algorithm detailed 
[here](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/f_BC_vel_tan.h) 

Going a step further, we demonstrate the *kinematic* ***invariance*** 
of such an inviscid simulation, by simply applying at a given event a 
velocity reversal.

This test case is then a generalization of the [simple test of a 2D oblique 
interface placed in a tangent flow](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_corrected.c) 
to non-uniform 3D-AXI flows with random contact angles 
(the most general case possible).



## Results

![Kinematic Invariance of an AXI-bubble growth](rd_CA_bubble/reverse_bubble_growth_AXI_.mp4)(width="400" height="800")


## Code
*/

#define LEVEL 6
#define MAXLEVEL 7
#define SIZE 1.e2
#define T_INV 2.e2
#define T_END 2.*T_INV
#define R0 .6*SIZE
#define Y_OFFSET .5*SIZE 
#define X_OFFSET -.2*SIZE // .7*SIZE with Y_OFFSET =.5*SIZE for bottom right loc.

#include "axi.h"
#include "centered_bubble.h" // removing lines uf.n[bottom] = 0.;  
                             // and uf.t[bottom] = dirichlet(0); 
                             // as we do not have an axis in the simulation
#include "contact.h"
#include "two-phase_bubble.h"
#include "tension.h"
#include "u_BC_bubble.h" // analytical functions for the velocity field
#include "f_BC_bubble.h"


// Vectors, scalar fields...
scalar omega[];
scalar f0[];
scalar curv_viz[]; // to store curvature 
scalar alpha_front[];
vector h[];
vector n_front[];

/**
### Boundary Conditions

We need to compute correctly for both $\mathbf{n}$ (normal vector to the 
interface) and $\alpha$ (related intercept) their direction and 
value due to the presence of fluxes on top and left boundaries. 
Indeed, they are set by default to *symmetry conditions*, 
inducing errors when complex inflows/outflows exist, as demonstrated in the 
[advection of an oblique interface in a tangent flow](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/vel_tan_no_correction.c). 

For a bubble in expansion, the normal (directed from the liquid 
phase towards the gas one) reads in axisymmetric coordinates:
$$
n_z = - \dfrac{\cos \theta}{|\cos \theta| + |\sin \theta|}
\quad ; \quad
n_r = - \dfrac{\sin \theta}{|\cos \theta| + |\sin \theta|}
$$
where $\theta = \arctan(r/z)$ and using the $L_1-$norm convention.
*/



              /* Normal vectors at the interface (ghost cells) */

// pay attention to the "foreach_dimension()" declaration for the normal components
n_front.y[bottom] = -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
n_front.x[bottom] = -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
alpha_front[bottom] = plane_alpha (f[ghost], (coord){
  -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x)))),
  -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))))
  });


n_front.x[left] = -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
n_front.y[left] = -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
alpha_front[left] = plane_alpha (f[ghost], (coord){
  -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x)))),
  -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))))
  });


// pay attention to the "foreach_dimension()" declaration for the normal components
n_front.y[top] = -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
n_front.x[top] = -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
alpha_front[top] = plane_alpha (f[ghost], (coord){
  -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x)))),
  -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))))
  });


n_front.x[right] = -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
n_front.y[right] = -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))));
alpha_front[right] = plane_alpha (f[ghost], (coord){
  -cos(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x)))),
  -sin(atan2(y,x))/(fabs(-cos(atan2(y,x))) + fabs(-sin(atan2(y,x))))
  });


/**
Then, we pass to the volume fraction field of ghost cells an auxiliary one 
that stores the correct computation of volume fraction, thanks to the 
applied algorithm from `f_BC_bubble.h`:
*/

// Volume Fractions:
f[left] = f0[-1,0];
f[bottom] = f0[0,-1];
f[right] = f0[1,0];
f[top] = f0[0,1];

/**
After that, we apply the proper contact angles, found analytically, 
but varying along with the expansion/collapse of the bubble, depending on 
the contact point of $(r,z)$ coordinates:
*/
  // Random Contact Angles (but analytically determined):
h.t[left] = contact_angle (pi - atan2(y,x)) ;
h.t[bottom] = contact_angle (pi/2. + atan2(y,x)) ; // if bottom right contact
// h.t[bottom] = contact_angle ((3.*pi)/2. - atan2(y,x)) ; // if bottom left contact

h.t[right] = contact_angle (atan2(y,x)) ;
h.t[top] = contact_angle (atan2(y,x) - pi/2.) ;

/**
Pressure and velocities BCs have to be set properly considering 
the inflow/outflow created, to avoid the divergence of the *Poisson solver*:
*/

  // Inflow Conditions: 
p[left] = neumann (0.);
pf[left] = neumann (0.);

p[bottom] = neumann (0.);
pf[bottom] = neumann (0.);

  // Outflow Conditions:
u.n[top] = neumann (0.); // to disable if we choose having dirichlet BC for top

p[top] = dirichlet (0.); // should be neumann BC if u.n/t[top] are dirichlet BC
pf[top] = dirichlet (0.);

p[right] = neumann (0.);
pf[right] = neumann (0.);

/**
For the *kinematic reversal* testing the kinematic **invariance**, 
please refer to the BCs defined at the section 
[*Kinematic Invariance*](http://basilisk.fr/sandbox/cailler/test-cases/rd_CA_bubble/rd_CA_bubble.c#kinematic-invariance).

### Generic Events
*/

int main() {
  size(SIZE) ;
  init_grid(1 << LEVEL) ;
  origin(X_OFFSET, Y_OFFSET) ;

  rho1 = 1., rho2 = 1.e-3 ; // here 1/ is the liquid and 2/ is the gas
  f.height = h ;
  f.sigma = 1.; // enable surface tension
  mu1 = 0, mu2 = 0;

  run() ;
}


event init (t = 0){
  // Initial free-surface definition:
  #if TREE
    do{
    fraction (f0, (sq(x) + sq(y) - sq(R0))); 
    } while (adapt_wavelet ({f0}, (double[]){1e-2}, MAXLEVEL).nf); 
  #else // !TREE
    fraction (f0, (sq(x) + sq(y) - sq(R0))); 
  #endif

/**
Correct initialization of contact angles at the intersected borders:
*/
  foreach_boundary(bottom)
    f0[0,-1] = f_BC_bottom (x, y, Delta, f0[], f0[1,0], f0[-1,0]);
  foreach_boundary(left)
    f0[-1,0] = f_BC_left (x, y, Delta, f0[], f0[0,1], f0[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (x, y, Delta, f0[], f0[1,0], f0[-1,0]);
  foreach_boundary(right)
    f0[1,0] = f_BC_right (x, y, Delta, f0[], f0[0,1], f0[0,-1]);

  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0}); // for boundary conditions on levels

/**
At first, we need to initialize the volume fraction field for ghost cells with 
the auxiliary one. The velocity field is then defined, and we control the 
update of *height functions* for minimizing curvature computation errors.
*/
  foreach(){
    f[] = f0[];

    u.x[] =  uz_A(x,y); 
    u.y[] =  ur_A(x,y);
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
  #if AXI
  foreach()
    omega[] = (u.y[1,0]-u.y[-1,0]+u.x[0,-1]-u.x[0,1]) / (2.*Delta + SEPS);
  #else // !AXI
    vorticity (u, omega);
  #endif
  curvature (f, curv_viz, 1, add = false);
}


#if TREE
  event adapt (i++) {
    adapt_wavelet ({f,u}, (double[]){1e-2,1e-2,1e-2}, MAXLEVEL);

/**
To ensure stability and correct updating of volume fractions each time
AMR is applied, in the right order, we need to update again volume 
fraction values at boundaries:
*/

  foreach_boundary(bottom)
    f0[0,-1] = f_BC_bottom (x, y, Delta, f[], f[1,0], f[-1,0]);
  foreach_boundary(left)
    f0[-1,0] = f_BC_left (x, y, Delta, f[], f[0,1], f[0,-1]);
  foreach_boundary(top)
    f0[0,1] = f_BC_top (x, y, Delta, f[], f[1,0], f[-1,0]);
  foreach_boundary(right)
    f0[1,0] = f_BC_right (x, y, Delta, f[], f[0,1], f[0,-1]);

  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0}); // for boundary conditions on levels

  boundary({f});
  if (f.height.x.i){
    heights (f, f.height);
  }
  boundary({f});
  }
#endif

/**
### Kinematic Invariance 
*/
event reversal_BC (i++) {
  if (t < T_INV){
    u.n[left] = dirichlet ( uz_A(x,y) );
    u.t[left] = dirichlet ( ur_A(x,y) );
    u.n[right] = dirichlet ( uz_A(x,y) );
    u.t[right] = dirichlet ( ur_A(x,y) );

    u.n[bottom] = dirichlet ( ur_A(x,y) );
    u.t[bottom] = dirichlet ( uz_A(x,y) );
  }
  else {
    u.n[left] = dirichlet (- uz_A(x,y) );
    u.t[left] = dirichlet (- ur_A(x,y) );
    u.n[right] = dirichlet (- uz_A(x,y) );
    u.t[right] = dirichlet (- ur_A(x,y) );

    u.n[bottom] = dirichlet (- ur_A(x,y) );
    u.t[bottom] = dirichlet (- uz_A(x,y) );
  }
}

event reversal (t = T_INV) {
  foreach_face(){
    uf.x[] = 0.;
  }
  foreach(){
    u.x[] *= -1.;
    u.y[] *= -1.;
  }
} 

event end (t = T_END){}


/**
### Outputs 
*/

event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffer
}

#include "view.h"
#include "cmap_black-body.h"


event mov (t += 1.) {
  clear();
  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  view (fov = 60, tx = -0.28, width = 800, height = 1600);
  draw_vof ("f", lw = 1.5);
  squares ("p", spread=-1, linear=true, map=black_body);
  vectors ("u", scale = 0.3);
  mirror ({0.,1}) {
    draw_vof ("f", lw = 1.5);
    squares ("omega", spread=-1, linear=true, map=cool_warm);
    cells();
  }
  draw_string(legend, 1, size = 20., lw = 2.1);
  save ("reverse_bubble_growth_AXI_.mp4");
}