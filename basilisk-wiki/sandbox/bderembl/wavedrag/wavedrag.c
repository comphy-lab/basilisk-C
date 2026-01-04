/**
# Topographic waves

Implementation of the low-topography interaction described in Lott et
al. (2024).

The mountain is defined as

$$  h(x)=H \exp(-\frac{x^2}{2L^2}) $$

with $H$ the height of the mountain and $L$ the characteristic lateral scale.
We introduce the non dimension slope

$$S = \frac{H}{L} $$

To introduce the incoming flow, we first define

$$U_V(Z)=\frac{u_*}{\kappa} \log\left(\frac{\sinh\kappa(Z+z_0)/l_\infty}{\sinh\kappa z_0/l_{\infty}}\right)$$

where $Z$ is the terrain following vertical coordinate, $l_{\infty}$ is a mixing length, far from the topography, $\kappa$ the Von Karman constant, and $z_0$ a roughness length, and $u_*$ a velocity scale.


This velocity profile is the analytical solution of the log layer with specific
conditions see Lott et al.). However, since this profile diverges at $Z
\rightarrow \infty$ we introduce a boundary layer depth $d$ and define the incoming profile as

$$U(Z)=\frac{u_{{}*{}}d}{l_\infty} \tanh\left[ \frac{l_\infty}{u_{{}*{}} d } U_V(Z) \right]$$


The background buoyancy profile is defined as

$$B(Z)=\frac{b_*}{\kappa} \log\left(\frac{\sinh\kappa(Z+z_0)/l_{\infty}}{\sinh\kappa z_0/l_{\infty}}\right)$$

with $b_* = N^2 l_\infty$ a buoyancy scale.

The vertical viscosity/diffusivity is defined as

$$ \nu=l^2\left\|\frac{\partial u}{\partial Z}\right\|$$

with $l$ the mixing length:

$$l=l_{\infty} \tanh \left( \kappa\frac{Z+z_0}{l_{\infty}} \right)$$

In order to damp wave reflections at the top and right edge of the domain, we
use a sponge layer.


The code below is heavily influenced by Antoon's lee wave code, adapted to the
present flow and topography profile.

TODO: this should now use the layer formulation.


In typical conditions, we need to simulate a domain with a very small aspect
ratio Ly/Lx << 1 in order to see the downtream evolution of topographic
waves. This can be very easily done with the layer formulation (switch to
layers!). With the centered version of the Navier-Stokes solver, we need to run
this code in parallel. Here is a possible set of parameter

with 144 procs

npx = 12
npy = 3
N = 2048

I get an aspect ratio of 1/4 (npy/npx) the number of grid points for this value
of N is 3072x768 (each tile is 128x128)




### Compilation

local:

~~~bash
CC99='mpicc -std=c99' qcc -disable-dimensions -D_MPI=1 -lm -lpnetcdf -O3 wavedrag.c
mpirun -np 144 a.out 
~~~

Cluster:

~~~bash
qcc -D_MPI=1 -grid=multigrid -disable-dimensions -source wavedrag.c
mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_MPI=1 -O3 _wavedrag.c -lm -lpnetcdf
nohup mpirun -machinefile ~/MachineFile -np 144 a.out &
~~~
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#include "../libs/extra.h"
#include "../libs/pnetcdf_bas.h"

/**
A buoyancy tracer field is declared: 
*/
scalar b[];
scalar * tracers = {b};

/**
Here are the default parameters
*/

int N0 = 256;    // nominal number of grid point
int npx = 4;     // number of processor in the horizontal direction
int npy = 1;     // number of processor in the vertical direction
double S = 1;    // mountain slope
double N2 = 1;   // background stratification
double Re = 1e4; // inverse of background viscosity
double Pr = 1;   // Prandtl number
double tend = 1; // end of the simulation
double dtout = 1; // output frequency

double tau_sponge = 0;
double itau_sponge = 0;

double bg_u = 0;  // experiemntal: background u
double sh_u = 0;  // experimental: background shear
double D = 0;    // boundary layer depth
double xm = 0; // position of the mountain
double xs = 0; // x-position of the sponge
double ys = 0; // y-position of the sponge
int slip_bc = 0; //0: no slip, 1: free slip
int buoyancy_bc = 0; // 0: neumann, 1: dirichlet


// viscosity arrays
face vector av[];
face vector u0[];
face vector muc[];
scalar f[];

mgstats mgb;
face vector kappac[];

// log profile
double kappa_vk = 0.4;
double l_mix = 0.;
double z0 = 0.;
double z_bot;
double eps = 1e-7; // a small number to avoid division by zero
//#define l_mix0 ((kappa_vk*(y + z0)))


char* fileout = "vars.nc";


/**
## Definition of the topography, mean flow, stratification and mixing length (see definition above)

We put the topography slightly above the bottommost cell so that all cells are
embeded
*/

// offset so that the topography is exactly at a grid cell (even for refined cells) (z_bot)
#define MOUNTAIN (S*exp(-(sq(x - xm)/2.)) - y + 1.4987654321*L0/N)

#define SHEAR  (MOUNTAIN < 0?(bg_u + sh_u*y + D*tanh(l_mix/(kappa_vk*(D+eps))*log(sinh(kappa_vk*(y + z0)/(l_mix+eps))/sinh(kappa_vk*z0/(l_mix+eps)))  )):0)

#define STRATIFICATION  (MOUNTAIN < 0?(N2*l_mix/kappa_vk * log(sinh(kappa_vk*(y - 1.4987654321*L0/N+ z0)/(l_mix+eps))/sinh(kappa_vk*z0/(l_mix+eps))) ):0)

#define l_mix0 ((l_mix*tanh(kappa_vk*(- MOUNTAIN + z0)/(l_mix+eps))))


int main(int argc,char* argv[]) {

  params = array_new();

  add_param("N"    , &N0    , "int");
  add_param("npx"  , &npx   , "int");
  add_param("npy"  , &npy   , "int");
  add_param("L0"   , &L0    , "double");
  add_param("Re"   , &Re    , "double");
  add_param("Pr"   , &Pr    , "double");
  add_param("N2"   , &N2    , "double");
  add_param("S"    , &S     , "double");
  add_param("D"    , &D     , "double");
  add_param("l_mix", &l_mix , "double");
  add_param("z0"   , &z0    , "double");
  add_param("xm"   , &xm    , "double");
  add_param("xs"   , &xs    , "double");
  add_param("ys"   , &ys    , "double");
  add_param("tend" , &tend  , "double");
  add_param("dtout", &dtout , "double");
  add_param("bg_u" , &bg_u  , "double");
  add_param("sh_u" , &sh_u  , "double");
  add_param("tau_sponge", &tau_sponge , "double");
  
  add_param("slip_bc"     , &slip_bc , "int");
  add_param("buoyancy_bc" , &buoyancy_bc , "int");
  
    // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(file_param);

  N = N0;

  if (tau_sponge != 0) itau_sponge = 1/tau_sponge;

  /**
     Viscosity CFL = 0.5
   */
  if (Re  != 0) DT = 0.5*min(DT,sq(L0/N)*Re/4.);

/**
   Create an output directory named outdir_0001 or use a higher number if 0001
   already exists
*/

  create_outdir();

/**
   copy the configuration file in that directory
*/
  backup_config(file_param);

  periodic (left);
  a = av;
  mu = muc;
  dimensions(nx=npx, ny=npy);
  init_grid (N); 

  // I don't know why I get a segfault without this
  origin (-1e-7, 0);

  TOLERANCE = 1e-5; // low value (1e-5) for reliable energy budgets


  run();

/**
   Free the params variable
*/
  array_free (params);

}


/**
 We adjust the viscosity and diffusion coefficient according to the mixing
 length definition.
*/
event properties (i++) {
  /* foreach_face() */
  /*   muc.x[] = fm.x[]/Re; */
  foreach_face(x)
    muc.x[] = fm.x[]/Re;
  foreach_face(y)
    muc.y[] = fm.y[]/Re + fm.y[]*sq(l_mix0)*fabs(face_gradient_y (u.x, 0));
  boundary((scalar*){muc});
}

event tracer_diffusion (i++) {

  /* foreach_face() */
  /*   kappac.x[] = fm.x[]/(Pr*Re); */
  foreach_face(x)
    kappac.x[] = fm.x[]/(Pr*Re);
  foreach_face(y)
    kappac.y[] = fm.y[]/(Pr*Re) + fm.y[]*sq(l_mix0)*fabs(face_gradient_y (u.x, 0));

  boundary((scalar*){kappac});

  mgb = diffusion (b, dt, kappac);
  boundary ({b});
}



/**
 Specify topography and boundary conditions
*/

event init (t = 0) {

  solid (cs, fs, -MOUNTAIN);

  b[top] = dirichlet (STRATIFICATION);
//  b[top] = dirichlet(100.0); // verification von karman

  if (buoyancy_bc == 0){
    b[bottom] = dirichlet(0.);
    b[embed]  = dirichlet(0.);
  } else{
    b[bottom] = neumann(0.);
    b[embed]  = neumann(0.);
  }


  u.t[top] = neumann(0.0);
  u.n[top] = dirichlet(0.); // no flux
//  u.t[top] = dirichlet(1.0); // verification von karman

  if (slip_bc == 0){// no slip
    u.t[bottom] = dirichlet(0.);
    u.t[embed] = dirichlet(0.);
    u.n[embed] = dirichlet(0.); // see discussion here https://groups.google.com/g/basilisk-fr/c/NxPmr50E0KM
  }
  else{ // free slip
    u.t[bottom] = neumann(0.);
    u.t[embed] = neumann(0.);
    u.n[embed] = neumann(0.); // see discussion here https://groups.google.com/g/basilisk-fr/c/NxPmr50E0KM
  }

/**
Initial conditions
*/
  foreach() {
    u.x[] = SHEAR*cs[];
    b[] = STRATIFICATION*cs[];
  }

  
  // Store initial condition for sponge relaxation
  foreach_face()
    u0.x[] = u.x[];
  
  // create netcdf output
  char file_tmp[100];
  sprintf (file_tmp,"%s%s", dpath, fileout);
  create_nc({b,u,p}, file_tmp);  

}

/**
Add the acceleration due to boyancy on the z coordinate
*/ 

event acceleration (i++) {
  coord grav = {0, 1};
  foreach_face()
    av.x[] = grav.x*face_value(b, 0);
  /* foreach_face(y) */
  /*   av.y[] = fm.y[]*(b[] + b[0, -1])/2.; */
  /**
An inlet profile is mimicked by forcing the velocities and buoyancy
field at the left-hand-side of the domain.
  */

  foreach() {
    if ( x < X0 + xs || y > 0.8*L0*Dimensions.y/Dimensions.x) {
      u.x[] += fm.x[]*dt*(SHEAR - u.x[])*itau_sponge;
      u.y[] += fm.y[]*dt*(      - u.y[])*itau_sponge;
      b[] += cm[]*dt*(STRATIFICATION - b[])*itau_sponge;
    }
  }
}



event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event output (t = 0; t <= tend;  t += dtout) {
  write_nc();
}

/**
## Reference

~~~bib
@article{lott2024,
author = {Lott, Francois and Beljaars, Anton and Pauget, Lucile and Deremble, Bruno},
title = {Neutral and stratified turbulent boundary-layer flow over low mountains},
journal = {Quarterly Journal of the Royal Meteorological Society},
volume = {150},
number = {758},
pages = {195-212},
keywords = {flow blocking, mountain drag, mountain waves, neutral and stratified boundary layers, Reynolds stress, sheltering effect},
doi = {https://doi.org/10.1002/qj.4591},
url = {https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.4591},
eprint = {https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/qj.4591},
year = {2024}
}
~~~
*/ 

