/**
# Hybrid Level-set/Embedded boundary method

We suppose here that the embed.h has been previously included. The motion of the
interface is only driven by the diffusion of two scalars TL and TS
*/

/**
## Mesh adaptation

We need to define a secondary set of volume and face fractions that will be used
for restricting/prolongating the variables associated with both phases
simultaneously.
*/
#include "double_embed-tree.h"
/**
## Level-set functions and variables

We need a set of variables for the level-set method to be used, we give here an
exhaustive list:

- a level-set function **`dist`**
- two scalars for the temperature in each phase (two-fluid method) **`TL`** and **`TS`**
- two phase change velocity fields: **`vpc`** calculated on the centroids, **`vpcf`**
cell-centered and continuous for the advection of the level-set function
- two pointers that behave like tracers
- a face metric for the fluxes **`muv`**
- a maximum timestep specific to the level-set method
- number of iteration for the reconstruction of a continuous velocity field
`itrecons`
- two variables linked to the Narrow Band approach we have taken **`nb_cell_NB`**
the number of cells in each direction from the 0-level-set where the signed
distance function is correctly calculated, **`NB_width`** which is half of the
thickness of the narrow band

For the application of physical laws we define several parameters:

- the latent heat
- the thermal capacity of each phase **`lambda`**
- two constant for the Gibbs-Thomson law **`epsK`** and **`epsV`**
- one constant for a fourfold-anisotropy **`eps4`**, the user can define its own
anisotropy if needed, by defining a global **`ANISO`** constant
- an equilibrium temperature **`T_eq`**
- a curvature field

We also define various helpful functions:

- `invertcs`: changes a metric to its complementary to 1
- interpolate a cell-centered field, possible interpolations are: linear,
quadratic (default), cubic
- a face value operator that ignores the metric

*/
#include "level_set.h"

/**
## Build a continuous phase change velocity field

Also calculate a timestep accordingly if **`dtLS`** is defined
*/
#include "LS_speed.h"

/**
## Advect this continuous velocity field
*/

#include "LS_advection.h"

/**
## Test
*/
#include "alex_functions2.h"

/**
The advection of the level-set function must be done before the tracer diffusion
*/

event tracer_diffusion(i++){
  advection_LS(
  dist,
  cs,fs,
  TS,TL,
  vpcf,
  itredist = 10,
  s_clean = 1.e-10,
  NB_width,
  curve);

  boundary({TL});
  myprop(muv,fs,lambda[0]); // MANDATORY, the interface has moved !!!!!
  mgT = diffusion (TL,dt,muv);
  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[1]);
  boundary({TS});
  mgT = diffusion (TS,dt,muv);
  writefile(mgT);
  invertcs(cs,fs);
  myprop(muv,fs,lambda[0]);
}