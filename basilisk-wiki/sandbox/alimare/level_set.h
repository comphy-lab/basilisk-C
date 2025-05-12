/**
# Level-Set events

The level-set method relies on the use of $\phi$ a higher-dimensional function,
whose zero level-set $\Gamma$ (which is a hypersurface) is an interface between
two fluides.

In our cases, $\phi$ is initialized as a signed distance function.

The level set function is defined and initialized elsewhere (typically by the 
user), the face vector field `uf` and the timestep `dt` are defined by a
solver.

The interface is advected using:
$$
\partial_t\phi+\mathbf{u_f}\cdot\nabla \phi=0
$$
where $\mathbf{u_f}$ is the velocity field and $\phi$ is the level set 
function.
 */

/**
In this pointer we store the different components of $\mathbf{u_f}$
*/
scalar TL[], TS[], dist[];
vector vpc[],vpcf[];

scalar * tracers    = {TL};
scalar * tracers2 = {TS};

scalar * level_set  = {dist};
face vector muv[];
mgstats mgT;
double DT_LS;
int itrecons= 30;
double tolrecons = 1.e-10;
int     nb_cell_NB  ;  // number of cells for the NB
double  NB_width ;     // length of the NB


double  latent_heat = 1.;
double  lambda[2] = {1.,1.}; // thermal capacity of each material
double epsK , epsV;
double eps4;
double T_eq = 0.;

scalar curve[];

/**
These variables will be the different parameters for our level-set method. NB
stands for Narrow Band which is an approximation introduced by [Adalsteinsson and
Sethian, 1995](#Adalsteinsson1995).*/


#define LS_face_value(a,i)      ((a[i] + a[i-1])/2.)


/**
I use a few personal functions in the other level-set related functions.
*/
#include "alex_functions.h"

/**
This set of functions helps you extrapolate the value of the distance beyond a
wall
*/

foreach_dimension()
double distBeyondWall_x(Point point, scalar dist,int sign){
  double delt = dist[sign] - dist[];
  return sign2(delt)*sq(delt)/Delta;
}


#if TREE
event defaults (i = 0) {
  for (scalar s in tracers){
    s.refine = s.prolongation = refine_embed_linear;
    s.restriction = restriction_embed_linear;
  }
  for (scalar s in tracers2){
    s.refine = s.prolongation = refine_embed_linear2;
    s.restriction = restriction_embed_linear2;
  }
}
#endif


#if EMBED
/**
Inverts cs and fs fields.
*/
void invertcs(scalar cs, face vector fs){
  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face(){
    fs.x[]      = 1.-fs.x[];
  }

  boundary({cs});
  restriction({cs});
  fractions_cleanup(cs,fs);
}

/**
## Level-set function to volume and face fractions
*/
struct LS2frac{
  scalar dist;
  scalar cs;
  face vector fs;
  double s_clean;
};

void LS2fractions(struct LS2frac p){
  scalar dist = p.dist;
  scalar cs = p.cs;
  face vector fs = p.fs;
  double s_clean = p.s_clean;

  vertex scalar distn[];
  cell2node(dist,distn);
  fractions (distn, cs, fs);
  fractions_cleanup(cs,fs,smin = s_clean, opposite = true);

  boundary({cs});
  restriction({cs});
}
#endif



/**
~~~bib

@Article{Adalsteinsson1995,
  author    = {Adalsteinsson, David and Sethian, James A},
  title     = {A fast level set method for propagating interfaces},
  journal   = {Journal of computational physics},
  year      = {1995},
  volume    = {118},
  number    = {2},
  pages     = {269--277},
  publisher = {Elsevier},
}
~~~
*/

