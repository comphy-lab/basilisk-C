/**
# Phase change of mixtures

The file is to be used in addition to
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
to handle evaporation of mixtures.

A first step is to make diffuse the different compounds of the mixture in their
associated phase, whitout crossing the interface.

If not defined by the user we fixe a default value for F_ERR, the accepted error
over f to avoid division by zero, set FACE_FRACTION_REFINE to 1. */

#ifndef FACE_FRACTION_REFINE
  #define FACE_FRACTION_REFINE 0
#endif

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

/**
We need several header files:

* [diffusion.h](/src/diffusion.h): to make diffuse the diffusive tracer,
* [curvature.h](/src/curvature.h): to use *interfacial()* function,
* [my_functions.h](/sandbox/qmagdelaine/my_functions.h): to use some general
functions like compute the normal in every cell,
* [extend_restrict_fields.h](/sandbox/qmagdelaine/phase_change/extend_restrict_fields.h):
to use *compute_extended_mycs()* which computes the mycs normal on the interface
and extends it to the two first rows of neighbors. */

#include "diffusion.h"
#include "curvature.h"
#include "./../my_functions.h"
#include "./../../lopez/src/fracface.h"
#include "extend_restrict_fields.h"

/**
# Diffusion with an immersed no flux condition 

We consider the diffusion of a tracer $tr$ in a phase represented by a VOF
tracer $f$. We write the discretization of the diffusion equation, at an
interface, with a no flux condition, in 2D:

$$
\frac{\partial tr}{\partial t} = \nabla\cdot(D\, \nabla tr) \quad \text{with} 
\quad \iint{tr\, dS} \sim f\, tr\, \Delta^2 
$$
$$
\iint{\frac{\partial tr}{\partial t} \, dS} =
\int{D\, \nabla tr \cdot \mathbf{n}\, d\ell} \quad \text{then} \quad
\Delta^2\, f\, \frac{\partial tr}{\partial t} =
\sum_\text{faces}{\Delta\, f_f\, D\, \nabla tr \cdot \mathbf{n}}
$$
$$
f\, \frac{\partial tr}{\partial t} = \nabla\cdot(f_f\, D\, \nabla tr)
$$

Therefore we see that to embody the no flux condition in the diffusion coeffient,
we need to :

* multiplie by $f$ the left term, to take into account the fact that for the
tracer associated to $f$, the volume of the cell is only $f$ time its real
volume, it is a *volume correction*,
* multiplie by the face fraction $f_f$ (face value of $f$) the diffusion
coefficient, to ensure the no flux condition. 

We define here two functions: one uses directly the fully implicit scheme
written in [diffusion.h](/src/diffusion.h), the other use it in a way to 
implemente a Cranck-Nicholson scheme.

## Fully implicit scheme

The inputs of the function are:

* $tr$: diffusive tracer field,
* $f$: VOF tracer,
* $dt$: the time step. */

mgstats no_flux_diffusion (scalar tr, scalar f, double dt) {

  /**
  We allocate fields for the *volume correction*, the face fraction and the
  weighted diffusion coefficient. */

  scalar volume_correction[];
  face vector f_f[], diffusion_coefficient[];

  #if TREE && !AXI
    volume_correction.prolongation = volume_correction.refine = fraction_refine;
  #endif

  #if TREE && !AXI && FACE_FRACTION_REFINE // seems better without
    foreach_dimension() {
      diffusion_coefficient.x.prolongation = fraction_refine;
      diffusion_coefficient.x.refine = fraction_refine;
    }
  #endif
  
  /**
  If the tracer is associatied to the phase $f=0$ instead of the phase $f=1$,
  the *volume correction* becomes $1 - f$ and the face fraction has to be
  replaced by $1 - f_f$. */ 

  foreach() {
    f[] = clamp(f[], 0., 1.);
    volume_correction[] = cm[]*max(tr.inverse ? 1. - f[] : f[], F_ERR);
  }
  boundary ({f, tr, volume_correction});

  /**
  To compute the face fraction, we use the function of Jose-Maria Lopez Herrera
  defined in [fracface.h](/sandbox/lopez/fracface.h). */

  face_fraction (f, f_f);
  foreach_face()
    diffusion_coefficient.x[] = tr.D*fm.x[]*(tr.inverse ? 1. - f_f.x[] : f_f.x[]);
  boundary((scalar *){diffusion_coefficient});

  /**
  The diffusion equation is solved thanks to [diffusion.h](/src/diffusion.h): */

  return diffusion (tr, dt, D = diffusion_coefficient, theta = volume_correction);
}

/**
## Cranck-Nicholson scheme

i and e indexes refer to *initial* values, before diffusion, and to *end* values,
after diffusion, respectively. The fully implicit scheme is written:
$$
f\, \frac{c^e - c^i}{dt} = \nabla\cdot(f_f\, D\, \nabla c^e) + s
$$
and Cranck-Nicholson one:
$$
2\, f\, \frac{c^e - c^i}{dt} = \nabla\cdot(f_f\, D\, \nabla c^i)
+ \nabla\cdot(f_f\, D\, \nabla c^e) + 2\, s
$$
Therefore, to transform our fully implicit scheme into a Cranck Nicholson one,
we just have to multiply the *volume correction* by 2, and to add to the source
term $s$ the explicit term $\nabla\cdot(f_f\, D\, \nabla c^i)$. */

mgstats no_flux_diffusion_cranck (scalar tr, scalar f, double dt) {

  /**
  We allocate fields for the *volume correction*, the face fraction and the
  weighted diffusion coefficient. */

  scalar volume_correction[];
  face vector f_f[], diffusion_coefficient[];

  #if TREE && !AXI
    volume_correction.prolongation = volume_correction.refine = fraction_refine;
  #endif

  #if TREE && !AXI && FACE_FRACTION_REFINE // seems better without
    foreach_dimension() {
      diffusion_coefficient.x.prolongation = fraction_refine;
      diffusion_coefficient.x.refine = fraction_refine;
    }
  #endif
  
  /**
  If the tracer is associatied to the phase $f=0$ instead of the phase $f=1$,
  the *volume correction* becomes $1 - f$ and the face fraction has to be
  replaced by $1 - f_f$. */ 

  foreach() {
    f[] = clamp(f[], 0., 1.);
    volume_correction[] = 2.*cm[]*max(tr.inverse ? 1. - f[] : f[], F_ERR);
  }
  boundary ({f, tr, volume_correction});

  /**
  To compute the face fraction, we use the function of Jose-Maria Lopez Herrera
  defined in [fracface.h](/sandbox/lopez/fracface.h). */

  face_fraction (f, f_f);
  foreach_face()
    diffusion_coefficient.x[] = tr.D*fm.x[]*(tr.inverse ? 1. - f_f.x[] : f_f.x[]);
  boundary((scalar *){diffusion_coefficient});

  scalar cranck_term[];
  my_laplacian (tr, cranck_term, diffusion_coefficient);
  boundary({cranck_term});
  
  /**
  The diffusion equation is solved thanks to [diffusion.h](/src/diffusion.h): */

  return diffusion (tr, dt, D = diffusion_coefficient, theta = volume_correction,
                    r = cranck_term);
}

/**
## Improvements to do for adaptative grids with axisymetry

The minimum working examples 
[mixture_diffusion_cranck.c](/sandbox/qmagdelaine/phase_change/2_diffusion_in_a_domain/mixture_diffusion_cranck.c)
and
[mixture_diffusion_translation.c](/sandbox/qmagdelaine/phase_change/2_diffusion_in_a_domain/mixture_diffusion_translation.c)
are in fact really fragile: with a larger max level or another value of the ink
patch radius, large values (positive or negative) of the tracer concentration
appears on the interface.
With a constant grid, the problem disappears, and whitout axysymetry, applying
the refine function of the fractions to the *volume correction* seems to correct
the bug. Neverthess, I do not know how write to right refine function nor
boundary condition for a field equal to $f\, cm$. */

/**
# Advection with a non-solenoidal velocity

## Velocity field normal to the interface

For test cases or mwe, it is convenient to have a constant velocity field
defined just at the interface, just as a phase change velocity. This velocity
can also be locally proportional to a given field $c$. */

struct Normal_Velocity {
  // mandatory
  scalar f;
  face vector ev;
  // optional
  double speed; // default 1e-3.
  scalar c; // default is uniform 1
};

void normal_velocity (struct Normal_Velocity p) {

  /**
  We redefine input variables for convenience. */
 
  scalar f = p.f, c = automatic (p.c);
  face vector ev = p.ev, cf[];
  double speed;
  if (p.speed)
    speed = p.speed;
  else
    speed = 1e-3;
    
  vector n[];
  compute_normal (f, n);

  if (p.c.i) { // if c is given
    boundary({c});
    foreach_face()
      cf.x[] = min(c[], c[-1]);
    boundary((scalar*){cf});  
    foreach_face() {
      ev.x[] = 0.;
	    if (interfacial(point, f)) {
		    if (interfacial(neighborp(-1), f))
          ev.x[] = - fm.x[]*speed*(n.x[] + n.x[-1])/2.*(n.x[] > 0. ? cf.x[-1] : cf.x[1]);
        else
          ev.x[] = - fm.x[]*speed*n.x[]*(n.x[] > 0. ? cf.x[-1] : cf.x[1]);
	    }
	    else if (interfacial(neighborp(-1), f))
		    ev.x[] = - fm.x[]*speed*n.x[-1]*(n.x[-1] > 0. ? cf.x[-1] : cf.x[1]);
	    else if (fabs(f[] - f[-1]) > 1. - 2.*F_ERR)
		    ev.x[] = fm.x[]*speed*(f[-1] > f[] ? - cf.x[-1] : cf.x[1]);
	  }
	}
	else { // if c isn't given
    foreach_face() {
      ev.x[] = 0.;
	    if (interfacial(point, f)) {
		    if (interfacial(neighborp(-1), f))
          ev.x[] = - fm.x[]*speed*(n.x[] + n.x[-1])/2.;
        else
          ev.x[] = - fm.x[]*speed*n.x[];
	    }
	    else if (interfacial(neighborp(-1), f))
		    ev.x[] = - fm.x[]*speed*n.x[-1];
	    else if (fabs(f[] - f[-1]) > 1. - 2.*F_ERR)
		    ev.x[] = (f[] > f[-1] ? 1. : -1.)*fm.x[]*speed;
	  }
	}
	boundary((scalar*){ev});
}

/**
## Distribution of the solute lost in the dry cell

After the advection with respect to the phase change velocity, the solute
remaining in the dry cells has to be redistributed in the neighboor cells.
There is several ways to do this redistribution, I propose here to do it
following the phase change velocity.

The inputs of the function are:

* f: the VOF tracer describing the interface, after the advection corresponding
to the phase change,
* qt: the quantity field ($f\,c$) of the solute,
* velocity: the phase change velocity.
*/

void distribution (scalar f, scalar qt, face vector velocity) {

  /**
  To redistribute the solute remaining is the dry cell, we compute fluxes
  between the newly dry cells and their wet neighbors. These fluxes are weighted
  by the phase change velocity, but we need to normalize it as weight
  factors at the scale of a cell. For each cell, if a neighboor cell
  is wet and if the phase change velocity is outward, we add it postively to
  the norm. */

  scalar norm_distrib[];
  foreach() {
    norm_distrib[] = 0.;
    foreach_dimension() {
      norm_distrib[] += (f[-1] > F_ERR && velocity.x[] != nodata ?
                         max(-velocity.x[], 0.) : 0.)
                      + (f[1] > F_ERR && velocity.x[1] != nodata ?
                         max(velocity.x[1], 0.) : 0.);
    }
  }
  boundary({f, velocity, norm_distrib});
  
  /**
  We can now redistribute the solute. If a cell is wet, whereas one neighbor is
  dry, if the phase change velocity is inward and if the norm of the
  distribution is not null, we add to the wet cell its owed part of the solute
  remaning in the dry cell. */
  
  foreach() {
    foreach_dimension(){
      if (f[] > F_ERR && f[-1] < F_ERR && velocity.x[] > 0.
          && velocity.x[] != nodata && norm_distrib[-1] > 0.)
        qt[] += velocity.x[]*qt[-1]/norm_distrib[-1];
      if (f[] > F_ERR && f[1] < F_ERR && velocity.x[1] < 0.
          && norm_distrib[1] > 0.)
        qt[] -= velocity.x[1]*qt[1]/norm_distrib[1];
    }
  }
  foreach()
    qt[] = (f[] < F_ERR ? 0. : qt[]);
  boundary({qt});
}

/**
## Redistribution of the overloaded cells

If we consider the solute concentration as a fraction (molar, mass or volumetric
fraction), it can't exceed one. We define here a function to redistribute the
solute of the overloaded cells. 

The inputs of the function are:

* f: the VOF tracer describing the interface, after the advection corresponding
to the phase change,
* fp: the VOF tracer describing the interface, before the advection corresponding
to the phase change,
* qt: the quantity field ($f\,c$) of the solute,
* velocity: the phase change velocity.
*/

void distribution_over (scalar f, scalar fp, scalar qt, face vector velocity) {
  boundary({f, fp, qt});

  /**
  *extend_face_vector_w0()* is a function of
  [extend_restrict_fields.h](/sandbox/qmagdelaine/phase_change/extend_restrict_fields.h).
  It extends a velocity defined on a interface to the two first rows of
  neighbors. Since the phase change velocity has been defined with respect to
  the previous position of the interface, it has to be extend with the
  corresponding values of the VOF tracer. */

  face vector ev[];  
  extend_face_vector_w0 (velocity, ev, fp);
  
  scalar norm_distrib[];
  foreach() {
    norm_distrib[] = 0.;
    foreach_dimension()
      norm_distrib[] += max(-ev.x[], 0.) + max(ev.x[1], 0.);
  }
  boundary({norm_distrib});
  
  /**
  We compute the flux between the overloaded cells and their neighbors. A
  *concentration* field over 1 corresponds to a *quantity* field over $f$. */
  
  face vector flux[];
  foreach_face() {
    flux.x[] = - (qt[-1] > f[-1] && ev.x[] > 0. && norm_distrib[-1] > F_ERR ?
                  ev.x[]*(qt[-1] - f[-1])/norm_distrib[-1] : 0.)
               + (qt[] > f[] && ev.x[] < 0. && norm_distrib[] > F_ERR ?
                  ev.x[]*(qt[] - f[])/norm_distrib[] : 0.);
  }
  boundary((scalar *) {flux});

  /**
  We apply the computed flux.*/

  foreach() {
    foreach_dimension()
        qt[] += flux.x[] - flux.x[1];
  }
  boundary({qt});
}

/**
I propose a second function to do this redistribution, much simpler because
it does not used the VOF field of the previous step, but whose logic and results
are less convincing. The inputs of the function are:

* f: the VOF tracer describing the interface, after the advection corresponding
to the phase change,
* qt: the quantity field ($f\,c$) of the solute.
*/

void distribution_over_2 (scalar f, scalar qt) {
  boundary({f, qt});
  vector mycs_vector[];
  
  /**
  *compute_extended_mycs()* is a function of
  [extend_restrict_fields.h](/sandbox/qmagdelaine/phase_change/extend_restrict_fields.h).
  It computes the mycs normal on the interface and extends it to the two first
  rows of neighbors. */
  
  compute_extended_mycs (f, mycs_vector);
  
  /**
  We compute the flux between the overloaded cells and their neighbors. A
  *concentration* field over 1 corresponds to a *quantity* field over $f$. */
  
  face vector flux[];

  foreach_face() {
    flux.x[] = - (qt[-1] > f[-1] && mycs_vector.x[-1] < 0. ?
                  mycs_vector.x[-1]*(qt[-1] - f[-1]) : 0.)
               - (qt[] > f[] && mycs_vector.x[] > 0. && mycs_vector.x[] < 2. ?
                  mycs_vector.x[]*(qt[] - f[]) : 0.);
  }
  boundary((scalar *) {flux});

  /**
  We apply the computed flux.*/

  foreach() {
    foreach_dimension()
        qt[] += flux.x[] - flux.x[1];
  }
  boundary({qt});
}

/**
## Improvements to do

It works approximatly as it is on adaptive meshes, but it is far less accurate.
It would be great also to simplify it.
*/
