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
#include "./../qmagdelaine/my_functions.h"
#include "./../lopez/fracface.h"
#include "./../qmagdelaine/phase_change/extend_restrict_fields.h"

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

mgstats with_flux_diffusion (scalar tr, scalar source, scalar f, double dt) {

  /**
  We allocate fields for the *volume correction*, the face fraction and the
  weighted diffusion coefficient. */

  scalar volume_correction[], source_term[];
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
    source_term[] = source[];
  }
  boundary ({f, tr, volume_correction, source_term});

  /**
  To compute the face fraction, we use the function of Jose-Maria Lopez Herrera
  defined in [fracface.h](/sandbox/lopez/fracface.h). */

  face_fraction (f, f_f);
  foreach_face()
    diffusion_coefficient.x[] = tr.D*fm.x[]*(tr.inverse ? 1. - f_f.x[] : f_f.x[]);
  boundary((scalar *){diffusion_coefficient});

  /**
  The diffusion equation is solved thanks to [diffusion.h](/src/diffusion.h): */

  return diffusion (tr, dt, D = diffusion_coefficient, theta = volume_correction,
                    r = source_term);
}

/**
## Cranck-Nicholson scheme

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

