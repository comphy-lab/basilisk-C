/**
# Phase change of an elementary body

This file proposes functions to handle the phase change of an elementary body
limited by the diffusion of a diffusive tracer (its vapor for an evaporation or
the temperature for solidification). This file has to be included after
[advection.h](/src/advection.h) or [centered.h](/src/navier-stokes/centered.h)
and [tracers.h](/src/tracers.h).

If not defined by the user we fixe a default value for F_ERR, the accepted error
over f to avoid division by zero. */

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

/**
We use the VOF description of the interface, the Bell-Collela-Glaz advection
solver and the fully implicit reaction-diffusion solver with the following
header files:

* [vof.h](/src/vof.h): to advect the VOF tracer,
* [diffusion.h](/src/diffusion.h): to make diffuse the diffusive tracer,
* [curvature.h](/src/curvature.h): to use *interfacial()* function;
* [my_functions.h](/sandbox/qmagdelaine/my_functions.h): to use some general
functions like compute the normal in every cell. */

#include "vof.h"
#include "diffusion.h"
#include "curvature.h"
#include "my_functions.h"

/**
We add three attribute to the scalar field: a diffusion coefficient $D$, a 
equatilibrium value $t_eq$ and a Peclet number (for the evaporation this is the
ratio between the saturation vapour concentration and the liquid denstity, for
the solification is the ratio between the heat capacity times the a temperature
scale and the latent heat). */

attribute {
  double D;
  double tr_eq;
  double peclet;
}

/**
The modelisation of the evaporation of a pure liquid is devided here in two
parts:

* computation of the advection velocity field for the interface according to
Fick's law,
* diffusion of the diffusive tracer with an immersed dirichlet condition.

One function for each of these steps are defined.

## Phase change velocity 

This function evaluates the local gradient of the diffusive tracer in order to
compute the phase change velocity.

The inputs of the function are:

* $f$: VOF tracer,
* $tr$: diffusive tracer field,
* $\mathbf{v}_{pc}$: the phase change velocity. */

void phase_change_velocity (scalar f, scalar tr, face vector v_pc) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
  \mathbf{v}_{pc} = \mathrm{Pe}\, D\, \nabla tr
  $$
  
  we need the tracer gradient $\mathbf{gtr}$ in vapor. It is a priori not well
  defined in a cell crossed by the interface since there is liquid in it.
  Therefore we need to average the values of the gradients in the vapor neighbor
  cells.
  
  Thus, to use the right $\Delta$ in the computation of the gradient, we have to
  compute it before the main loop of the function: */
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});

  /**
  To find the vapor neighbor cells and weight the averaging between them, we
  compute the normal (normalized w.r.t. the norm-2) to the interface. We
  have to compute it before the main loop, and not locally, to apply *boundary()*
  to it and get consistent values in the ghost cells. */

  vector n[];
  compute_normal (f, n);

  /**
  With the concentration gradient and the normal vector we can now compute the
  evaporation velocity $\mathbf{v}_{pc}$, following the lines drawn in
  [meanflow.c](/sandbox/popinet/meanflow.c). We define it as the product between
  the density ratio and the diffusive flow. Note that $\mathbf{v}_{pc}$ is weighted
  by the face metric. */
  
  foreach_face() {
    v_pc.x[] = 0.;
    
    /**
    Foreach face, we first compute the average value of the normal on the face
    and renormalize it w.r.t the norm-1. We will use this normal to look where
    in the phase in which the tracer diffuse and to weight the average of its
    gradient. We inverse the normal if the tracer is associated the $f=1$ phase
    to take the gradient in the right phase. */
    
    if (interfacial(point, f) || interfacial(neighborp(-1), f)) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
        nf.x += n.x[-1];
        nf.y += n.y[-1];
      }
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= (tr.inverse ? norm : - norm);
      
      /**
      We compute the evaporation velocity. */
      
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] 
                    + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]
                    + fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
      v_pc.x[] *= fm.x[]*tr.D*tr.peclet;
    }
  }
  boundary((scalar *){v_pc});
}


/**
## Diffusion with an immersed dirichlet condition 

An important point is how the immersed Dirichlet boundary condition is
handled. To ensure that the diffusive tracer at the interface stays
at the equilibrium value $t_eq$ and is not washed away, we introduce a source
term in the diffusion equation. This term is activated all over the liquid.
This term is as large as the difference between the actual value and
$t_eq$ is. This modified equation therefore reads:
$$
\frac{\partial tr}{\partial t} = \Delta tr + p\, \delta_{op} \quad \text{with}
\quad p = \frac{tr_{eq} - tr}{\tau_e}
$$
where $\delta_{op}$ is 1 in the *other phase* (where the tracer does not diffuse),
$\tau_e$ is a time, which has to be shorter than the smallest diffusive time in
order for the control to be effective. The smallest diffusive time is $\tau_D =
\frac{\delta^2}{D}$ where $\delta$ is the smallest size in all the tree. Thus,
we choose $\tau_e = \frac{\tau_D}{\alpha}$ with $\alpha$ larger than 1.
Thereafter, $\alpha$ will be refered as *time factor*.

In 2D, the production terme $P$ intergrated over the cell is:
$$
P = \iint_S{p\, \delta_{op}} = p\, f\, \Delta^2
= \frac{tr_{eq} - tr}{\tau_e}\, f\, \Delta^2
$$

The inputs of the function are:

* $tr$: diffusive tracer field,
* $f$: VOF tracer,
* *max_level*: maximal level in the simulation,
* $dt$: timestep,
* $tr_{op}$: another tracer in the *other phase*, unused for elementary bodies,
but requied for mixtures (Raoult law). */

struct Dirichlet_Diffusion {
  // mandatory
  scalar tr;
  scalar f;
  int max_level;
  double dt;
  double time_factor;
  // optional
  scalar tr_op; // default uniform 1.
};

mgstats dirichlet_diffusion (struct Dirichlet_Diffusion p) {

  /**
  We redefine the inputs for convenience. */
  
  scalar tr = p.tr, f = p.f, tr_op = automatic (p.tr_op);
  int max_level = p.max_level;
  double time_factor = p.time_factor;

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  if (p.tr_op.i)
    boundary ({tr_op});

  /**
  We compute the volumic metric of the cells, the dirichlet source terms and
  the diffusion coefficient. If the diffusive tracer is associated to the $f=1$
  phase, we set the production term in the $f=0$ phase. */

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[];

  foreach() {
    volumic_metric[] = cm[];
    if (p.tr_op.i)  
      dirichlet_source_term[] = cm[]*tr.tr_eq*tr_op[]*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    else
      dirichlet_source_term[] = cm[]*tr.tr_eq*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    dirichlet_feedback_term[] = - cm[]*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});

  /**
  The diffusion equation is solved thanks to [diffusion.h](/src/diffusion.h): */
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
}

/**
## Further improvements

Reading the work of Leon Malan, we can list at least three possible improvements
of what it is done here. Two of them are based on the idea to take advantage of
our precise knowledge of the concentration value at the interface (saturation)
and its position (thanks to the height function):

* a better evaluation of concentration gradients at the interface,
* choose another strategy to implemente the dirichlet condition at the interface,
imposing a flux.

The third one is to account for the Stefan flow and the recoil pressure by
adding a term in the projection step of the resolution of the Navier-Stokes
equation. */
