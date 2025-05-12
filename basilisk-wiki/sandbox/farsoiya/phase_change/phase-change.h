


/**
# Advection/diffusion of a soluble tracer with phase change

We consider the transport and diffusion of a tracer $c$ couple with phase change. 
Thanks to Quentin and Arnaud for their open source code [elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h).
Numerical framework and implementation details can be found in [Magdelaine, 2019](#magdelaine2019hydrodynamics).


The diffusion coefficient in the phase where tracer diffuses is $D$ and
the jump in concentration at the interface is given by
$$
c_1 = \alpha c_2
$$
The advection/diffusion equation for $c$ can then be written
$$
\partial_t c + \nabla\cdot(\mathbf{u} c) = 
   \nabla\cdot\left(D\nabla c \right)
$$
with $f$ the volume fraction. 


The diffusion coefficient and solubility are attributes of each stracer.

The *stracers* list of soluble tracers must be defined by the calling code. */

attribute {
  scalar phi1, phi2; // private
  double  alpha;  // solubility
  double D; // Diffusivity of tracer in the diffusive phase
  double tr_eq; // Static concentration of tracer in non-diffusive phase
  double  mw;  // Molecular weight of the tracer
}

extern scalar * stracers;


/**
If not defined by the user we fixe a default value for F_ERR, the accepted error
over f to avoid division by zero. */

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

#include "diffusion.h"

/**

## Phase change velocity 

This function evaluates the local gradient of the diffusive tracer in order to
compute the phase change velocity.

The inputs of the function are:

* $f$: VOF tracer,
* $c$: diffusive tracer field,
* $\mathbf{v}_{pc}$: the phase change velocity. */

void phase_change_velocity (scalar f, face vector v_pc) {
 // Multicomponent capability doesn't work now but we will introduce in future
	// Currently use only one stracer
for (scalar c in stracers) {		                                  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, c});
  
  /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
  \mathbf{v}_{pc} = \frac{\mathrm{M}_w}{\rho}\, D\, \nabla c
  $$
  
  we need the tracer gradient $\mathbf{gtr}$ in vapor. It is a priori not well
  defined in a cell crossed by the interface since there is liquid in it.
  Therefore we need to average the values of the gradients in the vapor neighbor
  cells.
  
  Thus, to use the right $\Delta$ in the computation of the gradient, we have to
  compute it before the main loop of the function: */
  
  face vector gtr[];
  
  foreach_face(){
	gtr.x[] = (c[] - c[-1])/Delta;
  }
  
  boundary((scalar*){gtr});

  /**
  To find the vapor neighbor cells and weight the averaging between them, we
  compute the normal (normalized w.r.t. the norm-2) to the interface. We
  have to compute it before the main loop, and not locally, to apply *boundary()*
  to it and get consistent values in the ghost cells. */

  vector n[];
  
  foreach(){
	  coord no = mycs (point, f);
  double nn = 0.;
  foreach_dimension()
    nn += sq(no.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x[] = no.x/nn;
  }
boundary((scalar*){n});
  /**
  With the concentration gradient and the normal vector we can now compute the
  phase change velocity $\mathbf{v}_{pc}$, following the lines drawn in
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
      coord nf, nfm1;
      foreach_dimension(){
        nf.x = 0.;  
        nfm1.x = 0.;  
      }
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
	       foreach_dimension(){
		       nf.x += n.x[-1];
		  nfm1.x = n.x[-1]; 
		}
	}

      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      //FIXME
      if (norm < 1.e-10){
	foreach_dimension()
	    nf.x -= nfm1.x;
	 foreach_dimension()
	      norm += fabs(nf.x);
      }
      foreach_dimension()
        nf.x /= (c.inverse ? norm : - norm);
      
      /**
      We compute the phase change velocity. */

	if (nf.x > 0.) {
	  v_pc.x[] = fabs(nf.x)*gtr.x[1, 0] 
	    + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1] )
#if dimension > 2
	    + fabs(nf.z)*(nf.z > 0. ? gtr.x[1, 0,1] : gtr.x[1, 0, -1]) 
#endif
	    ;
	}
	else if (nf.x < 0.) {
	  v_pc.x[] = fabs(nf.x)*gtr.x[-1, 0]
	    + fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1])
#if dimension > 2
	    + fabs(nf.z)*(nf.z > 0. ? gtr.x[-1, 0, 1] : gtr.x[-1, 0, -1])
#endif
	    ;
	}
    

      v_pc.x[] *= fm.x[]*c.mw/(c.inverse ? rho1 : rho2)*c.D;

      }
    }
  } // for all stracers
  boundary((scalar *){v_pc});
}
/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

#if TREE
event defaults (i = 0)
{
  for (scalar s in stracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
  }
}
#endif // TREE

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 2.3 of [Farsoiya et al.,
2020](#farsoiya2020). */

static scalar * phi_tracers = NULL;

event vof (i++)
{
  phi_tracers = f.tracers;
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;
    
    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c \frac{\alpha f}{\alpha f + (1 - f)}
    $$
    $$
    \phi_2 = c \frac{1 - f}{\alpha f + (1 - f)}
    $$
    */
		  
    foreach() {
      double a = c[]/(f[]*c.alpha + (1. - f[]));
      phi1[] = a*f[]*c.alpha;
      phi2[] = a*(1. - f[]);
    }
    boundary ({phi1, phi2});
  }
}

/**
## Diffusion with an immersed dirichlet condition 

The advected concentration is computed from $\phi_1$ and $\phi_2$ as
$$
c = \phi_1 + \phi_2
$$
and these fields are then discarded. 

The diffusion equation for $c$ is then solved using the implicit
discretisation
$$
\frac{c^{n + 1} - c^n}{\Delta t} = 
\nabla\cdot (D \nabla c^{n + 1}) - \frac{c - \Xi}{\tau_c}\mathcal{H} $$

where $\mathcal{H}$ is heavyside function for the non diffusive phase. $\Xi$ is $\alpha c_{sat}$ and $c_{sat}$ in the interfacial and non-interfacial cells respectively.
$$\tau_{\mathrm{c}}=\frac{1}{10^6} \tau_{\mathrm{D}, \Delta} \quad \text { with } \quad \tau_{\mathrm{D}, \Delta} \equiv \frac{\Delta^{2}}{D_{\mathrm{a}, \mathrm{s}}}$$
see [Magdelaine, 2019](#magdelaine2019hydrodynamics). */
    
/**

The inputs of the function are:

* $c$: diffusive tracer field,
* $f$: VOF tracer,
* *max_level*: maximal level in the simulation,
* $dt$: timestep,
 */

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;
  for (scalar c in stracers) {
    scalar phi1 = c.phi1, phi2 = c.phi2;
    foreach()
      c[] = phi1[] + phi2[];
    delete ({phi1, phi2});

#if TREE
  int max_level = MAX_LEVEL;
#else
  int max_level = LEVEL;
#endif
    
  double dirichlet_time_factor = 1.e6;

    
    scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
    face vector diffusion_coef[];

    double time_factor = dirichlet_time_factor;
  foreach() {
    volumic_metric[] = cm[];
	
    if (interfacial(point, f)){
      dirichlet_source_term[] = cm[]*c.tr_eq*c.alpha*c.D*time_factor*(c.inverse ? f[] : 1. - f[])*sq((1<<max_level)/L0);
    }
    else{
      dirichlet_source_term[] = cm[]*c.tr_eq*c.D*time_factor*(c.inverse ? f[] : 1. - f[])*sq((1<<max_level)/L0);
     }
    dirichlet_feedback_term[] = - cm[]*c.D*time_factor*(c.inverse ? f[] : 1. - f[])*sq((1<<max_level)/L0);

  }
 foreach_face()
    diffusion_coef.x[] = fm.x[]*c.D;

    boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});
    diffusion (c, dt, D = diffusion_coef, r = dirichlet_source_term,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
    boundary({c});

  }
}


/**
## Phase change velocity

The velocity due to evaporation is computed in the *stability()* event to take
into account this velocity in the CFL condition. Before to modify $\mathbf{u}_f$
we save it in another face vector field, in order to recover it just after the
advection of the interface. */

face vector uf_save[];

event stability (i++) {

	face vector ev[];

  phase_change_velocity (f, ev);
 
  foreach_face() {
    uf_save.x[] = uf.x[];
    uf.x[] += ev.x[];
  }
  boundary((scalar*){uf});
  
  
}

/**
After the *vof()* event, the phase change velocity has to be set back to the
*real* velocity. The phase change velocity is not a flow velocity but just a 
displacement of the interface. */

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar*){uf});
}



/**
## References

~~~bib

@phdthesis {magdelaine2019hydrodynamics,
  title = {Hydrodynamics of liquid films h {\ 'e} t {\' e} rog {\ `e} nes},
  author = {Magdelaine-Guillot de Suduiraut, Quentin},
  year = {2019},
  school = {Sorbonne university {\ 'e}}
}

@article{farsoiya2021,
  title = {Bubble mediated single component gas transfer in homogeneous isotropic turbulence},
  author = {P. K. Farsoiya and Q. Magdelaine and A. Antkowiak and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2021},
  note = {submitted}
}
~~~
*/
