/**
# Transport of a species concentration field with phase change

This module can be used to compute the concentration field of a soluble species in a two-phase system with mass transfer; the volume change due to the exchange of mass between the phases is taken into account.

Here we implement the two-scalar method of [Fleckenstein and Bothe,2015](#fleckenstein2015), where two transport equations are solved for each species. However, we consider pure disperse (gas) phases only (i.e. no mixtures), where the concentration field is uniform (and known). Therefore, we need to solve only one transport equation in the continuous (liquid) phase:
$$
\partial_t c_c + \mathbf{u}\cdot\nabla c_c = \nabla\cdot\left(D\nabla c_c\right) - \frac{\dot{m}}{M}\delta_\Sigma 
$$
where:

* $c_c$ is the phase-average concentration field in the continuous (liquid) phase
* $D$ is the diffusion coefficient in the continuous phase
* $\dot{m}$ is the mass transfer rate 
* $M$ is the molar mass
* $\delta_\Sigma$ is the interfacial delta function

We assume the primary phase to be the continuous one (i.e. $f=1$ in the liquid domain). The cell-average concentration ($tr_c$) is related to the phase average one via: $tr_c=fc_c$.

The interface is considered to be at the equilibrium (i.e. saturated) and the concentration jump across the interface is given by Henry's Law:
$$
(c_c)_\Sigma=\frac{(c_d)_\Sigma}{He}
$$
where $c_d$ is the (known) concentration in the disperse region and $He$ is Henry's constant.

For the details of the following implementation, the reader is referred to [Gennari et al.,2022](#gennari2022)
*/

#include "mass_refine_prolongation-gennari.h"
#include "../fracface-gennari.h"
#include "gradients-gennari.h"
  
/**
We define the attributes for the concentration field.
*/

attribute {
  double D; //diffusion coefficient in the liquid phase
  double mol_mass; //molar mass
  double cd; //uniform concentration in the gas phase
  double He; //Henry's coeff
  bool diff; //we can switch off (on) the diffusion term if diff is false (true)
}

/**
F_ERR is the accepted error over f to avoid division by zero; it is also used to identify interfacial cells.
*/

#ifndef F_ERR
#define F_ERR 1e-10 //tolerance on the interfacial cells
#endif

/**
We need a scalar *m* to store the redistributed mass transfer field. The list *species* contains the (only) transferable tracer and *m_list* is used to store the interfacial (not redistributed) mass transfer term $\dot{m}$. We have the option to switch on the volume change at time *tt*.
*/

scalar m[]; //overall (redistributed) mass transfer term
extern scalar * species; //list containing the only transferable species
extern scalar * m_list; //mass transfer in the interfacial cells

double tt = 0.; //used to start the volume change at t = tt

/**
## Defaults
*/

//We set the refine/prolongation/restriction operators for Tree grids
event defaults (i = 0)
{
#if TREE
  for (scalar m_i in m_list) {
    m_i.c = f;
    m_i.refine = m_i.prolongation = mass_refine;
    m_i.restriction = no_restriction;
    m_i.dirty = true;
  }

  m.refine = m.prolongation = refine_injection;
  m.dirty = true;
#endif
}

/**
## Tracer advection
The advection of the soluble tracer is performed in the vof event (along with the interface). Here we take into account the contribution of the mass transfer $\dot{m}$ on the species concentration (last term on the RHS of the transport equation) and volume of fluid. We also compute the redistributed mass transfer $m$.
*/

event tracer_advection (i++)
{
  scalar m_i;
  
  /**
  MASS_TRANSFER_RATE is used to set a constant $\dot{m}$ (only useful for debugging purposes).
  */
#ifndef MASS_TRANSFER_RATE
  scalar tr_c;
  for (tr_c, m_i in species, m_list) {
    bool third = false; //we use the scheme of Fleck. and Bothe (VOF avg.)
    
    /**
    We recover the phase-average concentration field $c_c$.
    */
  
    //tr_c -> c_c
    foreach()
      tr_c[] = f[] > F_ERR ? tr_c[]/f[] : 0.;
    boundary({tr_c});
    
    /**
    To compute the gradient of concentration normal to the interface we need the face fraction of the vof (liquid) field.
    */

    face vector d[]; //d is the face fraction of f
    face_fraction (f, d);
    
    /**
    Here we compute the mass transfer rate:
    $$
    \dot{m} = -\frac{MD}{1-\rho_{c_c}/\rho_c}\frac{\partial c_c}{\partial \mathbf{n}_\Sigma}
    $$
    where $\rho_{c_c}$ is the local density of the soluble species on the liquid side of the interface and $\rho_c$ is the liquid density. Note the for large density ratios, the term $1-\rho_{c_c}/\rho_c \approx 0$.
    */
    
    //compute m_true and set m to zero
    foreach() {
      m_i[] = 0.; //set individual mass transfer term = 0
      m[] = 0.; //set global (redistributed) mass transfer term = 0
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
	coord n = interface_normal(point, f);
#if dimension == 2
	coord p = {0.,0.};
#else //dimension ==3
	coord p = {0.,0.,0.};
#endif
	double alpha = plane_alpha (f[], n);
	plane_area_center (n, alpha, &p);
	normalize (&n); // |n| = 1
	double bc = tr_c.cd/tr_c.He; //Henry's Law
	m_i[] = -concentration_gradient (point, tr_c, f, d, n, p, bc, third)*tr_c.mol_mass*tr_c.D; //Fix me. Add 1-\rho_{c_c}/\rho_c term
      }
    }
    
    /**
    We recover the cell-average concentration $tr_c$.
    */

    //c_c -> tr_c
    foreach() {
      tr_c[] *= f[];
    }
  }
  
  /**
  If we use a constant mass transfer rate, we set $\dot{m} = MASS\_TRANSFER\_RATE$.
  */
  
#else //MASS_TRANSFER_RATE
  
  for (m_i in m_list) {
    foreach() {
      m_i[] = 0.; //set individual mass transfer term = 0
      m[] = 0.; //set global (redistributed) mass transfer term = 0
      if (f[] > F_ERR && f[] < 1. - F_ERR)
	m_i[] = MASS_TRANSFER_RATE;
    }
  }
  
#endif
  
  /**
  The incompressible VOF scheme implemented in Basilisk ([vof.h](/src/vof.h)) ensures mass conservation provided the velocity field is divergence-free everywhere (always true for incompressible flows without mass transfer). When phase-change occurs, the velocity field at the interface is generally not divergence-free ($\nabla\cdot\mathbf{u} \neq 0$). Here we propose an algorithm to redistribute the original interfacial mass transfer rate $\dot{m}$ into the neighbouring pure gas cells, in order to have a divergence-free velocity field for the cells that contain the continuous phase (both pure liquid and interfacial cells).
  
  We first store in the field *avg* the number of neighbouring pure gas cells (in a $3\times3$ ($3\times3\times3$) stencil in 2D (3D)) for each interfacial cell. Note that we do this when we want to take into account the volume change, i.e. for $t>tt$.
  */
  
  scalar avg[];
  avg.c = f;

#if TREE
  avg.refine = avg.prolongation = refinement_avg;
  avg.restriction = no_restriction;
  avg.dirty = true;
#endif

  if (t > tt) {
    //Compute avg
    foreach() {
      avg[] = 0.;
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
	int count = 0;
	foreach_neighbor(1) {
	  if (f[] < F_ERR)
	    count ++;
	}
	avg[] = count;
      }
    }
    
    /**
    Here we compute the redistributed mass transfer term $m$. The idea is to replace the [projection](/sandbox/ggennari/phase_change/poisson.h#projection-of-a-velocity-field) step:
    $$
    \nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_f}{\Delta t} - \frac{\dot{m}}{\Delta t}\left(\frac{1}{\rho_d} - \frac{1}{\rho_c}\right)\delta_\Sigma
    $$
    with:
    $$
    \nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_f}{\Delta t} - \frac{m}{\Delta t}
    $$
    The field $m$ is $\neq 0$ only for pure gas cells close to the interface and takes into account the density jump $\left(\frac{1}{\rho_d} - \frac{1}{\rho_c}\right)$ and the Dirac function (approximated as $\delta_\Sigma = A_\Sigma/V$).
    */
      
    //Compute m
    foreach() {
      if (f[] < F_ERR) {
	double val = 0.;
	foreach_neighbor(1) {
	  if (f[] > F_ERR && f[] < 1. - F_ERR && avg[] > 0) {
	    coord n = interface_normal(point, f), p;
	    double alpha = plane_alpha (f[], n);
	    double area = plane_area_center (n, alpha, &p);
	    for (scalar m_i in m_list) {
#if AXI
	      val += cm[]*m_i[]*(area*(y + p.y*Delta)*(1./rho2 - 1./rho1)/(Delta*y)/avg[]);
#else
	      val += m_i[]*(area*(1./rho2 - 1./rho1)/Delta/avg[]);
#endif
	    }
	  }
	}
	m[] = val;
      }
      else {
	m[] = 0.;
      }
    }
  }
  
  /**
  The change in volume due to the mass transfer is computed. We assume a rigid displacement of the interface along the normal direction $\mathbf{n}_\Sigma$:
  $$
  \mathbf{h} = -\frac{\dot{m}}{\rho_c}\frac{\Delta t}{\Delta}\mathbf{n}_\Sigma
  $$
  The change in intercept after the displacement is:
  $$
  \Delta \alpha = -\frac{\dot{m}}{\rho_c}\frac{\Delta t}{\Delta}\sqrt{n_{\Sigma x}^2 + n_{\Sigma y}^2 + n_{\Sigma z}^2}
  $$
  */

//Update f according to mass transfer \dot{m}
  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha (f[], n);
      double val = 0.;
      for (scalar m_i in m_list)
	val += m_i[]; //total mass transfer
      double delta_alpha = -val*dt*sqrt(sq(n.x) + sq(n.y) + sq(n.z))/rho1/Delta; //total interface shift (in terms of alpha). You must be sure that |nx|+|ny|+|nz|=1. 
      //This is not true when contact.h is used. Fix me.
      
      if (t > tt) {
      	double ff = plane_volume (n, alpha + delta_alpha); //total interface shift (in terms of volume fraction)
      	if (ff > F_ERR && ff < 1. - F_ERR)
      	  f[] = ff;
      	f[] = clamp(f[], 0., 1.);
      }
      
      /**
      The concentration in the liquid domain is updated according to $\dot{m}$.
      */

      coord p;
      double area = plane_area_center (n, alpha, &p);
      scalar tr_c, m_i;
      for (tr_c, m_i in species, m_list) {
	
#if AXI
	tr_c[] -= m_i[]*dt*area*(y + p.y*Delta)/(Delta*y)/tr_c.mol_mass;	
#else
	tr_c[] -= m_i[]*dt*area/Delta/tr_c.mol_mass;
	
#endif
      }
    }
  }
  
  //Set tr_c = 0 in pure gas cells (according to F_ERR)
  foreach() {
    for (scalar tr_c in species) {
      tr_c[] *= (f[] > F_ERR);
    }
  }
}

/**
## Tracer diffusion

In this event we perform the integration of the diffusion term. The concentration field $c_c$ must be confined to the liquid side of the interface. It is important to prevent any diffusion across the interface, as this would represent an artificial mass transfer (since $\dot{m}$ has already been taken into account in the event before). Here we use the same approach proposed by Quentin Magdelaine in [mixtures.h](/sandbox/qmagdelaine/phase_change/mixtures.h#diffusion-with-an-immersed-no-flux-condition). Details of the numerical implementation can be found in [Magdelaine, 2019](#magdelaine2019).

The FV scheme for diffusion with an immersed no flux condition reads:
$$
f\, \frac{\partial c_c}{\partial t} = \frac{1}{V}\,
\sum_\text{cell faces}{f_f\, D\, \frac{\partial c_c}{\partial \mathbf{n}}\, A}
$$
where $f_f$ is the face fraction of the VOF field.
*/

event tracer_diffusion (i++) 
{ 
  for (scalar tr_c in species) {
    scalar theta[]; //volume correction (it takes into account the effective volume of species)
    face vector d[];

#if TREE
    theta.refine = theta.prolongation = fraction_refine;
    theta.dirty = true;
#endif

    foreach() {
      theta[] = cm[]*max(f[], F_ERR);
    }
    
    face_fraction (f, d);
    foreach_face() {
      d.x[] *= fm.x[]*tr_c.D;
    }
    boundary((scalar *){d});
    
    //tr_c -> c_c
    foreach() {
      tr_c[] = f[] > F_ERR ? tr_c[]/f[] : 0.;
    }
    
    if (tr_c.diff) {
      diffusion (tr_c, dt, D = d, theta = theta);
    }
    
    //c_c -> tr_c and set tr_c = 0 in pure gas cells
    foreach() {
      tr_c[] *= f[]*(f[] > F_ERR);
    }
  }
}

/**
## References

~~~bib
@article{gennari2022,
title = {A phase-change model for diffusion-driven mass transfer problems in incompressible two-phase flows},
journal = {Chemical Engineering Science},
volume = {259},
pages = {117791},
year = {2022},
issn = {0009-2509},
doi = {https://doi.org/10.1016/j.ces.2022.117791},
url = {https://www.sciencedirect.com/science/article/pii/S000925092200375X},
author = {Gabriele Gennari and Richard Jefferson-Loveday and Stephen J. Pickering}
}

@article{fleckenstein2015,
title = {A Volume-of-Fluid-based numerical method for multi-component mass transfer with local volume changes},
journal = {Journal of Computational Physics},
volume = {301},
pages = {35-58},
year = {2015},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2015.08.011},
url = {https://www.sciencedirect.com/science/article/pii/S0021999115005306},
author = {Stefan Fleckenstein and Dieter Bothe}
}

@phdthesis {magdelaine2019,
  title = {Hydrodynamics of liquid films h {\ 'e} t {\' e} rog {\ `e} nes},
  author = {Magdelaine-Guillot de Suduiraut, Quentin},
  year = {2019},
  school = {Sorbonne university {\ 'e}}
}
~~~
*/