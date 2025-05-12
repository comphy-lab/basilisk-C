/**
# A semi-implicit Saint-Venant solver

This solver is based on the paper of [Kwatra et al,
2009](references.bib#kwatra2009) and the implicit solver of Saint-Venant. 

A 1D model for arterial flows can be derived from the Navier-Stokes
equations, in terms of the cross sectional area $a$ and flow rate $q$,
we have
$$
\partial_t a +\partial_x q  = 0 
$$
$$
\partial_t q +\partial_x (q^2/a) = - a \partial_x p (a) /\rho - f_r 
$$ 
where $p(a)$ models the wall properties of the arteries, $\rho$ is the
blood density and $f_r$ stands for the wall shear stress. 

We start with the advection solver and add the Poisson solver. */

#include "advection.h"
#include "poisson.h"

/**
The primitive variables are the cross section areaA $a$ and the flux
$\mathbf{q}=a\mathbf{u}$. */

scalar q[];
scalar a[];

/**
Source terms */

scalar sa[];
scalar sq[];

/**
User defined functions */

double rho = 1. ;
double c2;

double pressurevalue (double h);
double derivativevalue(double h);
void update_source (scalar * conserved, scalar * updates);

/**
Both $\mathbf{q}$ and $a$ are advected using the conservative
advection scheme. */

scalar * tracers = {a,q};

/**
The default slope-limiting is the same as for the explicit
[Saint-Venant solver](saint-venant.h). */

event defaults (i = 0) {
  theta = 1.3;
  gradient = minmod2;
}

/**
The advection solver defines a face velocity field. We initialise it
by interpolation from the flux and the cross sectional area a. */

event init (i = 0) {
  boundary ({q,a});
  foreach_face()
    uf.x[] = (q[] + q[-1])/(a[] + a[-1]);
  boundary ((scalar *){uf});
}

/**
The equation for the pressure is a Poisson--Helmoltz problem which we
will solve with the [multigrid solver](poisson.h). The statistics for
the solver will be stored in *mgp*. */

mgstats mgp;

/**
Up to this point, the advection solver includes events which will
advect $h$ and $\mathbf{q}$ to the next timestep i.e. we now have
$$
\begin{eqnarray}
  a_{n + 1} & = & a_n - \Delta t \nabla \cdot (\mathbf{u}_n a_{n + 1 / 2})\\
  \mathbf{q}_{\star} & = & \mathbf{q}_n - \Delta t \nabla \cdot
  (\mathbf{u}_n \mathbf{q}_{n + 1 / 2})
\end{eqnarray}
$$
We now need to add the pressure gradient term to get
$\mathbf{q}_{n+1}$ from $\mathbf{q}_{\star}$.
*/

event pressure (i++, last)
{

  /**
  We first define a temporary face velocity field $\mathbf{u}_\star$
  using simple averaging from $\mathbf{q}_{\star}$ and $a_{n + 1}$. */

  foreach_face()
    uf.x[] = (q[] + q[-1])/(a[] + a[-1]);
  boundary ((scalar *){uf});

 
 /**
  The evolution equation for the pressure is $$p_t +\mathbf{u} \cdot
  \partial_x p = - \rho c^2 \partial_x \cdot \mathbf{u}$$ with $\rho$
  the density and $c$ the speed of sound. For the arterial model we
  have $p=p(a)$ then $$ \frac{Dp}{dt} = \partial_a p \frac{ Da}{dt} $$
  and using the mass conservation we found the result follows. We
  recall that $c^2 = \frac{ a }{\rho } \partial_a p$.

  Following the classical [projection
  method](navier-stokes/centered.h#approximate-projection) for
  incompressible flows, we set
  $$
  \mathbf{u}_{n + 1} = \mathbf{u}_{\star} - 
                         \Delta t \frac{\nabla p}{\rho^{n + 1}}
  $$
  The evolution equation for the pressure can then be discretised as
  $$
  \frac{p_{n + 1} - p_n}{\Delta t} +\mathbf{u}_n \cdot \nabla p_n = 
     - \rho c^2_{n + 1} \nabla \cdot \mathbf{u}_{n + 1}
  $$
  which gives, after some manipulations, the Poisson--Helmholtz equation
  $$
  \lambda_n p_{n + 1} + \nabla \cdot \left( \frac{\nabla p}{\rho_{}}
  \right)_{n + 1} = \lambda_n p_{\star} + \frac{1}{\Delta t} \nabla \cdot
  \mathbf{u}_{\star}
  $$
  with
  $$
  p_{\star} = p_n - \Delta t\mathbf{u}_n \cdot \nabla p_n
  $$
  and
  $$
  \lambda = \frac{- 1}{\Delta t^2 \rho c^2}
  $$ */
  
  scalar p[], lambda[], rhs[];
  foreach() {

    /**
    To compute the right-hand-side of the Poisson--Helmholtz
    equations, we first compute the divergence of the velocity
    field. */
    
    rhs[] = 0.;
    foreach_dimension()
      rhs[] += uf.x[1] - uf.x[];
    rhs[] /= dt*Delta;

    /**
    The pressure $p_\star$ is obtained simply using the equation of
    state ($p=p(a)$) and the value of $a$ at time $n+1$. */
    
    p[] = pressurevalue(a[]) ;
    
    c2 = a[] / rho * derivativevalue(a[]) ;

    /**
    The $\lambda$ coefficient is computed and the pressure term is
    added to the r.h.s. */
    
    lambda[] = -1./(sq(dt)*rho*c2);
    rhs[] += lambda[]*p[];
  }
  boundary ({p});

  /**
  We then compute the coefficients for the Laplacian operator. */
  
  face vector alpha[];
  foreach_face()
    alpha.x[] = 2./(a[] + a[-1]);
  boundary_flux ({alpha});
  
  /**
  The Poisson--Helmholtz solver is called with a [definition of the
  tolerance](poisson.h#377) identical to that used for incompressible
  flows. */
  
  mgp = poisson (p, rhs, alpha, lambda, tolerance = TOLERANCE/sq(dt));
     
  /**
  The pressure gradient is applied to $\mathbf{u}_\star$ to obtain the
  face velocity field at time $n + 1$. */
  
  foreach_face()
    uf.x[] -= dt*alpha.x[]*(p[] - p[-1])/Delta;
  boundary ((scalar *){uf});
  
  /**
  We define a face pressure using a density-weighted averaging. */
  
  face vector pf = alpha;
  foreach_face()
    pf.x[] = (p[]*a[] + p[-1]*a[-1])/(a[] + a[-1]);
  boundary_flux ({pf});
  
  /**
  And finally we apply the pressure gradient term to the flux/momentum. */
  
  foreach()
    q[] -= dt*a[]*(pf.x[1] - pf.x[])/Delta;
  boundary ((scalar *){q});
}

/**
The last event updates the source terms */

event event_source (i++,last)
{
  update_source ({a,q},{sa,sq});
  foreach() {
    a[] += dt *  sa[] ;
    q[] += dt *  sq[] ;
  }
  boundary ({a,q});
}
