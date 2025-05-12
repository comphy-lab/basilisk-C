/**
# Erodible bed

This adds erosion/deposition to the topography of the [multilayer
solver](/src/layered/README). 

The evolution of the bottom topography $z_b$ is assumed to be
described by the Exner equation ([Exner, 1925](#exner1925))
$$
\partial_t z_b = - \partial_x q
$$
with $q$ the sediment flux.

The sediment flux is assumed to be quasi-stationary (i.e. the bed
evolution is much slower than the flow) and is given by
$$
l_s \partial_x q + q = q_e
$$
where $l_s$ is a characteristic "settling distance" or "saturation
length" ([Andreotti et al., 2002](#andreotti)) and $q_e$ is the
erosion flux.

The erosion flux is assumed to depend only on the skin friction (or
bottom shear stress $\tau$) and is given by ([Kouakou & Lagrée,
2006](#koukaou2006))
$$
q_e = E\tau^a\bar{\omega}(\tau - \tau_s)^b
$$
where the threshold function $\bar{\omega}$ is such that if $(\tau
- \tau_s) > 0$ then $\bar{\omega}(\tau - \tau_s) = \tau - \tau_s$
else $\bar{\omega}(\tau - \tau_s) = 0$.

A non-erodible "bedrock" level $z_\text{br}$ can optionally be
specified. */

(const) scalar z_br = {-1};

/**
## Default parameters

With the default parameters the bed is not erodible (since $E = 0$). */

double l_s = 0., E = 0., tau_s = 0., e_a = 0., e_b = 1., TgsTf = 0.02;
double m,xmax=-100,xmax0=-100,zbmax=0,cdune=0,qmax=0;
/**
The deposition and erosion fluxes are zero initially. */

scalar q[], q_e[], frott[];

event defaults (i = 0) {
  reset ({q, q_e}, 0.);
}

/**
## Erosion/deposition event 

We will need to solve a linear system for $q$. */

#define BGHOSTS 2
#include "solve.h"

/**
The skin friction $\partial_zu|_{z=z_b}$ is computed using a third-order
accurate discretisation. */

#define dudz(u) (2.*((h[0,0,1]*(h[0,0,1]/h[] + 4.) + 4.*h[])*u.x[]	\
		     - h[]*u.x[0,0,1])/					\
		 (h[0,0,1]*(h[0,0,1] + 3.*h[]) + 2.*sq(h[])))

event erosion (i++)
{
  
  /**
  We first compute the erosion flux $q_e$. */

  foreach() {
    
    /**
    The erosion flux is
    $$
    q_e = E\tau^a\bar{\omega}(\tau - \tau_s)^b
    $$
    */
    
    double tau = dudz(u);
    frott[] = tau;
      
    if (tau > tau_s)
      q_e[] = E*(tau - tau_s);//E*pow(tau, e_a)*pow(tau - tau_s, e_b);
    else
      q_e[] = 0.;
    
    /**
    If the "bedrock level" $z_\text{br}$ is specified, the maximum
    erosion flux is given by */

    if (z_br.i >= 0) {
      double qmax = l_s/TgsTf*(zb[] - z_br[])/dt;
      if (q_e[] > qmax)
	q_e[] = qmax;
    }
  }
  boundary ({q_e});
 
  /**
  We get $q$ from the solution of the linear system
  $$
  l_s \partial_x q + q = q_e
  $$
  and a fourth-order upwind discretisation of the gradient. */

  solve (q, l_s*(q[-2] - q[1] + 8.*(q[] - q[-1]))/(6.*Delta) + q[], q_e);
  
  /**
  We then update the bottom topography using Exner's equation
  $$
  \partial_t z_b = - \partial_x q = \frac{q - q_e}{l_s}
  $$  
  */

  foreach() {
    zb[] += TgsTf*dt*(q[] - q_e[])/l_s;

    /**
       Since the bottom has changed, we need to update the free surface
       $\eta$. */
    
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];
  }
  boundary ({zb, eta});
    
    
    

    m=0;
    xmax=-100;
    zbmax=-1;
    foreach(){
        m+=Delta*zb[];
        if(zb[]>=zbmax){zbmax=zb[];xmax=x;qmax=q[];}
    }
    
   
   
    
    
 
    
}

/**
# References

~~~bib
@article{exner1925,
  title={Uber die wechselwirkung zwischen wasser und geschiebe in flussen},
  author={Exner, Felix M},
  journal={Akad. Wiss. Wien Math. Naturwiss. Klasse},
  volume={134},
  number={2a},
  pages={165--204},
  year={1925}
}

@article{charru2006,
  title={Selection of the ripple length on a granular bed sheared 
         by a liquid flow},
  author={Charru, Fran{\c{c}}ois},
  journal={Physics of fluids},
  volume={18},
  number={12},
  pages={121508},
  year={2006},
  publisher={American Institute of Physics}
}

@article{kouakou2006,
  title={Evolution of a model dune in a shear flow},
  author={Kouakou, Kouam{\'e} Kan Jacques and Lagr{\'e}e, Pierre-Yves},
  journal={European Journal of Mechanics-B/Fluids},
  volume={25},
  number={3},
  pages={348--359},
  year={2006},
  publisher={Elsevier},
  pdf = {http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/kouakoulagree06.pdf}
}

@article{andreotti2002,
  title={Selection of dune shapes and velocities Part 1: Dynamics of sand, 
        wind and barchans},
  author={Andreotti, Bruno and Claudin, Philippe and Douady, St{\'e}phane},
  journal={The European Physical Journal B-Condensed Matter and Complex Systems},
  volume={28},
  number={3},
  pages={321--339},
  year={2002},
  publisher={Springer},
  pdf={https://arxiv.org/pdf/cond-mat/0201103}
}
~~~
*/
