/**
# Erodible bed (non-steady-state)

This is the non-steady-state generalisation of [this
solver](erosion.h). We solve
$$
\begin{aligned}
  \partial_t z_b + \partial_x q & = 0\\
  \partial_t h_s + \partial_x q & = l_e \tau_s  \bar{\tau}^a  \bar{\omega} 
  (\bar{\tau} - 1)^b - Bh_s \tau\\
  q & = l_a h_s \tau\\
  \tau & = \partial_z u\\
  \bar{\tau} & = \frac{\tau}{\tau_s}
\end{aligned}
$$
where $l_a$ and $l_e$ are characteristic transport and erosion lengths
respectively and $B$ is without dimension. The length $h_s(x,t)$ can be
interpreted as the thickness of the moving sediment layer.

The correspondence can be made with the parameters of the
[steady-state solver](erosion.h) by considering the steady-state
solution of both systems i.e.
$$
\begin{aligned}
\partial_x q & = l_e \tau_s  \bar{\tau}^a  \bar{\omega} 
  (\bar{\tau} - 1)^b - \frac{B}{l_a} q \\
l_s\partial_x q &= E \tau_s  \bar{\tau}^a  \bar{\omega} 
  (\bar{\tau} - 1)^b - q
\end{aligned}
$$
which readily gives the relations
$$
\begin{aligned}
l_e & = \frac{E}{l_s} \\
B & = \frac{l_a}{l_s}
\end{aligned}
$$

A non-erodible "bedrock" level $z_\text{br}$ can optionally be
specified. */

(const) scalar z_br = {-1};

/**
## Default parameters

With the default parameters the bed is not erodible (since $l_a = 0$). */

double tau_s = 0., e_a = 0., e_b = 1.;
double l_a = 0., l_e = 0., B = 0.;

/**
The moving sediment layer has zero initial thickness. */

scalar h_s[], dz[];

event defaults (i = 0) {
  reset ({h_s, dz}, 0.);
#if TREE
  h_s.refine = refine_linear;
  h_s.restriction = restriction_volume_average;
#endif
}

/**
## Erosion/deposition event

The skin friction $\tau = \partial_zu|_{z=z_b}$ is computed using a third-order
accurate discretisation. */

#define dudz(u) (2.*((h[0,0,1]*(h[0,0,1]/h[] + 4.) + 4.*h[])*u.x[] \
		     - h[]*u.x[0,0,1])/				   \
		 (h[0,0,1]*(h[0,0,1] + 3.*h[]) + 2.*sq(h[])))

/**
We use the Bell-Collela-Glaz advection solver to transport the
sediment. */

#include "bcg.h"

event erosion (i++)
{
  scalar tau[];
  foreach()
    tau[] = dudz(u);

  /**
  The sediment advection velocity is $l_a\tau$. */
  
  face vector uf[];
  foreach_face()
    uf.x[] = fm.x[]*l_a*(tau[] + tau[-1])/2.;

  /**
  We compute the fluxes $q = h_s l_a \tau$ using the BCG scheme. */
  
  face vector q[];
  tracer_fluxes (h_s, uf, q, dt, zeroc);

  /**
  If a non-erodible level is defined, we limit the maximum flux. */
  
  if (z_br.i >= 0)
    foreach_face() {
      int i = q.x[] > 0 ? -1 : 0;
      double qmax = Delta*cm[]*(zb[i] - z_br[i])/dt;
      if (fabs(q.x[]) > qmax)
	q.x[] = qmax*sign(q.x[]);
    }

  /**
  The bathymetry and moving sediment layer thickness are updated using
  $$
  \begin{aligned}
  \partial_t z_b + \partial_x q & = 0\\
  \partial_t h_s + \partial_x q & = 0
  \end{aligned}
  $$
  */
  
  foreach() {
    dz[] = (q.x[] - q.x[1])/(Delta*cm[]);
    h_s[] += dt*dz[];
    zb[] += dt*dz[];

    /**
    Since the bottom has changed, we need to update the free surface
    $\eta$. */
    
    eta[] = zb[];
    foreach_layer()
      eta[] += h[];

    /**
    The sediment erosion and deposition are added using a
    time-implicit discretisation of the form
    $$
    \frac{h_s^{t + \Delta t} - h_s^t}{\Delta t} = 
    l_e \tau_s  \bar{\tau}^a  \bar{\omega} (\bar{\tau} - 1)^b 
    - B h_s^{t + \Delta t} \tau
    $$
    */
    
    double h_e = 0., tau_bar = fabs(tau[])/tau_s;
    if (tau_bar > 1.)
      h_e = l_e*tau_s*pow(tau_bar, e_a)*pow(tau_bar - 1., e_b);
    h_s[] = (h_s[] + dt*h_e)/(1. + dt*B*tau_s*tau_bar);
  }
}
