#include "./solute.h"
#include "./surface.h"
#include "./tension.h"

attribute {
//  scalar c;
  double omega; // default 0 -- To modify for activation
  double alpha; // default 0 -- To modify for activation
  double K;     // default 0 -- To modify for activation
  double max;   // default HUGE -- Modify to activate Langmuir
  double A;     // default 0 -- Modify to activate Frumkin
  double M0;    // default 0 -- For Henry
  double c0;    // default 0 -- For Henry
}

/**double D = 0;
event init (i = 0)
{
  c.D_hor = c.D_vert = D;
}*/

/**
## Surface concentration to surface tension
*/
double sigma_solv = 0.;
event pressure (i++)
{
  foreach() 
    sigma[] = sigma_solv;
    
  for (scalar M in surface)
    if (M.alpha) {
      foreach() {
        if (M.M0) {
          sigma[] -= M.alpha * (M[]/area_m(point) - M.M0);
        }
        else{
          double theta = M[]/(area_m(point) * M.max);
          sigma[] += M.alpha * M.max * (log(1 - theta) - M.A/2 * pow(theta,2));
        } 
      }
    }
  boundary({sigma});
}

/**
## Solute adsorption
Surfactant adsorption follows an adsorption/desorption kinetic, 
which can be modelled using Langmuir isotherm :
$$
\frac{d \Gamma}{dt} = \omega\, (\sigma\, c\, (1-\Gamma) - \Gamma)
\quad \textrm{then} \quad
\frac{dM}{dt} = \omega\, (\sigma\, c\, (\mathcal{A}-M) - M)
$$

The parameters $\omega$ and $\sigma$ are introduced. The adsorption rate $\omega$ 
must be initialized with a non-zero value to add adsorption. */


/**
This function calculate the adsorption term.
*/
scalar adsorption_term (scalar M)
{
//  scalar c = M.c;
  scalar ads[];
  foreach() {
    if (M.M0) {
      ads[] = M.omega * ((M.M0 + M.K * (c[0,0,nl-1] - M.c0)) * area_m(point) - M[]);
    }
    else {
      double theta = M[]/(area_m(point) * (M.max ? M.max : HUGE));
      ads[] = M.omega * (M.K * c[0,0,nl-1] * area_m(point) * (1 - theta) - M[] * exp(- M.A * theta));
    }
  }
  boundary({ads});  
  return ads;
}

/**
To ensure stability, the timestep must be maximized. The event set_dtmax is called to ensure that:
$$
dt<\frac{1}{5 * \omega}
\quad\textrm{and}\quad
dt<\frac{c * h}{10 \, \partial_t \Gamma}
$$ 
*/
event stability (i++) 
{
  for (scalar M in surface)
    if (M.omega && M.K) {
//      scalar c = M.c;
      double dt_max = 1/(5.*M.omega);
      scalar d_ads_t = adsorption_term(M);
      foreach()
        if (d_ads_t[])
          dt_max = min(dt_max, c[0,0,nl-1]*h[0,0,nl-1]/fabs(d_ads_t[])/10.);
      dtmax = min(dtmax, dt_max);
    }
}

/**
The function is called before the event remap (so as to have the new surface concentration for viscosity_event).
Improvements : implicit this step to remove the stability condition */
event remap (i++)
{
  for (scalar M in surface)
    if (M.omega && M.K) {
      scalar d_ads_t = adsorption_term(M);
//      scalar c = M.c;
      foreach(){
        M[] += d_ads_t[]*dt;
        c[0,0,nl-1] -= d_ads_t[]*dt/h[0,0,nl-1];
      }
      boundary({M, c});
  }
}
