/** 
# Mende-Thompson's Models 

This file passes a variable structure viscosity ($\mu_s(\Lambda)$) and a variable relaxation time ($\lambda(\Lambda)$) to the viscoelastic log-conform.h solver. The Mendes and Thompson's model (Souza Mendes and Thompson, 2013) is base on the viscoelastic Oldroyd-B model, where the fully unstructured material viscosity ($\mu_\infty$) is a constant, and the structural elastic modulus ($G(\Lambda)$) is a function of the structure parameter ($\Lambda$). Plasticity is introduced through $\mu_s(\Lambda)$. The structure parameter expresses the state of the structure. By definition, it ranges from 0 to 1, 0 corresponding to a completely unstructured state and 1 corresponding to a completely structured state.

From the  Oldroyd-B model we have:

$\mathbf{\tau} = \mathbf{\tau_s} + \mathbf{\tau_n}$.

Here, $\mathbf{\tau_n} = 2 \mu_\infty \mathbf{D}$, where $\mathbf{D}$ is the deformation tensor. For $\mathbf{\tau_s}$ we have,

$\frac{1}{\mu_s(\Lambda)}\mathbf{\tau_p} +  \frac{1}{G(\Lambda)}\stackrel{\bigtriangledown}{\mathbf{\tau_p}}  = \dot{\mathbf{\gamma}}$

$\mu_s(\Lambda) = \left[ \left( \frac{\mu_0}{\mu_\infty} \right)^\Lambda -1 \right] \mu_\infty$

$G(\Lambda) = G_0\exp(m \left( \frac{1}{\Lambda} - 1 \right))$

$\lambda(\Lambda) = \frac{\mu_s(\Lambda)}{G(\Lambda)}$

where $\mu_0$ and $G_0$ are the values of viscosity and elastic modulus in the fully structured state, respectively, and $m$ is a positive constant.

A kinetic equation that governs the evolution of the structure parameter is introduced in order to describe the non-equilibrium micro-structured states. This equation usually obeys the following backbone (material derivative):

$\frac{d\Lambda}{dt} = \dot{\Lambda}_f - \dot{\Lambda}_d$

where $\dot{\Lambda}_f$ represents the rate of structure formation, while $\dot{\Lambda}_d$ represents the rate of destruction of the bonds of the material. An equilibrium state is achieved when $\dot{\Lambda}_f = \dot{\Lambda}_d$. In the present file, another similar backbone is used to represent the evolution equation for the structure parameter. We use the difference between the current value of $\Lambda$ and the equilibrium value, $\Lambda_{eq}$, as a driving potential for structure changes.

$\frac{\partial \Lambda}{\partial t} + \nabla \cdot (\mathbf{u}\Lambda) = \frac{1}{t_{eq}} \left(1 - \frac{\Lambda}{\Lambda_{eq}(\tau_d)} \right)$

$t_{eq}$ is the thixotropic equilibrium time. When $t_{eq} \rightarrow 0$, $\Lambda \rightarrow \Lambda_{eq}$ faster, and the thixotropic behavior vanishes. $\Lambda_{eq}$ is a function of the norm of the deviatoric part of stress tensor. At steady state:

$\mu_s(\Lambda_{eq}) = \mu_{s eq} = \left[ \left( \frac{\mu_0}{\mu_\infty} \right)^{\Lambda_{eq}} -1 \right] \mu_\infty$

$\Lambda_{eq}$ is obtained through the inverse of the previous equation:

$\Lambda_{eq} = \frac{\ln (\mu_{seq}/mu_\infty)}{\ln (\mu_0/\mu_\infty)}$

$\mu_{seq} = \frac{||\mathbf{\tau_d}||}{||\mathbf{\dot{\gamma}_{eq}}||}$

where $||\mathbf{\dot{\gamma}_{eq}}||$ is the intensity of the strain rate tensor if the flow were at steady state and with a stress magnitude equal to $||\mathbf{\tau_d}||$. Taking the flow curve as, $||\mathbf{\dot{\gamma}_{eq}}||$ should satisfy:

$||\mathbf{\tau_d}|| =  \left[ 1 - \exp \left( - \frac{\mu_0 ||\mathbf{\dot{\gamma}_{eq}}||}{\tau_y} \right) \right]\left( \tau_y + K ||\mathbf{\dot{\gamma}_{eq}}||^n \right) + \mu_\infty ||\mathbf{\dot{\gamma}_{eq}}||$


## Code
*/
  
#include "log-conform.h"

/** The model can also be used with two-phase.h*/
#ifndef VOF
# define VOF 1
#endif
/** For axisymmetric simulations*/
#ifndef AXI
# define AXIS 1
#endif

double tau_y;    // Yield stress
double K;        // Consistency index
double fi;       // Flow index
//double lamb_c;   // Characteristic relaxation time
double mupp;     // Characteristic structure viscosity
double mui;      // Unstructured material viscosity (only for monophasic flows)
double mu0;      // Structured material viscosity
double rhoc;     // Density (only for monophasic flows)
double t_eqq;    // Thixotropic equilibrium time
double gc;       // Characteristic elastic modulus
double g0;       // Structured material elastic modulus
double mm;       // Positive constant of the structure elastic modulus model

tensor shear_t[], stress_t[];
scalar mus[], lamb[];
scalar yielded[], shear_norm[], shear_eq[], stress_norm_dev[], Lambda[], Lambda_eq[], eta_eq[], G[], stress_fc[], div_ul[], shear_tq[], stress_tq[];


event defaults (i = 0) 
{
  lambda = lamb;
  mup = mus;

  foreach()
  {
    shear_eq[] = 1.;
    stress_fc[] = 1.;
  }
  boundary ((scalar *){shear_eq, stress_fc});

#if !VOF
  const face vector muii[] = {mui,mui};
  mu = muii;
  const scalar rhocc[] = rhoc;
  rho = rhocc;
#endif
}

 /** Calculating $\mu_s, G$ and $\lambda$ in the event properties */
event properties (i++)
{
#if !VOF

  foreach()
  {
    mus[] = (pow(mu0/mui, Lambda[]) - 1)*mui;
    G[] = g0*exp(mm*(1/Lambda[] - 1));
    lamb[] = mus[] / G[];
  }
  boundary ({mus, G, lamb});

  double er = 1.e-5;
  /** Calculation of shear_eq, equilibrium viscosity and equilibrium structure parameter*/
  foreach()
  {
    shear_t.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear_t.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear_t.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear_t.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

#if !AXI
    shear_norm[] = sqrt(0.5)*sqrt(sq(shear_t.x.x[])+sq(shear_t.x.y[])+sq(shear_t.x.y[])+sq(shear_t.y.y[]));
#else
    shear_tq[] = 2*u.y[]/y;
    shear_norm[] = sqrt(0.5)*sqrt(sq(shear_t.x.x[])+sq(shear_t.x.y[])+sq(shear_t.x.y[])+sq(shear_t.y.y[])+sq(shear_eq[]));
#endif    

    stress_t.x.x[] = tau_p.x.x[] + shear_t.x.x[] * mui;
    stress_t.x.y[] = tau_p.x.y[] + shear_t.x.y[] * mui;
    stress_t.y.x[] = tau_p.y.x[] + shear_t.y.x[] * mui;
    stress_t.y.y[] = tau_p.y.y[] + shear_t.y.y[] * mui;

#if !AXI
    stress_norm_dev[] = sqrt(sq(stress_t.y.x[]) + sq((stress_t.x.x[] - stress_t.y.y[])/2.);
#else
    stress_tq[] = tau_qq[] + shear_eq[] * mui;
    stress_norm_dev[] = sqrt(0.5)*sqrt(sq(stress_t.y.x[]) + sq(stress_t.x.y[]) + sq((2*stress_t.x.x[]-stress_t.y.y[]-stress_tq[])/3.) + sq((2*stress_t.y.y[]-stress_t.x.x[]-stress_tq[])/3.) + sq((2*stress_tq[]-stress_t.y.y[]-stress_t.x.x[])/3.));
#endif


    while (stress_fc[] - stress_norm_dev[] > er || stress_fc[] - stress_norm_dev[] < -er)
    { 
      if ((stress_fc[] - stress_norm_dev []) > er)
      {
        shear_eq[] = shear_eq[] * 0.5;
      }
      if ((stress_fc[] - stress_norm_dev []) < -er)
      {
        shear_eq[] = shear_eq[] * 1.5;
      }

      stress_fc[] =  (1. - exp(-mu0*shear_eq[]/tau_y)) * (tau_y + K * pow(shear_eq[], fi)) + mui*shear_eq[];
    }
    
    if (shear_eq[] < 1e-5)
      eta_eq[] = mu0;
    else
      eta_eq[] = (stress_fc[]) / (shear_eq[]);  // Equilibrium viscosity

    if (eta_eq[] > mu0)
      eta_eq[] = mu0;

    Lambda_eq[] = log(eta_eq[]/mui) / log(mu0/mui);
  }
  boundary({shear_t, shear_norm, stress_t, stress_norm_dev, stress_fc, shear_eq, eta_eq, Lambda_eq, shear_tq, stress_tq});
  
  /** Advection of Lambda. Calculating Lambda after its advetion*/
  foreach()
  {
    {  
      /** Calculation of $\nabla \cdot (\mathbf{u}\Lambda)$ */
#if !AXI
    div_ul[] = (u.x[1,0]-u.x[-1,0])*Lambda[]/(2*Delta) + (Lambda[1,0]-Lambda[-1,0])*u.x[]/(2*Delta) + (u.y[0,1]-u.y[0,-1])*Lambda[]/(2*Delta) + (Lambda[0,1]-Lambda[0,-1])*u.y[]/(2*Delta);
#else
     div_ul[] = (u.x[1,0]-u.x[-1,0])*Lambda[]/(2*Delta) + (Lambda[1,0]-Lambda[-1,0])*u.x[]/(2*Delta) +(u.y[0,1]-u.y[0,-1])*Lambda[]/(2*Delta) + (Lambda[0,1]-Lambda[0,-1])*u.y[]/(2*Delta) + Lambda[]*u.y[]/y;
#endif

      /**  Calculation of the advected $\Lambda$*/
      Lambda[] = ((1/t_eqq) * (1 - Lambda[] / Lambda_eq[]) - div_ul[])*dt + Lambda[];

      if (Lambda[] > 1.)
        Lambda[] = 1.;
      if (Lambda[] < 0.)
        Lambda[] = 0.;
    }
  }
  boundary({div_ul, Lambda});

 // distinguishes the yielded and unyielded regions                         
 foreach()
 {
    if(stress_norm_dev[] > tau_y)
    {
      yielded[] = 1.;
    }
    else
    {
      yielded[] = 0.; 
    }
  }
  boundary ({yielded});
#endif




/** For two phases flows */
#if VOF

 double ff = 0.99;

  foreach()
  {
    if (f[] > ff)
    {
      mus[] = f[]*(pow(mu0/mui, Lambda_eq[]) - 1)*mui;
      G[] = (g0*exp(mm*(1/(Lambda_eq[]) - 1)));
      lamb[] = f[]*mus[] / G[];
    }

    else if (f[] < 1. - ff)
    {
      mus[] = 0.;
      lamb[] = 0.;
      G[] = g0;
    }

/** Calculation of $\mu_s$ and $\lambda$ on the interface 
*/
    else
    {
      int c1, c2, c3, c4;
      c1 = c2 = c3 = c4 = 1;
      if (f[-1,0] > 1. - ff && f[-1,0] < ff)
       c1 = 0;
      if (f[1,0] > 1. - ff && f[1,0] < ff)
       c2 = 0;
      if (f[0,-1] > 1. - ff && f[0,-1] < ff)
       c3 = 0;
      if (f[0,1] > 1. - ff && f[0,1] < ff)
       c4 = 0;

  /*   if (c1 == 0 && c2 == 0 && c3 == 0 && c4 ==0)
       mus [] = 0.;
    else*/
      mus[] = max((mus[-1,0]*c1 + mus[1,0]*c2 + mus[0,-1]*c3 + mus[0,1]*c4) * f[] / (c1 + c2 + c3 + c4), 0.);
      lamb[] = max((lamb[-1,0]*c1 + lamb[1,0]*c2 + lamb[0,-1]*c3 + lamb[0,1]*c4) * f[] / (c1 + c2 + c3 + c4), 0.);
      G[] = g0;
    }
  }
  boundary ({mus, G, lamb});


  double er = 1.e-5;
  /** Calculation of shear_eq, equilibrium viscosity and equilibrium structure parameter*/
  foreach()
  { 
    if (f[] > ff)
    {
    shear_t.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear_t.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear_t.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear_t.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

#if !AXI
    shear_norm[] = sqrt(0.5)*sqrt(sq(shear_t.x.x[])+sq(shear_t.x.y[])+sq(shear_t.x.y[])+sq(shear_t.y.y[]));
#else
    shear_tq[] = 2*u.y[]/y;
    shear_norm[] = sqrt(0.5)*sqrt(sq(shear_t.x.x[])+sq(shear_t.x.y[])+sq(shear_t.x.y[])+sq(shear_t.y.y[])+sq(shear_eq[]));
#endif    

    stress_t.x.x[] = tau_p.x.x[] + shear_t.x.x[] * mui;
    stress_t.x.y[] = tau_p.x.y[] + shear_t.x.y[] * mui;
    stress_t.y.x[] = tau_p.y.x[] + shear_t.y.x[] * mui;
    stress_t.y.y[] = tau_p.y.y[] + shear_t.y.y[] * mui;

#if !AXI
    stress_norm_dev[] = sqrt(sq(stress_t.y.x[]) + sq((stress_t.x.x[] - stress_t.y.y[])/2.);
#else
    stress_tq[] = tau_qq[] + shear_eq[] * mui;
    stress_norm_dev[] = sqrt(0.5)*sqrt(sq(stress_t.y.x[]) + sq(stress_t.x.y[]) + sq((2*stress_t.x.x[]-stress_t.y.y[]-stress_tq[])/3.) + sq((2*stress_t.y.y[]-stress_t.x.x[]-stress_tq[])/3.) + sq((2*stress_tq[]-stress_t.y.y[]-stress_t.x.x[])/3.));
#endif

    while (stress_fc[] - stress_norm_dev[] > er || stress_fc[] - stress_norm_dev[] < -er)
    { 
      if ((stress_fc[] - stress_norm_dev []) > er)
      {
        shear_eq[] = shear_eq[] * 0.5;
      }
      if ((stress_fc[] - stress_norm_dev []) < -er)
      {
        shear_eq[] = shear_eq[] * 1.5;
      }

      stress_fc[] =  (1. - exp(-mu0*shear_eq[]/tau_y)) * (tau_y + K * pow(shear_eq[], fi)) + mui*shear_eq[];
    }

    eta_eq[] = (stress_fc[]) / (shear_eq[]);  // Equilibrium viscosity

    if (eta_eq[] > mu0)
      eta_eq[] = mu0;

    Lambda_eq[] = f[] * log((eta_eq[])/mui) / log(mu0/mui);
  }

  else if (f[] < 1. - ff)
  {
    shear_norm[] = 0.;
    stress_norm_dev[] = 0.;
    stress_fc[] = 1.;
    shear_eq[] = 1.;
    eta_eq[] = 1.;
    Lambda_eq[] = 0.;
  }

  else
  {
      int c1, c2, c3, c4;
      c1 = c2 = c3 = c4 = 1;
      if (f[-1,0] > 1. - ff && f[-1,0] < ff)
       c1 = 0;
      if (f[1,0] > 1. - ff && f[1,0] < ff)
       c2 = 0;
      if (f[0,-1] > 1. - ff && f[0,-1] < ff)
       c3 = 0;
      if (f[0,1] > 1. - ff && f[0,1] < ff)
       c4 = 0;

      Lambda_eq[] = max((Lambda_eq[-1,0]*c1 + Lambda_eq[1,0]*c2 + Lambda_eq[0,-1]*c3 + Lambda_eq[0,1]*c4) * f[] / (c1 + c2 + c3 + c4), 0.);
  }

  }
  boundary({shear_t, shear_norm, stress_t, stress_norm_dev, stress_fc, shear_eq, eta_eq, Lambda_eq, shear_tq, stress_tq});


  foreach()
  {
    if (f[] > ff)
    {
#if !AXI
    div_ul[] = (u.x[1,0]-u.x[-1,0])*Lambda[]/(2*Delta) + (Lambda[1,0]-Lambda[-1,0])*u.x[]/(2*Delta) + (u.y[0,1]-u.y[0,-1])*Lambda[]/(2*Delta) + (Lambda[0,1]-Lambda[0,-1])*u.y[]/(2*Delta);
#else
     div_ul[] = (u.x[1,0]-u.x[-1,0])*Lambda[]/(2*Delta) + (Lambda[1,0]-Lambda[-1,0])*u.x[]/(2*Delta) +(u.y[0,1]-u.y[0,-1])*Lambda[]/(2*Delta) + (Lambda[0,1]-Lambda[0,-1])*u.y[]/(2*Delta) + Lambda[]*u.y[]/y;
#endif

    Lambda[] = f[]*(((1/(2*dt)) * (1 - Lambda[] / (Lambda_eq[])) - div_ul[])*dt + Lambda[]);

      if (Lambda[] > 1.)
        Lambda[] = 1.;
      if (Lambda[] < 0.)
        Lambda[] = 0.;
    }

    else if (f[] < 1 - ff)
    {
      div_ul[] = 0.;
      Lambda[] = 0.;
    }

    else
    {
      int c1, c2, c3, c4;
      c1 = c2 = c3 = c4 = 1;
      if (f[-1,0] > 1. - ff && f[-1,0] < ff)
       c1 = 0;
      if (f[1,0] > 1. - ff && f[1,0] < ff)
       c2 = 0;
      if (f[0,-1] > 1. - ff && f[0,-1] < ff)
       c3 = 0;
      if (f[0,1] > 1. - ff && f[0,1] < ff)
       c4 = 0;

      Lambda[] = max((Lambda[-1,0]*c1 + Lambda[1,0]*c2 + Lambda[0,-1]*c3 + Lambda[0,1]*c4) * f[] / (c1 + c2 + c3 + c4), 0.);
    }
  }
  boundary({div_ul, Lambda});


 foreach()
 {
    if(stress_norm_dev[] > tau_y*f[] || f[] < ff)
    {
      yielded[] = 1.;
    }
    else
    {
      yielded[] = 0.; 
    }
  }
  boundary ({yielded});

#endif
}
                             

/** 
## References

* [Cassio M. Oishi, Roney L. Thompson, Fernando P. Martins, Transient motions of elasto-viscoplastic thixotropic materials subjected to an imposed stress field and to stress-based free-surface boundary conditions, International Journal of Engineering Science, Volume 109, 2016, Pages 165-201, ISSN 0020-7225, https://doi.org/10.1016/j.ijengsci.2016.08.004.](https://www.sciencedirect.com/science/article/pii/S0020722516308175)

*/