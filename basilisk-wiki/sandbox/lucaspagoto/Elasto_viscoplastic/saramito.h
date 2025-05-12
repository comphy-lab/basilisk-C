/** 

# Saramito's Models 
 
This file passes a variable polymeric viscosity ($\mu_p$) and a variable relaxation time ($\lambda$) to the viscoelastic log-conform.h solver. The Saramito's model (Saramito, 2007 and 2009) is base on the viscoelastic Oldroyd-B model, where the solvent viscosity ($\mu_s$) and the elastic modulus ($G$) are kpet constant, and plasticity is introduced through $\mu_p$, which is based on the Bingham or Herschel-Bulkley models.

From the  Oldroyd-B model we have:

$\mathbf{\tau} = \mathbf{\tau_p} + \mathbf{\tau_s}$.

Here, $\tau_s = 2 \mu_s \mathbf{D}$, where $\mathbf{D}$ is the deformation tensor. For $\mathbf{\tau_p}$ we have

$\mathbf{\tau_p} + \lambda \stackrel{\bigtriangledown}{\mathbf{\tau_p}}  = \mu_p \dot{\mathbf{\gamma}}$
or
$\frac{1}{\mu_p}\mathbf{\tau_p} +  \frac{1}{G}\stackrel{\bigtriangledown}{\mathbf{\tau_p}}  = \dot{\mathbf{\gamma}}$

But,

$||\mathbf{\tau_{pd}}|| = \tau_y + K ||\mathbf{\dot{\gamma}_v}||^{n}$

$||\mathbf{\dot{\gamma}_v}|| = \left( \frac{||\mathbf{\tau_{pd}}|| - \tau_y}{K} \right)^{1/n}$

$||\mathbf{\tau_{pd}}|| = \mu_p ||\mathbf{\dot{\gamma}_v}||$

$\mu_p = \left( \frac{K||\mathbf{\tau_{pd}}||^n}{||\mathbf{\tau_{pd}}|| - \tau_y} \right)^{1/n}$

The subscript $d$ means the deviatoric part of the polymeric stress tensor. Plugging the last equation in the equation for $\mathbf{\tau_p}$ (Upper-convected Maxwell model)

$\left( \frac{||\mathbf{\tau_{pd}}|| - \tau_y}{K||\mathbf{\tau_{pd}}||^n} \right)^{1/n} \mathbf{\tau_p} +  \frac{1}{G}\stackrel{\bigtriangledown}{\mathbf{\tau_p}}  = \dot{\mathbf{\gamma}}$

So we have that:

$\lambda = \frac{\mu_p}{G} = \frac{\left( \frac{K||\mathbf{\tau_{pd}}||^n}{||\mathbf{\tau_{pd}}|| - \tau_y} \right)^{1/n}}{G}$

We define the characteristic relaxation time:

$\lambda_c = \lambda \frac{\mu_{pc}}{\mu_{p}} = \frac{\mu_{p}}{G} \frac{\mu_{pc}}{\mu_{p}} = \frac{\mu_{pc}}{G}$

# Code
*/
  
#include "log-conform.h"

/** MODEL = 0 corredponds to the Bingham model, while MODEL = 1 corresponds to the Herschel-Bulkley model*/
#ifndef MODEL
# define MODEL 1
#endif
/** The model can also be used with two-phase.h*/
#ifndef VOF
# define VOF 1
#endif
/** For axisymmetric simulations*/
#ifndef AXI
# define AXIS 1
#endif

double tau_y;      // Yield stress
double K;          // Consistency index
double fi;         // Flow index
double epsilon;    // Regularization parameter
double lamb_c;     // Characteristic relaxation time
double mupp;       // Characteristic polymeric viscosity
double mus;        // Solvent viscosity (only for monophasic flows)
double rhoc;       // Density (only for monophasic flows)

scalar tau_norm[], deviatoric[], yielded[], mupc[], lamb[];

event defaults (i = 0) 
{
  lambda = lamb;   // Relaxation time 
  mup = mupc;      // Polymeric viscosity

#if !VOF  // For monophasic flows
  const face vector muss[] = {mus,mus};
  mu = muss;      // Solvent viscosity
  const scalar rhocc[] = rhoc;
  rho = rhocc;    // Density
#endif
}




 /* Calculating $\mu_p$ and $\lambda$ in the event properties */
event properties (i++)
{
#if !VOF
  foreach()
  {
#if !AXIS
    tau_norm[] = sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]));   //Frobernoius norm of the deviatoric polymeric stress tensor. Obs.: I still have to verify the axisymmetric case.
#else
    tau_norm[] = sqrt(0.5)*sqrt( sq(tau_p.x.y[]) + sq(tau_p.y.x[]) + sq((2*tau_p.x.x[]-tau_p.y.y[]-tau_qq[])/3) + sq((2*tau_p.y.y[]-tau_p.x.x[]-tau_qq[])/3)+ sq((2*tau_qq[]-tau_p.y.y[]-tau_p.x.x[])/3) );
#endif


    if (tau_norm[] <= tau_y + epsilon)   // epsilon is the regularization parameter
    {
     deviatoric[] = tau_y + epsilon;
     yielded[] = 0.;
    }
    else
    {
      deviatoric[] = tau_norm[];
      yielded[] = 1.;
    }
  }
  boundary ({tau_norm, deviatoric, yielded});
  
  foreach()
  {
#if defined  MODEL
    mupc[] = pow( (K * pow(deviatoric[], fi)/(deviatoric[] - tau_y)), 1/fi);  // Herschel-Bulkley model
#else
    mupc[] = deviatoric[]/(deviatoric[] - tau_y);   // Bingham model
#endif
    lamb[] = lamb_c * mupc[] / mupp;
  }
  boundary ({mupc, lamb});
  
 foreach()  // distinguishes the yielded and unyielded regions 
 {
    if(tau_norm[] > tau_y)
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



#if VOF    // For two phases flows
  double ff = 0.99;
  foreach()
  {
#if !AXIS
    tau_norm[] = sqrt(sq((tau_p.x.x[]-tau_p.y.y[])/2) + sq(tau_p.x.y[]));   //Frobernoius norm of the deviatoric polymeric stress tensor. Obs.: I still have to verify the axisymmetric case.
#else
    tau_norm[] = sqrt(0.5)*sqrt( sq(tau_p.x.y[]) + sq(tau_p.y.x[]) + sq((2*tau_p.x.x[]-tau_p.y.y[]-tau_qq[])/3) + sq((2*tau_p.y.y[]-tau_p.x.x[]-tau_qq[])/3)+ sq((2*tau_qq[]-tau_p.y.y[]-tau_p.x.x[])/3) );
#endif

    if (tau_norm[] <= tau_y + epsilon && f[] > ff)
    {
#if MODEL
     mupc[] = max(pow( (K * pow(tau_norm[] + epsilon, fi)/(epsilon)), 1/fi), 0.);   // Herschel-Bulkley fluid
#else  
      mupc[] = max((tau_y + epsilon)/(epsilon), 0.);  // Bingham fluid
#endif    
    }

    else if (f[] < 1. - ff)
      mupc[] = 0.;

    else if (tau_norm[] > tau_y + epsilon && f[] > ff)
    {
#if MODEL
      mupc[] = max(pow( (K * pow(tau_norm[], fi)/(tau_norm[] - tau_y)), 1/fi), 0.);    // Herschel-Bulkley model
#else
      mupc[] = max(tau_norm[]/(tau_norm[] - tau_y), 0.);  // Bingham model
#endif   
    }

/**
It averages the viscosity on the interface using neighboring cells. When the interface reachs a new cell, $\mathbf{\tau_{pd}}$ grows from zero in that cell. Since $\mu_p$ depends on $\mathbf{\tau_{pd}}$, and a small value of $\mathbf{\tau_{pd}}$ gives a very high viscosity (because $\mathbf{\tau_{pd}}  < \tau_y$), we calculate $\mu_p$ on the interface as the average of $\mu_p$ in the neighboring cells (not using neighboring cells also on the interface), instead of calculating $\mu_p$ as a function of $\mathbf{\tau_{pd}}$
*/
    else if (f[] > 1. - ff && f[] < ff)  
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

     if (c1 == 0 && c2 == 0 && c3 == 0 && c4 ==0)
    {
//      printf("%d %d %d %d %g\n", c1, c2, c3, c4, f[]);
      mupc [] = 0.;
    }
    else
      mupc[] = max((mupc[-1,0]*c1 + mupc[1,0]*c2 + mupc[0,-1]*c3 + mupc[0,1]*c4) * f[] / (c1 + c2 + c3 + c4), 0.);
    }

    lamb[] = max(lamb_c*mupc[]/mupp, 0.);
  }
  boundary ({tau_norm, mupc, lamb});


 foreach()
 {
    if(tau_norm[] > tau_y*f[] || f[] < ff)
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

* [Pierre Saramito, A new constitutive equation for elastoviscoplastic fluid flows, Journal of Non-Newtonian Fluid Mechanics, Volume 145, Issue 1,
2007, Pages 1-14, ISSN 0377-0257, https://doi.org/10.1016/j.jnnfm.2007.04.004.](https://www.sciencedirect.com/science/article/pii/S0377025707000869)

* [Pierre Saramito, A new elastoviscoplastic model based on the Herschel–Bulkley viscoplastic model, Journal of Non-Newtonian Fluid Mechanics,
Volume 158, Issues 1–3, 2009, Pages 154-161, ISSN 0377-0257, https://doi.org/10.1016/j.jnnfm.2008.12.001.](https://www.sciencedirect.com/science/article/pii/S0377025708002267)

*/