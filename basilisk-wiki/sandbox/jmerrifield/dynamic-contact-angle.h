/**
# Dynamic contact angle
This header computes two models of dynamic contact angle on embeded boundaries for a sessile : Simplified Cox-Voinov law and a model proposed by Afkhami.
*/

/**
## Simplified Cox-Voinov law
If $\theta_d < 3\pi/4$, Cox-Voinov law simplifies into :
$$
\theta_d^3 = \theta_S^3 + 9Ca\ln\left(\frac{L}{\lambda}\right)
$$
with Ca the capillary number, L a macroscopic length and $\lambda$ a microspic length. 
*/
scalar cox-voinov(SessileData c, double theta0, double L, double lambda)
{
  scalar theta_d;
  
  double Ca_l = mu1*c.U_left/f.sigma;
  double Ca_r = mu1*c.U_right/f.sigma;
  
  double theta_l = pow(pow(theta0,3) + Ca_l*log(L/lambda),1./3.);
  double theta_r = pow(pow(theta0,3) + Ca_r*log(L/lambda),1./3.);
  
  foreach()
  {
    // Left side
    if(fabs(x - c.x_left) < 4.*Delta)
    {
      theta_d[] = theta_l;
    }
    
    // Right side
    if(fabs(x - c.x_right) < 4.*Delta)
    {
      theta_d[] = theta_r;
    }
  }
  return theta_d;
}






/**
## Afkhami's model
Afkhami proposed a model based for 2D simulations and Cox expression that reduced to :
$$
\cos\theta_d = \cos\theta_S + 5.63 Ca \ln\left(\frac{K}{\Delta/2}\right)
$$
with Ca the capillary number, K a macroscopic length based on the contact zone and $\Delta$ the spatial step.
This model is valid if $\left|cos\theta_d\right|<0.6$.
*/
scalar afkhami(SessileData c, double theta0, double K)
{
  scalar theta_d;
  
  double Ca_l = mu1*c.U_left/f.sigma;
  double Ca_r = mu1*c.U_right/f.sigma;
  
  double theta_l = acos(clamp(cos(theta0) + 5.63*Ca_l*log(K/(L0/(2.*N))),-1.,1.));
  double theta_r = acos(clamp(cos(theta0) + 5.63*Ca_r*log(K/(L0/(2.*N))),-1.,1.));
  
  foreach()
  {
    // Left side
    if(fabs(x - c.x_left) < 4.*Delta)
    {
      theta_d[] = theta_l;
    }
    
    // Right side
    if(fabs(x - c.x_right) < 4.*Delta)
    {
      theta_d[] = theta_r;
    }
  }
  return theta_d;
}