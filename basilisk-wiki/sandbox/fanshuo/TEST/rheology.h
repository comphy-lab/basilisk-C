/** # Rheological and geometry parameters of flows  
## Overview:

Example of local "non Newtonian" rheologies :

*** bla bla bla 

turbulent 
$$\ell=y$$
$$\ell=\sqrt{1-y/h}$$
$$\tau = \kappa \ell^2 \frac{\partial u}{\partial y}^2$$
granular
 $$\tau = \tau_c+\mu(I)P$$
 Bingham 
 $$\tau = \tau_y + \eta   \frac{\partial u}{\partial y}$$
 Herschell-Bulkley 
$$\nu_{eq}= \kappa y^2 \frac{\partial u}{\partial y}$$

*** bla bla bla 

need of $\nu_{eq}$
$$\nu_{eq}=\frac{\tau}{|\frac{\partial u}{\partial y}|}$$
*/
/** 
## To be done
It would be good to add dimensions to certain variables to avoid confusions and conflicts.

For exemples, the "nu" for Poiseuille is the cinematic viscosity, the "mu" for Bingham flow is the dynamic viscosity (also called "eta" in some textbooks). While the "mu0" used by granular flows is Coulomb friction coefficient, which is dimensionless.

It is ok to mix dynamic and cinematic viscosity at current stade because the fluide density is set to one by default. One should be careful if trying to vary densities in the layers.

A global slope is imposed, one has to adapt with `zb` 
*/

/**## Basic rheological parameters */
/* For Poiseuille*/
double nu = 0.;
/*For Bingham*/
double tauy = 0.0;
double mu = 1;

/*For cohesive Bagnold and dry Bagnold*/
double I0 = 0.3 [0];
double mu0 = 0.1 [0];
double deltamu = 0.26 [0];
double dg = 0.4 [1];
double rho = 1.0;
double tauc = 0.0 ; // cohesive stress, =0 if BAGNOLDDRY


/* For turbulent viscous */
double kappa = 0.4;
double rugo = 0.01;
double ll = 0.0;

/** the slope angle is expressed in rad  */
/* up to now constante slope, to be changed */
double slope = 0.25;

/* new parameter for test case with different viscosity*/
int bilayer = 0;
/** 
## Bottom/Top viscosity used in Thomas Algorithm  */

/**In newtonian flow, their values should equal to nu*/
double D_t = 0.;
double D_b = 0.;
/** 
## Conditions at two boundaries (top + bottom) */
(const) vector lambda_b[] = {0.0,0,0}, dub[] = {0,0,0}, u_b[] = {0.0,0,0};
(const) vector lambda_t[] = {-0.0,0,0},u_t[] = {0.0,0,0}, dut[] = {0,0,0};


/** 
# Some functions to compute the equivalent viscosity */
/**Overall, we consider the quantities as point values located at the centers of the layers. When values at an interface are required, they are obtained by averaging the center values of the two adjacent layers.

In particular, as it will be shown in the "shear function", the velocity $u_x$ is treated as a pointwise velocity (vitesse ponctuelle) located at the center of each layer. Strictly speaking, this is not correct, since $u_x$ actually represents the layer-averaged velocity. This approximation introduces a discretization error. However, when the velocity profile is smooth ($C_2$), the resulting error is expected to be second-order accurate.*/


/** ## The shear function (for Bingham, $\mu(I)$ and turbulent viscosities)*/

/** it would be reused by other functions defined later. To ensure the numerical consistency, the numerical scheme is of order two.
*/

double shear(Point point, scalar s, scalar h, int layer, int layercoef){

  double shear=0;
/** Apart from the bottom and top layers, we use the central difference scheme
to compute the shear of each layer. Attention: not the shear at interface, the shear at interface is computed from the average of adjacent layers.
$$
\dfrac{du}{dy}\vert_l = \dfrac{(du/dy)_{l+1}-(du/dy)_{l-1}}{0.5(h_{l-1}+h_{l+1})+h_l}
$$
*/
  if (layercoef>0&&layercoef<(nl-1)){ 
    shear = (s[0,0,layer+1]-s[0,0,layer-1])/(h[0,0,layer]+0.5*(h[0,0,layer+1]+h[0,0,layer-1]));
  }
/** For the layers at boundaries, the center scheme fails. Instead of computing the value of "goste cell" where a profil shape has to be prescribed, we compute the sheat at boudary using  second-order forward/backward difference scheme. This necessitates number of layer $nl\ge 3$, otherwise one should use first order scheme.

To compute the shear of the bottom layer, we write first the Taylor expansion of $u_1$ and $u_2$ at $u_0$
$$
\begin{aligned}
u_1 &= u_0  + \dfrac{du}{dy}\vert_{l=0}\dfrac{1}{2}(h_0+h_1) + \dfrac{1}{2}\dfrac{d^2u}{dy^2}\vert_{l=0}\dfrac{1}{4}(h_0+h_1)^2+O(h^3)\\
u_2 &= u_0  + \dfrac{du}{dy}\vert_{l=0}(h_1+\dfrac{1}{2}(h_0+h_2)) + \dfrac{1}{2}\dfrac{d^2u}{dy^2}\vert_{l=0}(h_1+\dfrac{1}{2}(h_0+h_2))^2+O(h^3)
\end{aligned}
$$
solve it, one gets
$$
\dfrac{du}{dy}\vert_{l=0} = ...+O(h^2)
$$


*/
  else if(layercoef==0){
    // shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));
    shear = -s[0,0,0]*((h[0,0,0]+1.5*h[0,0,1]+0.5*h[0,0,2])/(0.5*(h[0,0,0]+h[0,0,1])*(h[0,0,1]+0.5*h[0,0,0]+0.5*h[0,0,2])))
							+s[0,0,1]*(h[0,0,1]+0.5*h[0,0,0]+0.5*h[0,0,2])/(0.25*(h[0,0,0]+h[0,0,1])*(h[0,0,1]+h[0,0,2]))
							-s[0,0,2]*(0.5*h[0,0,0]+0.5*h[0,0,1])/(0.5*(h[0,0,1]+h[0,0,2])*(h[0,0,1]+0.5*h[0,0,0]+0.5*h[0,0,2]));
  }
  /** Similarily, to compute the shear of the top layer, we write first the Taylor expansion of $u_{nl-2}$ and $u_{nl-3}$ at $u_{nl-1}$
$$
\begin{aligned}
u_{nl-2} &= u_{nl-1}  - \dfrac{du}{dy}\vert_{l=nl-1}\dfrac{1}{2}(h_{nl-1}+h_{nl-2}) + \dfrac{1}{2}\dfrac{d^2u}{dy^2}\vert_{l=nl-1}\dfrac{1}{4}(h_{nl-1}+h_{nl-2})^2+O(h^3)\\
u_{nl-3} &= u_{nl-1}  - \dfrac{du}{dy}\vert_{l=nl-1}(h_{nl-2}+\dfrac{1}{2}(h_{nl-1}+h_{nl-3})) + \dfrac{1}{2}\dfrac{d^2u}{dy^2}\vert_{l=nl-1}(h_{nl-2}+\dfrac{1}{2}(h_{nl-1}+h_{nl-3}))^2+O(h^3)
\end{aligned}
$$
solve it, one gets
$$
\dfrac{du}{dy}\vert_{l=nl-1} = ...+O(h^2)
$$


*/
  else if(layercoef==nl-1){
   // shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
   shear = s[0,0,nl-1]*(h[0,0,nl-1]+1.5*h[0,0,nl-2]+0.5*h[0,0,nl-3])/(0.5*(h[0,0,nl-1]+h[0,0,nl-2])*(h[0,0,nl-2]+0.5*(h[0,0,nl-1]+h[0,0,nl-3])))
					-s[0,0,nl-2]*(h[0,0,nl-2]+0.5*(h[0,0,nl-1]+h[0,0,nl-3]))/(0.25*(h[0,0,nl-1]+h[0,0,nl-2])*(h[0,0,nl-2]+h[0,0,nl-3]))
					+s[0,0,nl-3]*(0.5*(h[0,0,nl-1]+h[0,0,nl-2]))/(0.5*(h[0,0,nl-2]+h[0,0,nl-3])*(h[0,0,nl-2]+0.5*(h[0,0,nl-1]+h[0,0,nl-3])));
    
  } 
/** we choose this kind of regularisation preventing shear from decreasing to zero*/
 // return fabs(shear);
  return sqrt(sq(shear) + .1e-4);  
}

/** ## The pressure field (for $\mu(I)$ rheology)*/
/** For now, we use hydrostatic pressure for the pressure field 
$$
P(y)= \rho g \cos{\theta} (\eta - z_b -y)
$$
Note the $\cos{\theta}$ represents the projection of gravity, and $y\in [0, H]$
*/
double pressionHydro(Point point, scalar h,int layer){

  double H = 0.;
  double zc = 0.;
    for (int l = 0; l < layer; l++) {
    H+=h[0,0,l];
  }
  zc = H + 0.5*h[0,0,layer];
  return rho*G*cos(slope)*(eta[]-zb[]-zc);
}

/** ## The inertial number $I$ (for $\mu(I)$ rheology)*/
/** This is a dimensionless number
$$
I = d\dfrac{du/dy}{\sqrt{P/\rho}},
$$
where $d$ is grain diameter, $P$ is the pressure, $\rho$ is the density.
*/

double nombreInertie(Point point, scalar s, scalar h, int layer){

  double ans;
  ans = dg*shear(point,s,h,layer,layer)/sqrt(pressionHydro(point,h,layer)/rho);
  return ans;
}
/** ## $\mu(I)$ rheology*/
/**
It is shown in litteratue, the rheology of dense granular flow follows Coulombien-like friction law, obeying $\mu(I)$ dependency
$$
\mu(I) = \mu_0 + \dfrac{\Delta \mu}{1+ \frac{I_0}{I}},
$$
where $\mu_0$, $\Delta \mu$ and $I_0$ are rheological parameters, fit from experiments
$$
\tau = {\mu(I)P}=\nu_{eq}\dfrac {\partial u}{\partial y}
$$
for cohesive granualr flows
$$
\tau = \tau_c+ {\mu(I)P}=\nu_{eq}\dfrac {\partial u}{\partial y}
$$
*/
double coeffFrotte(Point point, scalar s, scalar h, int layer){

  double _rapport;
  _rapport = I0/nombreInertie(point, s, h, layer);
  return mu0 + deltamu/(_rapport + 1);

}
/** ## Turbulent flow viscosity (1/2) */
/** A turbulent viscosity is derived from Prandtl mixing theory, permitting to recover the linear + logarithmic velocity profile :
$$\nu_{eq}= \kappa y^2 \frac{\partial u}{\partial y}$$
*** absolute value ?
*/
double NuturbulentPrandtl(Point point, scalar s, scalar h, int layer){

  double nu_eq = 0.0;
  double _y = 0.0;

    for (int l = 0; l < layer; l++) {
    _y+=h[0,0,l];
  }
  _y = _y + 0.5*h[0,0,layer];

  ll = kappa*_y;
  nu_eq = 0.5 + shear(point,s,h,layer,layer)*ll*ll;
  return nu_eq;
}
/** ## Turbulent flow viscosity (2/2) */
/** We adapt the above formula to describe turbulent flows in rivers. 

*** write the exact equation (beware y or y/h or h-y ???)
$$\ell=\sqrt{1-y/h}$$
$$\nu_{eq}= \kappa \ell^2 \frac{\partial u}{\partial y}$$
*/
 
double Nuturbulent(Point point, scalar s, scalar h, int layer){

  double nu_eq = 0.0;
  double _y = 0.0;

    for (int l = 0; l < layer; l++) {
    _y+=h[0,0,l];
  }
  _y = _y + 0.5*h[0,0,layer];

  ll = kappa*(_y+rugo)*sqrt(1-(_y));
  nu_eq = shear(point,s,h,layer,layer)*ll*ll;
  return nu_eq;
}

/**## The equivalent viscosity !*/


double Nueq(Point point, scalar s, scalar h, int layer){

  double ans=0;
/** 
### If dry granular flow
$$
\nu_{eq} = \dfrac{\mu(I)P}{du/dy}
$$
*/  
#if BAGNOLDDRY
  ans =  coeffFrotte(point,s,h,layer)*pressionHydro(point,h,layer)/shear(point,s,h,layer,layer);
/** 
### If cohesive granular flow,
we consider that the shear stress follows a linearized $\mu(I)$ rheology, complement with a cohesive shear stress $\tau_c$, then one gets
$$
\nu_{eq} = \dfrac{(\mu_0 + (\Delta \mu /I_0) I)P + \tau_c}{du/dy}
$$
*/ 
#elif BAGNOLDCOH
  double term1 = 0.;
  double term2 = 0.;
  term1 = deltamu*dg*pow(rho*pressionHydro(point,h,layer),0.5)/I0;
  term2 = (tauc+mu0*pressionHydro(point,h,layer))/shear(point,s,h,layer,layer);
  
  ans = term1 + term2;
/** 
### If Bingham flow
$$
\nu_{eq} = \dfrac{\tau_y}{du/dy} + \mu
$$
*/
#elif BINGHAM
    ans = mu + tauy/shear(point,s,h,layer,layer);
/** 
### If turbulent flows
*/  
#elif TURBULENTPRANDTL
  ans = NuturbulentPrandtl(point,  s, h, layer);
  
#elif TURBULENT
  ans = Nuturbulent(point,  s, h, layer);
/**
Otherwise, a constant viscosity suffices
*/
#else 
  ans = nu;
#endif

  return ans;
}

//with bilayer
// double Nueq(Point point, scalar s, scalar h, int layer, double D, int bilayer){
//   double ans=0;
//   if (bilayer == 0) ans = D;
//   else if (bilayer == 1) ans = coeffFrotte(point,s,h,layer)*pressionHydro(point,h,layer)/(shear(point,s,h,layer,layer)); // other regularization
//   else  ans = mu + tauy/(shear(point,s,h,layer,layer)+dry);
//   return ans;
// }


/** ## Régularisation for $\nu_{eq}$ (should choose another name..) */
/** 
To note that this is the only function being called in diffusionRheo.h, the aim is to compute accordingly the $\nu_{eq}$ for each layer (nueq[nl]), then to be used for vertical diffusion scheme. The reason why this function still exists is that sometimes we want to test another regularisation methode, the critical layer one... 
*/
void regularization(Point point, scalar s, scalar h, double nueq[], double D)
{
  int l;

  //for(l=0 ; l<nl;l++){ // Find the limiting layer shear
   // if (shear(point,s,h,l,l) <=1e-3) {
    //  nlc = l-1;
   // break;
   // }
 // }

  for( l=0 ; l<nl;l++){
    //if ( l<=nlc ) nueq[l] = Nueq(point,s,h,l);
  //  else nueq[l] = nueq[l-1];
    nueq[l] = Nueq(point,s,h,l);
   
  }

  // for(l=0 ; l<nl;l++){ // Potential former bilayer
  //   if ( l<=(nl/2.) )  nueq[l] = Nueq(point,s,h,l, D, bilayer);
  //   else nueq[l] = Nueq(point,s,h,l, D, bilayer+2);
  // }

// Nueq with bilayer need to be tested and compared to analytical solutions
  // for( l=0 ; l<nl;l++){ // Apply regularization
  //   if ( l<=nlc && l <=nl/2 ) nueq[l] = Nueq( point,s,h,l, D, bilayer+1 );
  //   else if ( l<=nlc && l > nl/2 ) nueq[l] = Nueq( point,s,h,l, D, bilayer  );
  //   else nueq[l] = nueq[l-1];
  // }
}
/**
## bibliography

*/