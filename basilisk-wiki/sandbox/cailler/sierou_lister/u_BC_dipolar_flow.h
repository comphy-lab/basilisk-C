/** 
# Far-field velocities in the presence of a *dipolar flow*

From [Sierou \& Lister, (2004)](#sierou2004), *Appendix A*, the consideration 
of a far-field *dipolar flow* can be described at dominant order by:

$$
\mu_d(R) = \sqrt{\dfrac{\sigma}{\rho}} \, \widetilde{\mu}_0 \, {R}^{1/2}
$$

with $\sigma$ the surface tension, $\rho$ the liquid density, 
$\widetilde{\mu}_0$ the dipolar flow intensity, and $R = \sqrt{r^2 + z^2}$ 
the radial distance in spherical-polar coordinates to the cone apex (see 
figure below).

![Recoiling cone in 3D--AXI](img_sierou/sierou_recoil.png){width="20%"}

Tedious theoretical developments have shown that velocity potentials 
interior to the liquid cone $\phi^-$ and outside of it $\phi^+$ are directly 
proportional to the far-field dipolar intensity $\widetilde{\mu}_0$:

$$
\phi^-(R,\theta,t) = A_0^{-} \, P_{1/2}(\cos \theta) \, R^{1/2} 
\quad ; \quad 
\phi^+(R,\theta,t) = A_{0,Q}^{+} \, Q_{1/2}(\cos \theta) \, R^{1/2}
$$

where $\theta = \arctan(r/z)$ and:
$$
\left\{
  \begin{array} {rcl}
    {A_0}^- 
    &=& 
    \sqrt{\dfrac{\sigma}{\rho}} \, 
    \dfrac{\cos \theta_0 \, Q_{1/2}(\cos \theta_0) - Q_{3/2}(\cos \theta_0)}
      {P_{3/2}(\cos \theta_0) Q_{1/2}(\cos \theta_0) 
        - P_{1/2}(\cos \theta_0) Q_{3/2}(\cos \theta_0)} 
        \, \widetilde{\mu}_0 \\
    & &\\
    {A_{0,Q}}^+ 
    &=& 
    \sqrt{\dfrac{\sigma}{\rho}} \,
    \dfrac{\cos \theta_0 \, P_{1/2}(\cos \theta_0) - P_{3/2}(\cos \theta_0)}
    {P_{3/2}(\cos \theta_0) Q_{1/2}(\cos \theta_0) 
      - P_{1/2}(\cos \theta_0) Q_{3/2}(\cos \theta_0)} 
      \, \widetilde{\mu}_0
  \end{array}
\right.
$$


with $P_{i/2}$ and $Q_{i/2}$ the *Legendre* polynomials of the first and second 
kind of fractional degree $i/2$. 
These special functions are defined in the following file:
*/

#include "elliptic.h"

/** 
Then we are able to compute the coefficients of $\phi^-$ and $\phi^+$: 
*/

#define PHI_0_A phi_0_A(THETA_0, mu_0[k])
#define PHI_0_OQ phi_0_OQ(THETA_0, mu_0[k])

/** 
where `mu_0[k]` is a user-defined constant defined in the user file to model 
the dipolar flow intensity.

Now, the far-field velocities in each phase
$\mathbf{u}_\infty^\pm = \boldsymbol{\nabla} \phi^\pm$ are readily calculated:

$$
\mathbf{u}_{\infty}^-(R,\theta,t) 
= \begin{pmatrix}
    A_0^- \, P_{1/2} (\cos \theta)/2 \\
    A_0^- \, \partial_\theta P_{1/2} (\cos \theta)
\end{pmatrix} R^{-1/2}
\,\, ; \,\, 
\mathbf{u}_{\infty}^+(R,\theta,t) 
= \begin{pmatrix}
    {A_{0,Q}}^+ \, Q_{1/2} (\cos \theta)/2 \\
    {A_{0,Q}}^+ \, \partial_\theta Q_{1/2} (\cos \theta)
\end{pmatrix} R^{-1/2}
$$

To pass them into the axisymmetric cylindrical coordinate system of `Basilisk`, 
one needs to project $\mathbf{u}_\infty^\pm$ onto $(\mathbf{e}_z, \mathbf{e}_r)$: 

$$
u_{z,\infty}^{\pm} = \cos \theta \, u_{R,\infty}^{\pm} 
  - \sin \theta \, u_{\theta,\infty}^{\pm} 
\quad ; \quad
u_{r,\infty}^{\pm} = \sin \theta \, u_{R,\infty}^{\pm} 
  + \cos \theta \, u_{\theta,\infty}^{\pm}
$$

Along with the previous formulae, we also compute derivatives 
used in the definition of the *Neumann* pressure boundary condition that 
takes into account the far-field velocity with non-zero normal components:
*/


                            /* LIQUID PHASE */

double uxA (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*cos(th)*legendreP_half(th);
  double B = sin(th)*legendreP_prime_half(th);

  return ( PHI_0_A*(A - B)/sqrt(R) );
}

double uyA (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*sin(th)*legendreP_half(th);
  double B = cos(th)*legendreP_prime_half(th);

  return ( PHI_0_A*(A + B)/sqrt(R) );
}

// ∂uxA/∂x
double duxA_dx (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 8.*sq(R)*legendreP_half(th);
  double B = 24.*x*R*legendreP_three_half(th);
  double C = 15.*sq(R)*legendreP_five_half(th);

  return ( 0.25*PHI_0_A*( A - B + C )/pow(R, 7./2.) );
}

// ∂uxA/∂y
double duxA_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 4.*x*sq(R)*legendreP_half(th);
  double B = R*( 9.*sq(x) + sq(y) )*legendreP_three_half(th);
  double C = 5.*x*sq(R)*legendreP_five_half(th);

  return ( -0.75*PHI_0_A*( A - B + C )/( y*pow(R, 7./2.)) );
}

// ∂uyA/∂x
double duyA_dx (double x, double y){
  return ( duxA_dy(x,y) );
}

// ∂uyA/∂y
double duyA_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = ( 15.*pow(x, 4.) + 14.*sq(x)*sq(y) - pow(y, 4.) )*legendreP_half(th);
  double B = -2.*R*( 5.*sq(x) + sq(y) )*legendreP_three_half(th);
  double C = 5.*x*sq(R)*legendreP_five_half(th);

  return ( 0.25*PHI_0_A*( A + 3.*x*( B + C ) )/( sq(y)* pow(R, 7./2.)) );
}


                              /* GAS PHASE */

double uxO (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*cos(th)*legendreQ_half(th);
  double B = sin(th)*legendreQ_prime_half(th);

  return ( PHI_0_OQ*(A - B)/sqrt(R) );
}

double uyO (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);
  double A = 0.5*sin(th)*legendreQ_half(th);
  double B = cos(th)*legendreQ_prime_half(th);

  return ( PHI_0_OQ*(A + B)/sqrt(R) );
}

// ∂uxO/∂x
double duxO_dx (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 8.*sq(R)*legendreQ_half(th);
  double B = 24.*x*R*legendreQ_three_half(th);
  double C = 15.*sq(R)*legendreQ_five_half(th);

  return ( 0.25*PHI_0_OQ*( A - B + C )/pow(R, 7./2.) );
}

// ∂uxO/∂y
double duxO_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = 4.*x*sq(R)*legendreQ_half(th);
  double B = R*( 9.*sq(x) + sq(y) )*legendreQ_three_half(th);
  double C = 5.*x*sq(R)*legendreQ_five_half(th);

  return ( -0.75*PHI_0_OQ*( A - B + C )/( y*pow(R, 7./2.)) );
}

// ∂uyO/∂x
double duyO_dx (double x, double y){
  return ( duxO_dy(x,y) );
}

// ∂uyO/∂y
double duyO_dy (double x, double y){
  double R = sqrt( sq(x) + sq(y) );
  double th = atan2(y,x);  
  double A = ( 15.*pow(x, 4.) + 14.*sq(x)*sq(y) - pow(y, 4.) )*legendreQ_half(th);
  double B = -2.*R*( 5.*sq(x) + sq(y) )*legendreQ_three_half(th);
  double C = 5.*x*sq(R)*legendreQ_five_half(th);

  return ( 0.25*PHI_0_OQ*( A + 3.*x*( B + C ) )/( sq(y)* pow(R, 7./2.)) );
}



/**
## References
~~~bib
@article{sierou2004,
  author = {Sierou, A. and Lister, J. R.},
  title = {Self-similar recoil of inviscid drops},
  journal = {Physics of Fluids},
  volume = {16},
  number = {5},
  pages = {1379-1394},
  year = {2004},
  month = {05},
  issn = {1070-6631},
  doi = {10.1063/1.1689031},
  url = {https://doi.org/10.1063/1.1689031},
  eprint = {https://pubs.aip.org/aip/pof/article-pdf/16/5/1379/19153172/1379\_1\_online.pdf},
}
~~~
*/