
/**
# $k$-$\epsilon$ implementation based on Umlauf and Burchard (2005)

This version of the $k$-$\epsilon$ model is taken from the GOTM implementation.


We use the $k$-$\epsilon$ approach to compute the eddy viscosity $\nu_t$ and
diffusivity $\nu'_t$. We have two prognostic equations for turbulent kinetic
energy ($k$) and dissipation ($\epsilon$). The equation for the TKE $k\equiv
\overline{ \mathbf{u'}^{2}}/2$ can be written in the form 


$$\frac{\partial k}{\partial t}+\frac{\partial T_k}{\partial z}=P+B-\epsilon $$

where the eddy flux $T_k$, the production rate $P$, and the dissipation rate
$\epsilon$ can be expressed exactly from correlations of turbulent
fluctuations. The flux $T_k$ tends to smooth the inhomogeneities of the
production $P$ while preserving the integral of $k$. It is modelled by a
diffusion law $T_k=-\nu_t \partial k/\partial z $. $P$ and $B$ are described below.



The $k$-$\epsilon$ approach considers $\epsilon$ as the additional variable
instead of the turbulence scale $l$. It takes a similar form as the TKE
equation,


$$ \frac{\partial\epsilon}{\partial t}+\frac{\partial}{\partial z}\left(\frac{\nu_{t}}{\sigma_{\epsilon}}\frac{\partial \epsilon}{\partial z}\right)=\frac{\epsilon}{k}\left(c_{\epsilon 1}\,P+c_{\epsilon 3}\,\overline{w'b'}-c_{\epsilon 2}\,\epsilon\right) $$

This is a pseudo-empirical equation, and the values chosen for the coefficients
are $\sigma_\epsilon = 1.3$, $c_{\epsilon 1} = 1.44$ and $c_{\epsilon 2} =
1.92$. Those are the classical choices for neutrally buoyant fluids, resulting
from fits to simple turbulent flow configurations \citep{rodi_examples_1987}.

*/

#define KEPS_STAB_D 1
#include "k-epsilon_stability.h"

scalar tke, epsilon;
scalar nuh, num;
int nlp1;

event defaults (i = 0)
{
  nlp1 = nl + 1;
  tke = new scalar[nlp1];
  epsilon = new scalar[nlp1];
  nuh = new scalar[nlp1];
  num = new scalar[nlp1];

}

event cleanup (i = end) {
  delete ({tke, epsilon});
  delete ({nuh, num});
}

void vertical_mixing()
{

/**
   I think ideally we want to keep this routine as a 1d process: the main
   foreach should encapsulate the entire vertical mixing.
*/



  foreach() {


/**
   Declaration of 1d fields.  

   TODO: this is probably not consistent with the existing vertical diffusivity
   solver.

   This mostly depends if we re-use the vertical diffusivity solver or if we
   write a 1d vertical solver.
*/

  // save u and v (only needed to compute the shear)

/**
   TODO: uequation, vequation, temperature and salinity in GOTM.

   This looks very much like vertical_diffusion  in src/layered/diffusion.h

   in GOTM, this the implicit step is done with diff_center

   but with a couple of differences: 
   - the diffusivity coefficient is function of height.
   - the diffusivity coeff is defined at cell interface
   - double check boundary conditions
   - use of so called Patankar trick for the drag
   - there are source terms in the rhs of

  // time step for u and v
  // time step for T and S    
    

*/

/**
##   Compute shear and stratification at layer interface
   
   $$ S^2 = \left( \frac{\partial u}{\partial z} \right)^2 + \left( \frac{\partial v}{\partial z} \right)^2 $$
   
   This is the ad-hoc discretization of Burchard (2002).
   https://www.sciencedirect.com/science/article/pii/S1463500302000094

   $$N^2_T = g \alpha_T \frac{\partial T}{\partial z}$$

   $$N^2_S = g \alpha_S \frac{\partial S}{\partial z}$$

   $$N^2 = N^2_T + N^2_S$$

*/

    double SS[nlp1];
    double NN[nlp1];
    double NNT[nlp1];
    double NNS[nlp1];
    
    for (int l = 1; l < nl; l++) {

      double dz = 0.5*(h[0,0,l] + h[0,0,l-1]);
      double du =  u.x[0,0,l] - u.x[0,0,l-1];
      double dv =  u.y[0,0,l] - u.y[0,0,l-1];
      double dT =  T[0,0,l] - T[0,0,l-1];
      double dS =  S[0,0,l] - S[0,0,l-1];
      double Tf =  0.5*(T[0,0,l] + T[0,0,l-1]);
      double Sf =  0.5*(S[0,0,l] + S[0,0,l-1]);

      SS[l] = sq(du/dz) + sq(dv/dz);

      NNT[l] = G*eos_alpha(Tf,Sf)*dT/dz;
      NNS[l] = -G*eos_beta(Tf,Sf)*dS/dz;

      NN[l] = NNT[l] + NNS[l];
    }

    // fix upper and lower faces
    SS[0] = SS[1];
    SS[nl] = SS[nl-1];

    NNT[0]  = 0.;
    NNT[nl] = 0.;
    NNS[0]  = 0.;
    NNS[nl] = 0.;
    NN[0]   = 0.;
    NN[nl] = 0.;


/**
##   Compute RHS of TKE and dissipation equation:

   Production and conversion from TKE to Potential Energy

$$   P = - \overline{u'w'}\frac{\partial u}{\partial z} - \overline{v'w'}\frac{\partial v}{\partial z} = \nu_{t}\left|\frac{\partial\mathbf{u}}{\partial z}\right|^{2} $$



The conversion rate to potential energy $\overline{w'b'}$ also represents the
vertical flux of buoyancy. Within the eddy diffusivity hypothesis, it is
expressed as 

$$ \overline{w'b'}=-\nu'_t \frac {\partial b}{\partial z}\equiv -\nu'_t N^2 $$

with $N^2 = \partial b / \partial z$ the local buoyancy frequency.

*/

    double Prod[nlp1];
    double Bconv[nlp1];
    
    for (int l = 0; l < nl+1; l++) {
      Prod [l] = num[l]*SS[l];
      Bconv[l] = -nuh[l]*NN[l];
    }
    

/**
##   Compute the stability functions

The turbulent viscosity and turbulent diffusivities are computed as 

$$    \nu_t = c_{\mu}\frac{k^2}{\epsilon} \qquad \nu'_t = c'_{\mu}\frac{k^2}{\epsilon} $$

For homogeneous fluids, the standard values $c_\mu\simeq c'_\mu \simeq 0.09$ are
used to reproduce data from simple shear flows.  In the stratified case, they
are assumed to depend on the local shear and stratification, expressed as
functions of the respective non-dimensional parameters $\alpha_S=k^2
\left|{\partial\mathbf{u}}/{\partial z}\right|^{2} /\epsilon^2$ and
$\alpha_N=k^2N^2 /\epsilon^2$. These so-called 'stability functions' fulfil
appropriate physical and mathematical constraints, as discussed by
\cite{burchard_comparative_2001} and \cite{umlauf2005second}.  Note that $\sqrt
{\alpha_N}$ can be viewed as the inverse of a turbulent Froude number from the
scaling $\epsilon\sim k^{3/2}/l$ , which yields $\sqrt{\alpha_N}=N
l/k^{1/2}$. Moreover, the Richardson number is the ratio $Ri=\alpha_N/\alpha_S$.


*/

    double alpha_N[nlp1];
    double alpha_S[nlp1];
    double c_mu[nlp1];
    double c_mup[nlp1];

    for (int l = 0; l < nl+1; l++) {
      double tau2 = sq(tke[0,0,l]/epsilon[0,0,l]);
      alpha_N[l] = tau2 * SS[l];
      alpha_S[l] = tau2 * NN[l];

      // clip negative values for shear only
      alpha_S[l] = max(alpha_S[l], 1e-10);

    }

/**
   There are several versions of the stability functions.
*/

    for (int l = 1; l < nl; l++) {
      stability_function(alpha_N[l], alpha_S[l], c_mu[l], c_mup[l]);
    }
    

/**
   Step turbulent quantities: tke and epsilon.  

!! These are coupled equations. especially the epsilon equation is full of
   $\epsilon/k$.

*/

/**
   TODO: tkeeq, dissipationeq in GOTM.

   This looks very much like vertical_diffusion  in src/layered/diffusion.h

   in GOTM, this the implicit step is done with diff_face

   but with a couple of differences: 
   - turbulent fields are defined at cell interface (instead of cell center for diffusion.h)
   - the diffusivity coefficient is function of height.
   - the diffusivity coeff is defined at cell interface
   - double check boundary conditions
   - use of so called Patankar trick for the production term (positivity)
   - there are source terms in the rhs of    

*/



  // call do_tke(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)
  // call do_kb(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)
  // call do_lengthscale(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,NN,SS)
  // call do_epsb(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)


/**
   Compute turbulent viscosity and turbulent diffusivity coefficients
*/
    
    for (int l = 0; l < nl+1; l++) {
      
      double dissip_len = cube(cmu0)*sqrt(cube(tke[0,0,l]))/epsilon[0,0,l];
      double dim_visc = sqrt(tke[0,0,l])*dissip_len;
      
      num[0,0,l] = dim_visc*c_mu[l];
      nuh[0,0,l] = dim_visc*c_mup[l];
      
    }
    
  }// end of 1d model
}

