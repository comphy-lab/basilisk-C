
/**
# Stability functions for the $k$-$\epsilon$ model

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


/**
   We need to find a proper terminology to name these stability functions.The
   GOTM naming convention "cmu_a", "cmu_b", etc... is not fully satisfactory

*/

// default value for c_mu0
double cmu0 = 0.5477;
double cm3_inv;


// cmue_d from gotm with Canuto A parameters
#if KEPS_STAB_D


double cm_N, cm_Nb;
double cm_d0, cm_d1, cm_d2, cm_d3, cm_d4, cm_d5;
double cm_n0, cm_n1, cm_n2, cm_nb0, cm_nb1, cm_nb2;

double anMin, asMax; 

event init(i = 0) {


/**
   These constants are from Canuto et al 2001 (sometimes called the Canuto A
   parameter set). 
*/

  double cc1 = 5.0000;
  double cc2 = 0.8000;
  double cc3 = 1.9680;
  double cc4 = 1.1360;
  double cc5 = 0.0000;
  double cc6 = 0.4000;
  double cb1 = 5.9500;
  double cb2 = 0.6000;
  double cb3 = 1.0000;
  double cb4 = 0.0000;
  double cb5 = 0.3333;
  double cbb = 0.7200;


/**
   These constants are then combined to get the Algebraic stress Model. The
   notations are from Umlauf and Burchard (2005)
*/

  double a1   =  2./3. - cc2/2.;
  double a2   =  1.    - cc3/2.;
  double a3   =  1.    - cc4/2.;
  double a4   =          cc5/2.;
  double a5   =  1./2. - cc6/2.;

  double ab1  =           1. - cb2 ;
  double ab2  =           1. - cb3 ;
  double ab3  =  2. *   ( 1. - cb4);
  double ab4  =  2. *   ( 1. - cb5);
  double ab5  =  2.*cbb*( 1. - cb5);



/**
   Compute numerator and denominator of stability functions.
*/

  cm_N    =   0.5*cc1;
  cm_Nb   =   cb1;
  
  cm_d0   =   36.* cube(cm_N) * sq(cm_Nb);
  cm_d1   =   84.*a5*ab3 * sq(cm_N) * cm_Nb  + 36.*ab5 * cube(cm_N) * cm_Nb;
  cm_d2   =   9.*(sq(ab2)-sq(ab1)) * cube(cm_N) - 12.*(sq(a2)-3.*sq(a3)) * cm_N * sq(cm_Nb);
  cm_d3   =   12.*a5*ab3*(a2*ab1-3.*a3*ab2) * cm_N + 12.*a5*ab3*(sq(a3)-sq(a2)) * cm_Nb      \
    + 12.*ab5*(3.*sq(a3)-sq(a2)) * cm_N * cm_Nb;
  cm_d4   =   48.*sq(a5)*sq(ab3) * cm_N + 36.*a5*ab3*ab5 * sq(cm_N);
  cm_d5   =   3.*(sq(a2)-3.*sq(a3))*(sq(ab1)-sq(ab2)) * cm_N;


  cm_n0   =   36.*a1 * sq(cm_N) * sq(cm_Nb);
  cm_n1   = - 12.*a5*ab3*(ab1+ab2) * sq(cm_N) + 8.*a5*ab3*(6.*a1-a2-3.*a3) * cm_N * cm_Nb   \
    + 36.*a1*ab5 * sq(cm_N) * cm_Nb;
  cm_n2   =   9.*a1*(sq(ab2)-sq(ab1)) * sq(cm_N);

  cm_nb0  =   12.*ab3 * cube(cm_N) * cm_Nb;
  cm_nb1  =   12.*a5*sq(ab3)  * sq(cm_N);
  cm_nb2  =   9.*a1*ab3*(ab1-ab2) * sq(cm_N) + (  6.*a1*(a2-3.*a3)                         \
                                             - 4.*(sq(a2)-3.*sq(a3)) )*ab3 * cm_N * cm_Nb;

  cm3_inv = 1./(cube(cmu0));


  //   mininum value of "an" to insure that "as" > 0 in equilibrium
  double anMinNum  = -(cm_d1 + cm_nb0) + sqrt(sq(cm_d1+cm_nb0) - 4.*cm_d0*(cm_d4 + cm_nb1));
  double anMinDen  = 2.*(cm_d4 + cm_nb1);
  anMin  = anMinNum / anMinDen;

}

/**
   The stability function itself is a fraction with power of $\alpha_N$ and $\alpha_S$.
*/


void stability_function(double an, double as, double c_mu, double c_mup){

  double anLimitFact = 0.5;

  an = max(an,anLimitFact*anMin);

  // compute the equilibrium value of as
  double tmp0  = -cm_d0 - (cm_d1 + cm_nb0)*an - (cm_d4 + cm_nb1)*sq(an);
  double tmp1  = -cm_d2 + cm_n0 + (cm_n1-cm_d3-cm_nb2)*an;
  double tmp2  =  cm_n2 - cm_d5;
    
  as = (abs(tmp2) < 1.e-10)? -tmp0/tmp1: (-tmp1 + sqrt(tmp1*tmp1-4.*tmp0*tmp2) ) / (2.*tmp2);

  // compute stability function
  double dCm  = cm_d0  +  cm_d1*an +  cm_d2*as + cm_d3*an*as + cm_d4*an*an + cm_d5*as*as;
  double nCm  = cm_n0  +  cm_n1*an +  cm_n2*as;
  double nCmp = cm_nb0 + cm_nb1*an + cm_nb2*as;
    
  // divide by cmu^3 only if use the length scale framework.
  c_mu  =  cm3_inv*nCm /dCm;
  c_mup =  cm3_inv*nCmp/dCm;

}


#endif
