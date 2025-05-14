/**
# Equation of state for the ocean

This module implements several equations of state for the ocean: a linear
equation of state and EOS-80

Units for temperatures are in degC and Salinity are in PSU (g/kg)

The user can use 
#define EOS_LIN 1 
for a linear equation of state
#define EOS_UNESCO 1
for EOS-80


## Linear equation of state

The linear equation of state is given by

$$\rho = \rho_0 (1 - \alpha_T (T - T_0) + \beta_S (S - S_0))$$

where $\rho_0$, T_0$, $S_0$ are given constants of a reference state. $\alpha_T$
and $\beta_S$ are the thermal expansion coefficient and haline contraction
coefficient defined as

$$
\alpha_T ={\frac {1}{\rho_0 }}{\frac {\partial \rho }{\partial T}}{\Bigg |}_{T ,p}
$$


$$
\beta_S ={\frac {1}{\rho_0 }}{\frac {\partial \rho }{\partial S}}{\Bigg |}_{T ,p}
$$

The values below are commonly used values in oceanography.
*/

#if EOS_LINEAR

double alphaT = 2.e-4; // K^-1
double betaS  = 8.e-4; // psu^-1
double T_ref = 0.;     // degreeC
double S_ref = 0.;     // psu
double rho_ref   = 1000.; // kg/m^3


/**
   Density and buoyancy
*/

#define eos_rho(T,S) (rho_ref*(1 - alphaT*(T - T_ref) + betaS*(S - S_ref)))

#define eos_b(T,S) (G*(rho_ref - eos_rho(T,S))/rho_ref)


/**
   Thermal expansion coefficient: alpha = -1/rho_0 drho/dT
   and Haline contraction coefficient: beta = 1/rho_0 drho/dS

   For the linear equation of state, these functions are trivial. We still keep
   it for the general case.
*/


#define eos_alpha(T,S) (alphaT)

#define eos_beta(T,S) (betaS)
  

#elif EOS_UNESCO
/**
##   EOS-80: Sea Water equation of state

Unesco (1983); (Fofonoff and Millard, 1983).
https://www.jodc.go.jp/jodcweb/info/ioc_doc/UNESCO_tech/059832eb.pdf

*/


double eos_QR  = 999.842594;
double eos_Q01 = 6.793952e-2;
double eos_Q02 = -9.095290e-3; 
double eos_Q03 = 1.001685e-4;
double eos_Q04 = -1.120083e-6;
double eos_Q05 = 6.536332e-9; 
double eos_Q10 = 0.824493;
double eos_Q11 = -4.08990e-3;
double eos_Q12 = 7.64380e-5; 
double eos_Q13 = -8.24670e-7;
double eos_Q14 = 5.38750e-9;
double eos_QS0 = -5.72466e-3; 
double eos_QS1 = 1.02270e-4;
double eos_QS2 = -1.65460e-6;
double eos_Q20 = 4.8314e-4; 


/**
   eos_rho: compute rho at P = 0
 */

#define eos_rho(T,S) (eos_QR                                          \
                        + (T)*(eos_Q01 + (T)*(eos_Q02 + (T)*(eos_Q03 + (T)*(eos_Q04 + (T)*eos_Q05)))) \
                        + (S)*(eos_Q10 + (T)*(eos_Q11 + (T)*(eos_Q12 + (T)*(eos_Q13 + (T)*eos_Q14))) \
                        + sqrt(S)*(eos_QS0 + (T)*(eos_QS1 + (T)*eos_QS2)) + (S)*eos_Q20))

/**
   eos_b: compute buoyancy at P = 0: b = g*(rho0 - rho)/rho0
 */

#define eos_b(T,S) (G*(eos_QR - eos_rho(T,S))/eos_QR)

/**
   Thermal expansion coefficient: alphaT = -1/rho_0 drho/dT
*/

#define eos_alpha(T,S) (-1/eos_QR*                                     \
                          (eos_Q01 + (T)*(2*eos_Q02 + (T)*(3*eos_Q03 + (T)*(4*eos_Q04 + (T)*5*eos_Q05))) \
                           + (S)*(eos_Q11 + (T)*(2*eos_Q12 + (T)*(3*eos_Q13 + (T)*4*eos_Q14)) + sqrt(S)*(eos_QS1+(T)*2*eos_QS2))))

/**
   Haline contraction coefficient: betaS = 1/rho_0 drho/dS
*/

#define eos_beta(T,S) (1/eos_QR*                                      \
                         (eos_Q10 + (T)*(eos_Q11 + (T)*(2*eos_Q12 + (T)*(3*eos_Q13 + (T)*4*eos_Q14))) \
                          + sqrt(S)*(1.5*eos_QS0 + (T)*(1.5*eos_QS1 + (T)*1.5*eos_QS2)) + (S)*2*eos_Q20))

#endif
