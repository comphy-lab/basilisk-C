/**
   EOS-80: Sea Water equation of state

   Unesco (1983); (Fofonoff and Millard, 1983).
   https://www.jodc.go.jp/jodcweb/info/ioc_doc/UNESCO_tech/059832eb.pdf
   
   Temperatures are in degC and Salinity are in PSU (g/kg)
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

double eos_Gravity = 9.81;

/**
   eos_rho_s: compute rho at P = 0
 */

#define eos_rho_s(T,S) (eos_QR                                          \
                        + (T)*(eos_Q01 + (T)*(eos_Q02 + (T)*(eos_Q03 + (T)*(eos_Q04 + (T)*eos_Q05)))) \
                        + (S)*(eos_Q10 + (T)*(eos_Q11 + (T)*(eos_Q12 + (T)*(eos_Q13 + (T)*eos_Q14))) \
                        + sqrt(S)*(eos_QS0 + (T)*(eos_QS1 + (T)*eos_QS2)) + (S)*eos_Q20))

/**
   eos_b_s: compute buoyancy at P = 0: b = g*(rho0 - rho)/rho0
 */

#define eos_b_s(T,S) (eos_Gravity*(eos_QR - eos_rho_s(T,S))/eos_QR)

/**
   Thermal expansion coefficient: alpha = -1/rho_0 drho/dT
*/

#define eos_alpha_s(T,S) (-1/eos_QR*                                     \
                          (eos_Q01 + (T)*(2*eos_Q02 + (T)*(3*eos_Q03 + (T)*(4*eos_Q04 + (T)*5*eos_Q05))) \
                           + (S)*(eos_Q11 + (T)*(2*eos_Q12 + (T)*(3*eos_Q13 + (T)*4*eos_Q14)) + sqrt(S)*(eos_QS1+(T)*2*eos_QS2))))

/**
   Haline expansion coefficient: beta = 1/rho_0 drho/dS
*/

#define eos_beta_s(T,S) (1/eos_QR*                                      \
                         (eos_Q10 + (T)*(eos_Q11 + (T)*(2*eos_Q12 + (T)*(3*eos_Q13 + (T)*4*eos_Q14))) \
                          + sqrt(S)*(1.5*eos_QS0 + (T)*(1.5*eos_QS1 + (T)*1.5*eos_QS2)) + (S)*2*eos_Q20))
