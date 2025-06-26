double HAM = 1e-10;

#define a_VdW(eta,i) (-HAM*(pow(2.*eta[i],-3)-pow(2.*eta[i-1],-3))/Delta )

#include "./hydro-tension.h"
