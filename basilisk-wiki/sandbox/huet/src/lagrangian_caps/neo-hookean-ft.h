/**
# Neo-Hookean law for elastic membranes

In this file we define the neo-Hookean strain energy function, used to compute
the elastic force on each Lagrangian node of the membrane.
*/

#ifndef E_S
  #define E_S 1.
#endif

#define DWDL1(L1, L2) (E_S/(3.*L1)*(sq(L1) - 1./(sq(L1*L2))))
#define DWDL2(L1, L2) (E_S/(3.*L2)*(sq(L2) - 1./(sq(L1*L2))))

#include "elasticity-ft.h"
