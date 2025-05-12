/**
# Higher-order neo-Hookean law for elastic membranes

In this file we define the higher-order neo-Hookean strain energy function, used
to compute the elastic force on each Lagrangian node of the membrane.
*/

#ifndef E_S
  #define E_S 1.
#endif
#ifndef C_3
  #define C_3 (E_S/30.)
#endif

#define DWDL1(L1, L2) ((L1 - 1./(cube(L1)*sq(L2)))*(E_S + 2*C_3*\
  sq(sq(L1) + sq(L2) + 1./(sq(L1*L2)) - 3)))
#define DWDL2(L1, L2) ((L2 - 1./(cube(L2)*sq(L1)))*(E_S + 2*C_3*\
  sq(sq(L1) + sq(L2) + 1./(sq(L1*L2)) - 3)))

#include "elasticity-ft.h"
