/**
# Functions $f_s$ and $f_r$ for the FENE-P model

See [log-conform.h](log-conform.h). */

#include "log-conform.h"

double L2 = 1.;

static void fenep (double trA, double * nu, double * eta) {
  *eta = 1;
  *nu = 1./(1. - trA/L2);
} 

event defaults (i = 0) {
  f_s = fenep;
  f_r = fenep;
}

event init (i = 0) {
#if AXI
  double dim = 3;
#else
  double dim = dimension;
#endif  
  scalar trac = trA;
  foreach()
    trac[] = dim*L2/(dim + L2);
}
