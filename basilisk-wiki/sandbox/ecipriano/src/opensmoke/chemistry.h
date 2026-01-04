#include "opensmoke.h"
#include "common-chemistry.h"
#include "reactors.h"

void opensmoke_ode_solver (ode_function batch,
    unsigned int NEQ, double dt, double * y0, void * data)
{
  OpenSMOKE_ODESolver (batch, NEQ, dt, y0, data);
}

event defaults (i = 0) {
  OpenSMOKE_InitODESolver();
  stiff_ode_solver = opensmoke_ode_solver;
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver();
}

