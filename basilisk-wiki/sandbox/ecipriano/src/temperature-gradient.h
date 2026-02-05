#include "intgrad.h"

scalar mEvap[], * mEvapList = {mEvap};

// fixme: put these operations in a function, such that the used can easily move
// it from the phase change event
event phasechange (i++) {
  phase_tracers_to_scalars (liq, f, tol = F_ERR);
  phase_tracers_to_scalars (gas, f, tol = F_ERR);

  scalar fl[], fg[];
  foreach() {
    fl[] = f[];
    fg[] = 1. - f[];
  }

  face vector fsl[], fsg[];
  face_fraction (fl, fsl);
  face_fraction (fg, fsg);

  //// fixme as:
  //// face vector fsl[], fsg[];
  //// face_fraction (f, fsl, inverse = false);
  //// face_fraction (f, fsg, inverse = true);

  foreach_interfacial (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TIntVal, false);
    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true, TIntVal, false);

    scalar lambdal = liq->lambda, lambdag = gas->lambda;
    scalar deltahev = liq->dhev;

    mEvap[] += (deltahev[] > 0.) ? lambdal[]*ltrgrad/deltahev[] : 0.;
    mEvap[] += (deltahev[] > 0.) ? lambdag[]*gtrgrad/deltahev[] : 0.;
  }

  foreach_interfacial_plic (f, F_ERR) {
    scalar TL = liq->T, TG = gas->T;
    double ltrgrad = ebmgrad (point, TL, fl, fg, fsl, fsg, false, TIntVal, false);
    double gtrgrad = ebmgrad (point, TG, fl, fg, fsl, fsg, true, TIntVal, false);

    scalar slT = liq->STexp, sgT = gas->STexp;
    scalar lambdal = liq->lambda, lambdag = gas->lambda;

    slT[] += lambdal[]*ltrgrad*dirac;
    sgT[] += lambdag[]*gtrgrad*dirac;
  }

  phase_scalars_to_tracers (liq, f);
  phase_scalars_to_tracers (gas, f);
}

