/**
# Fixed Flux Phase Change Model

This phase change model does not compute the vaporization
rate from a mass or energy balance, instead, it uses a constant
value *mEvapVal* defined by the user.
*/

double mEvapVal = 0.;

/**
## Fields

We define the field *mEvap* with the vaporization rate
per unit of surface in every interfacial cell, and we
decleare the *mEvapList* used by [evaporation.h](evaporation.h).
*/

scalar mEvap[];
scalar * mEvapList = {mEvap};

/**
## Phase Change Event

We loop over the interfacial cells and we assign the
value of *mEvapVal* in the field *mEvap*.
*/

event phasechange (i++)
{
  foreach() {
    mEvap[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      mEvap[] = mEvapVal;
  }
}

