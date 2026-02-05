/**
## Default OpenSMOKE++ header file
Need the OpenSMOKE++ interface to link the library functions see 
[OpenSMOKE++ Interface](https://github.com/edocipriano/OpenSMOKEppInterface) 
and [OpenSMOKE++](https://www.opensmokepp.polimi.it/).

This library provides functions to read thermodynamic, transport properties and reaction rates, 
all key necessity for computing reacting flows.
*/

#include "OpenSMOKE_Interface.h"

#define OPENSMOKE 1

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

char* kinfolder;

event cleanup (t = end)
{
  OpenSMOKE_Clean ();
}
