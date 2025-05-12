/**
# OpenSMOKE++

Basilisk interface for the [OpenSMOKE++ library](#cuoci2015opensmoke++),
which is a C++ framework for numerical simulations of reacting systems
with detailed kinetic mechanisms, including thousands of chemical
species and reactions. It provides the functionalities for managing
the calculation of thermodynamics and transport properties of mixtures
and for handling the chemical kinetics.

## OpenSMOKEInterface
This module requires an external library, which contains both the
OpenSMOKE++ library itself, and the C interface which contains the
functions that can be called from Basilisk. It can be downloaded
from [github](https://github.com/edocipriano/OpenSMOKEppInterface)
and it can be easily compiled using cmake, provided the following
external dependencies: 
[OpenSMOKE++](#cuoci2015opensmoke++), 
[boost](https://www.boost.org/), 
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), 
[libconfig++](https://github.com/hyperrealm/libconfig).
*/

#include "OpenSMOKE_Interface.h"

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

/**
## User Data

The variable *kinfolder*, with the path of the kinetic scheme folder
relative to the path $OPENSMOKE_INTERFACE/kinetics, must be set by the
user.
*/

char * kinfolder, * liqfolder;

event defaults (i = 0)
{
  char kinfolder_root[120];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  char liqfolder_root[120];
  sprintf (liqfolder_root, "%s/kinetics/LiquidProperties/%s",
      getenv ("OPENSMOKE_INTERFACE"), liqfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  OpenSMOKE_ReadLiquidKinetics (kinfolder_root);
  OpenSMOKE_ReadLiquidProperties (liqfolder_root);
}

event cleanup (t = end)
{
  OpenSMOKE_Clean ();
}

/**
## References

~~~bib
@article{cuoci2015opensmoke++,
  title={OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms},
  author={Cuoci, Alberto and Frassoldati, Alessio and Faravelli, Tiziano and Ranzi, Eliseo},
  journal={Computer Physics Communications},
  volume={192},
  pages={237--264},
  year={2015},
  publisher={Elsevier}
}
~~~
*/

