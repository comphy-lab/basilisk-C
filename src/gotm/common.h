/**
# A C interface to the General Ocean Turbulence Model (GOTM)

[GOTM](http://gotm.net) is a library implementing a large number of
turbulence models and auxilliary functions adapted to the description
of vertical mixing in ocean models.

This file and other files in this folder provide the C functions
definitions which can be used to call the corresponding Fortran
functions in the GOTM libraries.

Note that these files are entirely independent from Basilisk.

All the header files (except this one) were generated automatically
from the Fortran 90 source code of GOTM (stable version 5.2.1). The
[Makefile]() provided can be used to regenerate them if necessary.

## Installation of GOTM

The GOTM libraries must be installed separately. Only version 5.2.1 is
supported at the moment. A [patch](gotm.patch) must be applied before
compilation.

See the [GOTM installation
instructions](https://gotm.net/software/linux/) for details but
otherwise follow this recipe:

~~~bash
wget https://github.com/gotm-model/code/archive/v5.2.1.tar.gz
tar xzvf v5.2.1.tar.gz
cd code-5.2.1/src
wget http://basilisk.fr/src/gotm/gotm.patch?raw -O gotm.patch
patch -p0 < gotm.patch 
cd ..
mkdir build
cd build
cmake ../src -DGOTM_USE_FABM=off
make
~~~
*/

typedef int integer;
typedef int logical;
typedef double realtype;
typedef long timestepkind;
typedef struct {
  realtype * a;
} realtype_1d;

int strlencheck (const char * s) {
  return (s != NULL ? strlen(s) : 0);
}

extern integer __gotm_MOD_nlev;
#define gotm_nlev __gotm_MOD_nlev
extern realtype __gotm_MOD_dt;
#define gotm_dt __gotm_MOD_dt
extern realtype __gotm_MOD_cnpar;
#define gotm_cnpar __gotm_MOD_cnpar
extern integer __gotm_MOD_buoy_method;
#define gotm_buoy_method __gotm_MOD_buoy_method

#define airsea_tx __airsea_MOD_tx
#define airsea_ty __airsea_MOD_ty
#define airsea_evap __airsea_MOD_evap

#define time_unixseconds() ((long)((time_julianday - 2440587.5)*86400. + \
				   time_secondsofday))
