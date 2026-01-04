#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef MAXLEVEL
# define MAXLEVEL 10
#endif

#ifndef YTOL
# define YTOL 1.e-1
#endif

#ifndef TTOL
# define TTOL 1.e0
#endif

#ifndef UTOL
# define UTOL 1.e-1
#endif

#ifndef COMBUSTION
# define COMBUSTION 1
#endif

#ifndef SIGMA
# define SIGMA 0.03
#endif

#ifndef EMISSIVITY
# define EMISSIVITY 0.0
#endif

#ifndef TEMPERATURE
# define TEMPERATURE 300.
#endif

#ifndef TEMPERATURE_DROPLET
# define TEMPERATURE_DROPLET 300.
#endif

#ifndef PRESSURE
# define PRESSURE 1.
#endif

#ifndef DIAMETER
# define DIAMETER 1e-3
#endif

#ifndef KINFOLDER
# define KINFOLDER two-step/methanol
#endif

#ifndef LIQFOLDER
# define LIQFOLDER LiquidProperties
#endif

#ifndef LIQUID
# define LIQUID CH3OH_1.
#endif

#ifndef GAS
# define GAS N2_0.79_O2_0.21
#endif

#ifndef GRAVITY
# define GRAVITY 0.
#endif

#ifndef FIBER
# define FIBER 0.1
#endif

#ifndef ENVIRONMENT
# define ENVIRONMENT 7.986462e+01
#endif

#ifndef MAX_DD02
# define MAX_DD02 0.005
#endif

#ifndef SPARK_DIAMETER
# define SPARK_DIAMETER 0.2
#endif

#ifndef SPARK_START
# define SPARK_START 0.0
#endif

#ifndef SPARK_TIME
# define SPARK_TIME 0.01
#endif

#ifndef SPARK_VALUE
# define SPARK_VALUE 1e7
#endif

#ifndef DUMP_EVERY
# define DUMP_EVERY 0.005
#endif

#ifndef MOVIE_EVERY
# define MOVIE_EVERY 0.001
#endif

