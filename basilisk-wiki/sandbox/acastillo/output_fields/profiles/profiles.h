#define _mindelta (L0 / (double)(1 << (grid->maxdepth + 1)))

#if dimension == 2
#define SETUP_PROFILE_PLANE(hprof)\
  coord box[2], nsamples;\
  box[0] = (coord){xmin, hprof};\
  box[1] = (coord){xmax, hprof};\
  nsamples = (coord){m1, 1};
#else
#define SETUP_PROFILE_PLANE(hprof)\
  coord box[2], nsamples;\
  box[0] = (coord){xmin, ymin, hprof};\
  box[1] = (coord){xmax, ymax, hprof};\
  nsamples = (coord){m1, m2, 1};
#endif

#define PROFILE_PARAMS \
  scalar w = unity, \
  const char * filename = "profiles.asc", \
  double hmin = Z0 + _mindelta/2., \
  double hmax = Z0 + L0 - _mindelta/2., \
  double xmin = X0, \
  double xmax = X0 + L0, \
  double ymin = Y0, \
  double ymax = Y0 + L0, \
  int n = N, \
  int m1 = N, \
  int m2 = N, \
  const char * mode = "a"

#include "profiles_scalar.h"
#include "profiles_product.h"
#include "profiles_dissipation.h"
#include "profiles_dissipation_scalar.h"

#undef SETUP_PROFILE_PLANE
#undef PROFILE_PARAMS
#undef PROFILE_OPEN
#undef PROFILE_CLOSE
#undef _mindelta