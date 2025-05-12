#define BGHOSTS 2

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>

#define INTERFACIAL (interfacial_point[] > 0.1)
#define NOT_INTERFACIAL (interfacial_point[] < 0.1)

// If not defined by the user we fixe a default value for F_ERR, the accepted error
// over f to avoid division by zero. */
#define F_ERR 1.0e-10
#define F_COMPARE(a,b,tol) (fabs(a-b)>tol)

#if dimension == 2
#define CENTER_CELL (_i == point.i && _j == point.j)
#define NOT_CENTER_CELL (!CENTER_CELL)
#elif dimension == 3
#define CENTER_CELL (_i == point.i && _j == point.j && _k == point.k)
#define NOT_CENTER_CELL (!CENTER_CELL)
#endif

#define NOT_BOUNDARY (!is_boundary(cell))

#define is_LOCAL_and_ACTIVE (is_local(cell) && is_active(cell))
#define NOT_LOCAL_and_ACTIVE (!(is_local(cell) && is_active(cell)))

#define is_LEAF (is_leaf(cell))
#define NOT_LEAF (!is_LEAF)

#define is_VAPOR (f[] < 0.5)
#define is_LIQUID (f[] > 0.5)

// for debuging
#if _MPI
#define STOP(i) MPI_Abort(MPI_COMM_WORLD, i)
#else
#define STOP(i) fprintf(stderr, "stop error\n");fflush(stdout); exit(i)
#endif

// constants for the thermal conduction
double DIFFUSION_TOL = 1.0e-3;
#define TIME_FACTOR
double time_factor = 1.0e3;

// increase projection stability - not required
#define FILTERED

// default refinement levels
int MIN_LEVEL = 6;
int LEVEL = 6;
int MAX_LEVEL = 6;

// number of refined cells near interface using LS refinement 
#define REFINED_CELLS 5.0

// useful for refining in the init event
#define DX_MIN (L0 / (1 << MAX_LEVEL))
#define DX_START (L0 / (1 << LEVEL))
#define dR_refine (max(REFINED_CELLS*DX_MIN, 2.0*DX_START))
//
#define DT_MAX (T_END / (200.0+1.0e-4))
#define DELTA_T (T_END / (200.0+1.0e-4))

// simplification
#ifdef RHOCP1
#define rhocp_L (1.0)
#define rhocp_V (1.0)
#else
#define rhocp_L (rho_L*cp_L)
#define rhocp_V (rho_V*cp_V)
#endif

#include "utils.h"

#include "tag.h"
#if AXI
#include "axi.h"
#endif // AXI

#include "centered_evap.h"

#include "two-phase-evap.h"

#include "reduced.h"

#include "navier-stokes/conserving.h"

#include "tension.h"

// to determine the global minimum cell dimension
double get_delta_min() {
#if TREE
  int max_level = MAX_LEVEL;
#else
  int max_level = LEVEL;
#endif  // TREE
  return L0 / (1 << max_level);
}

// find minlevel - mainly use for post-processing
double get_min_level(){
  int min_level = 1e9;
  foreach(reduction(min:min_level)){
    min_level = min(point.level, min_level);
  }
  return min_level;
}

// find maxlevel - mainly use for post-processing
double get_max_level(){
  int max_level = -1e9;
  foreach(reduction (max:max_level)){
    max_level = max(point.level, max_level);
  }
  return max_level;
}

// check for nan volume fraction
void f_nan_check(int call) {
  int QUIT = 0;
  foreach(reduction(max:QUIT)){
    static bool once = true;
    if (once) {
      if (isnan(f[]) || isinf(f[])) {
        printf("f[] = %g call %d\n", f[], call);
        once = false;
        QUIT = 1;
      }
    }
  }
  if (QUIT) {
    STOP(2);
  }
}

// LS_dist scalar for adapt refinement near the interface
scalar LS_dist[];

#define SCALE_LS (1.0/get_delta_min())

// init volume fraction and LS_dist at the same time
#define fraction_LS(f, func)  \
    {                         \
        fraction(f, func);    \
        foreach ()            \
            LS_dist[] = SCALE_LS*func; \
    }

#include "temperature-phase-change.h"

#ifdef PERFS
#include "navier-stokes/perfs.h"
#endif  // PERFS

// for post-processing and debugging MPI cell distribution
#if _MPI
scalar PID_info[];
#endif

scalar T_L[], T_V[];
scalar fE_L[], fE_V[];

void setup_evap() {
  TOLERANCE = 1.0e-5;
  
  f.sigma = sigma_lv;

  rho1 = rho_L;
  rho2 = rho_V;

  mu1 = mu_L;
  mu2 = mu_V;

  DT = DT_MAX;

  evap_const = (1.0 / rho_V - 1.0 / rho_L);

  fE_L.inverse = false;
  fE_V.inverse = true;

  T_L.inverse = false;
  T_V.inverse = true;
}

