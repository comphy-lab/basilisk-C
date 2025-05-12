
/**
# LEVEL-SET FUNCTION ADAPT
*/

#include "../LS_funcs/LS_reinit.h"

void reinit_LS(scalar NEAR) {
  int it_MAX = 3;

  static bool FIRST = true;
  if (FIRST){
    double LS_max = -1.0e9;
    double LS_min = 1.0e9;
    foreach(reduction(max:LS_max) reduction (min:LS_min)){
      LS_max = max(LS_dist[], LS_max);
      LS_min = min(LS_dist[], LS_min);
    }
    if (LS_max == LS_min){
      printf("Error: LS_dist has not been initialized!\n");
      printf("Use fraction_LS marco instead of fraction in the init event\n");
      exit(1);
    }
    FIRST = false;
  }
  double delta_min = get_delta_min();
  foreach () {
    if (NEAR[] > 0.5) {
      LS_dist[] = delta_min * SCALE_LS * (f[] - 0.5);
    }
  }
  LS_reinit(LS_dist, it_max = it_MAX);
}

void clamp_f() {
  foreach () { f[] = (f[] < F_ERR ? 0. : f[] > 1. - F_ERR ? 1. : f[]); }
  return;
}

static double rand_double() { return rand() / (double)RAND_MAX; }

#if TREE

#define NEAR_CHECK(a) (a[-1] > 0.1 || a[1] > 0.1 || a[0] > 0.1)

#if dimension == 2
#define Uemax uemax, uemax
#elif dimension == 3
#define Uemax uemax, uemax, uemax
#endif

event adapt(i++) {

  scalar near_interface[];
  scalar impose_refine[];
  clamp_f();

  foreach () {
    if (interfacial(point, f)) {
      near_interface[] = 1.0;
    } else {
      near_interface[] = 0.0;
    }
  }

  reinit_LS(near_interface);

  double delta_min = get_delta_min();

  // exponential decay
  double cells = REFINED_CELLS;
  double decay = 10.0;
  foreach () {
    double dist = max(fabs(LS_dist[]) / (delta_min*SCALE_LS), cells) - cells;
    impose_refine[] = rand_double() * exp(-decay * dist);
    // made this for the film boiling case
#ifdef REFINE_VAPOR_REGION
    if (f[]<0.5){
      impose_refine[] = rand_double();
    }
#endif
  }
  impose_refine.dirty = true;
  boundary({impose_refine});
  double refine_tol = 1.0e-6;

  adapt_wavelet({impose_refine, f, T_L, T_V, u},
              (double[]){refine_tol, femax, Temax, Temax, Uemax},
              minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif  // TREE

