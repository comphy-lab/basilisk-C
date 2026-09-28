#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/diffusion.h"
#if !NO_REMAP && !NOREP
# include "layered/remap.h"
#endif

scalar T[];
double H0;

#if T_CST_LIKE
double Ts = 20.;
#else
double Ts = 1.;
#endif

int main(int argc, char *argv[])
{
  L0 = 1.;
  N = 1;      
  nl = 30;
  H0 = 100.;
  G = 9.81;
  periodic (top);
  periodic (left);

  run();

}

event init (i = 0) {
  T = new scalar[nl];
  /**
  Enrol T in the tracers list, exactly like layered/dr.h does in the T_cst
  test: from then on hydro's advection event and (if included) remap.h's
  vertical_remapping() act on T at every step. */
  tracers = list_append (tracers, T);

  foreach(){
    zb[] = -H0;
    eta[] = 0.;
    foreach_layer() {
      h[] = H0/nl;       
      T[] = Ts;       
      foreach_dimension()
        u.x[] = 0.;
      w[] = 0.;
    }
  }

#if NOREP
  /**
  ## 100 repeated calls of the stock vertical_diffusion()
  on the constant field, inside init, then exit. No advection, no remap,
  no pressure: only the tridiagonal solve runs. Exact answer: T stays Ts */
  for (int it = 0; it < 100; it++) {
    foreach() {
      vertical_diffusion (point, // point
                          h,     // h
                          T,     // scalar
                          1.,    // dt
                          1.,    // D
                          0.,    // dst
                          0.,    // sb
                          HUGE); // lambdab
    }
    foreach(cpu){
      foreach_layer()
        fprintf(stderr, "%d %d %.12e\n", it, point.l, T[]-1.);
    }
    fprintf(stderr,"\n");
  }
  exit(0);
#endif
}

event viscous_term (i++){
  foreach() {
    vertical_diffusion (point, // point
                        h,     // h
                        T,     // scalar
                        1.,    // dt
                        1.,    // D
                        0.,    // dst
                        0.,    // sb
                        HUGE); // lambdab
  }

  foreach(cpu){
    foreach_layer()
      fprintf(stderr, "%d %d %.12e\n", i, point.l, T[]-Ts);
  }
  fprintf(stderr,"\n");
}

event stop (i=200){
  delete ({T});
  exit(0);
}



