/**
## Isotropic turbulence (2d)
*/
#include "navier-stokes/centered.h"

int kVerbose = 1;
#define MU (1.0e-6)
#define RHO (1.0)
#include "stress_omega.h"

#define MIN_LEVEL (4)
int MAX_LEVEL = 7;

face vector muc[];

/**
We set periodic boundary for all directions
*/
int main(int argc, char * argv[]) {
  if (argc > 1)
    MAX_LEVEL = atoi(argv[1]);

  L0 = 1.0;
  N = 512;
  kTurbConstants.k_0 = 1e-4;
  kTurbConstants.omega_0 = 100.0;

  CFL = 0.2;

  foreach_dimension()
    periodic (right);

  run();
}

event init (t = 0){
  mu = muc;
do{
foreach(){
  u.x[] = cos(4*2.*pi*y) + noise()*0.01;
  u.y[] = sin(4*2.*pi*x) + noise()*0.01;
}

  boundary(all);
}
#if TREE
    while (
#if dimension == 2
adapt_wavelet((scalar *){u},(double[]){2e-3,2e-3},MAX_LEVEL,MIN_LEVEL).nf
#elif dimension == 3
adapt_wavelet((scalar *){u},(double[]){5e-3,5e-3,5e-3},MAX_LEVEL,MIN_LEVEL).nf
#endif
);
#else
    while (0);
#endif

}

#if TREE
event adapt (i++) {
#if dimension == 2
  adapt_wavelet((scalar *){u},(double[]){2e-3,2e-3},MAX_LEVEL,MIN_LEVEL);
#elif dimension == 3
  adapt_wavelet((scalar *){u},(double[]){5e-3,5e-3,5e-3},MAX_LEVEL,MIN_LEVEL);
#endif
}
#endif

event properties(i=0){
  foreach_face()
    muc.x[] = fm.x[]*(MU);
  boundary((scalar*){muc});
}

event output (t += 0.05;t<=10) {
  output_ppm (nu_t, min = 1e-6, max = 0.01, file = "nu_t.mp4");
}

/**
Visualization of the results

![Turbulence effective viscosity ranging from 1 - 1000x of molecule viscosity](isotropic/nu_t.mp4)
*/