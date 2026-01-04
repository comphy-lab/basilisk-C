/**
# Small-amplitude oscillations of a compressible gas bubble due to surface tension */

#if SPHERICAL
# include "spherisym.h"
#else
# include "axi.h"
#endif
#include "bubble.h"

int main()
{
  pinf = pg0;
  f.sigma = 1.;
  gamma2 = 1.4;
  gamma1 = 7.14;
  PI1 = 30000*pg0;
  tend = 0.4;
  L0 = 40.;
  
  MINLEVEL = 5;
  for (MAXLEVEL = 9; MAXLEVEL <= 12; MAXLEVEL++) {
#if TREE    
    N = 1 << MINLEVEL;
#else
    N = 1 << MAXLEVEL;
#endif
    run();
    fprintf (stderr, "\n\n");
  }
}

/**
Comparison with the Keller-Miksis solution.

~~~gnuplot Bubble radius as a function of time
set xlabel 't'
set ylabel 'R/R_0'
set key bottom
plot for [i=0:3] "log" index i u 1:($2*3.)**(1./3.) w l t sprintf("%d", 2**(i + 5)), \
                 "RP.dat" u 1:3 w l lw 2 t 'Keller-Miksis'
~~~

For adiabatic gas transformations $P_b V_b^{\gamma}$ should remain
constant inside the bubble.

~~~gnuplot Entropy errors
set ylabel 'p V^{/Symbol g}'
plot for [i=0:3] "log" index i u 1:(($2*3.)**(1.4)*$4) w l t sprintf("%d", 2**(i + 5))
~~~
*/
