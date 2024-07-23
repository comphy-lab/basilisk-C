/**
# Rayleigh collapse of a compressible gas bubble

A gas bubble at lower pressure collapses in a viscous fluid in the
absence of surface tension forces.

The fluids are either viscous or inviscid and we compare the evolution
of the radius with the [Rayleigh-Plesset](#brennen2014) and
[Keller-Miksis](#keller1980) models.

See also section 4.2.4 of [Fuster & Popinet,
2018](/src/compressible/two-phase.h#fuster2018).

## Results for an inviscid fluid

~~~gnuplot Bubble radius as a function of time
set xlabel 't/tR'
set ylabel 'R/R_0'
set key bottom
pg0 = 100.
pinf = 5.*pg0
tR = 0.915*sqrt(1./(pinf - pg0))
plot "../collapse-inviscid/log" u ($1/tR):($2*3.)**(1./3.) every 4 w p pt 6 t "Basilisk", \
     "RPinviscid.dat" u ($1/tR):2 w l lw 2 t 'Rayleigh-Plesset',			    \
     "" u ($1/tR):3 w l lw 2 t 'Keller-Miksis'
~~~

For adiabatic gas transformations $P_b V_b^{\gamma}$ should remain
constant inside the bubble.

~~~gnuplot Entropy errors
set ylabel 'p V^{/Symbol g}'
plot "../collapse-inviscid/log" u ($1/tR):(($2*3.)**(1.4)*$4) not w l 
~~~

The bubble does not remain spherical.

~~~gnuplot Interfaces
set xlabel 'x'
set ylabel 'y'
set size ratio -1
plot "../collapse-inviscid/out" u 1:2 not w l 
~~~

## Results for a viscous fluid

~~~gnuplot Bubble radius as a function of time
reset
set xlabel 't/tR'
set ylabel 'R/R_0'
set key bottom
plot "log" u ($1/tR):($2*3.)**(1./3.) every 3 w p pt 6 t "Basilisk", \
     "../collapse-spherical/log" u ($1/tR):($2*3.)**(1./3.) every 3 w p pt 8 t "Basilisk (1D)", \
     "RP.dat" u ($1/tR):2 w l lw 2 t 'Rayleigh-Plesset',             \
     "" u ($1/tR):3 w l lw 2 t 'Keller-Miksis',                      \
     "RPinviscid.dat" u ($1/tR):3 w l lw 2 t "Keller-Miksis (inviscid)"
~~~

For adiabatic gas transformations $P_b V_b^{\gamma}$ should remain
constant inside the bubble.

~~~gnuplot Entropy errors
set ylabel 'p V^{/Symbol g}'
plot "log" u ($1/tR):(($2*3.)**(1.4)*$4) not w l 
~~~

The bubble remains more spherical than in the inviscid case.

~~~gnuplot Interfaces
set xlabel 'x'
set ylabel 'y'
set size ratio -1
plot "out" u 1:2 not w l 
~~~

## References

~~~bib
@book{brennen2014,
  title={Cavitation and bubble dynamics},
  author={Brennen, Christopher E},
  year={2014},
  publisher={Cambridge university press}
}

@article{keller1980,
  title={Bubble oscillations of large amplitude},
  author={Keller, Joseph B and Miksis, Michael},
  journal={The Journal of the Acoustical Society of America},
  volume={68},
  number={2},
  pages={628--633},
  year={1980},
  publisher={Acoustical Society of America}
}
~~~

We run both the (1D) spherically-symmetric and the (2D) axisymmetric
versions. */

#if SPHERICAL
# include "spherisym.h"
#else
# include "axi.h"
#endif
#include "bubble.h"

int main()
{
  pinf = 5.*pg0;
  tend = 2.50001*0.915*sqrt (1./(pinf - pg0)); // fixme: stops too early with 2.5

#if INVISCID
  MAXLEVEL = 12;
#else // !INVISCID
  MAXLEVEL = 11;
  double Reynolds = 10.;
  mu1 = sqrt(pinf - pg0)/Reynolds;

  /** 
  Rayleigh-Plesset like models neglect the gas viscosity. Here we put
  a small value. */
  
  mu2 = mu1*0.0001;
#endif // !INVISCID
  
  CFLac = 10.;
  gamma2 = 1.4;
  gamma1 = 7.14;
  PI1 = 300*pg0;
     
  L0 = 50;
#if TREE  
  N = 1 << MINLEVEL;
#else
  N = 1 << MAXLEVEL;
#endif
  run();
}
