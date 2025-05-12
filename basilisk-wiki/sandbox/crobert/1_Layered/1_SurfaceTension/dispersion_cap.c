/**  
# Dispersion relations of pure capillary waves with PhiS implementation

This test case measures the dependence of the oscillation frequency of
a capillary wave of wavenumber $k$ on the water depth $h$. The
corresponding phase velocity $c=\omega/k$ is then compared to the
exact linear dispersion relation
$$
c_e = \sqrt{\frac{\gamma}{\rho}k\tanh(kh)}
$$ 

We use a 1D grid and the multi-layer solver. */

#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "../nh_PhiS.h"
// necessary for stability at high h0 of the box scheme
#include "layered/remap.h"
#include "../tension.h"

double emax = 1.;
double h0 = 0.1;

/**
We initialise a wave of small relative amplitude ($10^{-3}$) to be in
the linear regime. Surface tension is initialized to 1. */

event init (i = 0)
{
  foreach(){
    sigma[] = 1.;
    foreach_layer()
      h[] = beta[point.l]*h0*(1. + 0.001*cos(x));
  }
}

/**
The code below locates the time of zero "up-crossing" of the amplitude
in the middle of the wave and computes the corresponding averaged
wave period. */

double Tm = 0., Es = 1., Ee = 1.;
int nm = 0;

event logfile (i += 20) {
  double pe = 0., ke = 0.;
  
  foreach() {
    foreach_layer()
      ke += h[]*Delta*(sq(u.x[]) + sq(w[]))/2.; //Add the h[] : slight correction for error calculation
    pe += sigma[]*sq(eta[1]-eta[-1])/(8.*Delta);
  }
  
  Point point = locate (L0/2.);
  double dh = (eta[] - h0)/h0;

  static double told = 0., hold = 0., eold = 0., tsold = -1;
  if (i == 0) {
    Tm = 0., nm = 0;
    Es = 0.;
    told = 0., hold = 0., tsold = -1;
  }
  else {
    if (i > 1 && hold < 0. && dh > 0.) {
      // this is an (upward) zero-crossing at time ts
      double ts = told - hold*(t - told)/(dh - hold);
      Ee = eold - hold*(ke + pe - eold)/(dh - hold);
      if (Es == 0.)
	Es = Ee;
      //      fprintf (fp, "ts %g 0 %g\n", ts, Ee);
      if (tsold > 0. && nm < 4) {
	// we average the periods
	Tm += ts - tsold;
	nm++;
      }
      tsold = ts;
    }
    told = t;
    hold = dh;
    eold = ke + pe;
  }
}

/**
After 9.25 (exact) wave periods, we dump the vertical velocity
profiles. */

event profile (t = 9.25*2.*pi/sqrt(tanh(h0)))
{
  if (nl == 5) {
    char name[80];
    sprintf (name, "profile-%g", h0);
    FILE * fp = fopen (name, "w");
    Point point = locate (L0/2.);
    double zl = zb[];
    foreach_layer() {
      fprintf (fp, "%g %g %g\n", zl + h[]/2., w[], phi[]);
      zl += h[];
    }
    fclose (fp);
  }
}

/**
After ten (exact) wave periods, we stop and dump the solution. */

event dumps (t = 10.*2.*pi/sqrt(tanh(h0))) // ten wave periods
{
  FILE * fp = fopen ("prof", "w");
  fprintf (fp, "x");
  for (scalar s in all)
    fprintf (fp, " %s", s.name);
  fprintf (fp, "\n");
  foreach() {
    fprintf (fp, "%g", x);
    for (scalar s in all)
      fprintf (fp, " %g", s[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "\n");
  fclose (fp);
}

#if 0
event movie (i += 1) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fputs ("set term x11\n", fp);
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "unset key\n"
	   "p []'-' u 1:3 w l, '' u 1:5 w l, '' u 1:7 w l, '' u 1:9 w l\n", t);
  foreach() {
    fprintf (fp, "%g", x);
    foreach_layer()
      fprintf (fp, " %g", h[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);

  fprintf (stderr, "%g %g\n", t, statsf(eta).max); 
}
#endif

/**
We output the wave period and relative energy variation and compute the
relative error `emax`. */

event end (t = end)
{
  fprintf (stderr, "%g %g %g %g %d %d %d\n",
	   h0, sqrt(tanh(h0)), 2.*pi/(Tm/nm), (Ee - Es)/Es, mgp.i, mgp.nrelax, i);
  fflush (stderr);
  emax = 2.*pi/(Tm/nm)/sqrt(tanh(h0));
}

/**
We run for several numbers of layers and many water depths, stopping
when the error becomes too large. */

int main()
{
  nh_tension = true;
  periodic (right);
  size (2.*pi);
  N = 64;
  TOLERANCE = 1e-6;
  G = 0.;
  for (nl = 1; nl <= 5; nl++) {
    char name[80];
    sprintf (name, "log-%d", nl);
    freopen (name, "w", stderr);
    emax = 1.;
    DT = 2.*pi/sqrt(tanh(h0))/100.;
    for (h0 = 0.1; h0 <= 100. && emax > 0.96; h0 *= 1.3)
      run();
  }
}

/**
Note that the standard Green-Naghdi model has a [2,2] Pade approximant
dispersion relation (see e.g. [Clamond et al. (2017)](#clamond2017)). The
same as [Nwogu (1993)](#nwogu1993). Here we used the optimized version
of [Chazel et al. (2011)](#chazel2011).

The discrete dispersion relations can be computed using [this Maxima
script](dispersion.mac).

~~~gnuplot Dispersion relation
g = 1.
k = 1.

omega1_keller(h_0)=sqrt((4*h_0*g*k**2)/(h_0**2*k**2+4))

omega2_keller(h_0,h_1)=sqrt(((4*h_0*h_1**2+4*h_0**2*h_1)*g*k**4+(16*h_1+16*h_0)*g*k**2)/(h_0**2*h_1**2*k**4+(4*h_1**2+16*h_0*h_1+4*h_0**2)*k**2+16))

omega3_keller(h_0,h_1,h_2)=sqrt((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*g*k**6+((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*g*k**4+(64*h_2+64*h_1+64*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*k**6+((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*k**4+(16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*k**2+64))

omega4_keller(h_0,h_1,h_2,h_3)=sqrt(((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*h_3**2+4*h_0**2*h_1**2*h_2**2*h_3)*g*k**8+(((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*h_3**2+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*g*k**6+((64*h_2+64*h_1+64*h_0)*h_3**2+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3+(64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*g*k**4+(256*h_3+256*h_2+256*h_1+256*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*h_3**2*k**8+(((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*h_3**2+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3+4*h_0**2*h_1**2*h_2**2)*k**6+((16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*h_3**2+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*k**4+(64*h_3**2+(256*h_2+256*h_1+256*h_0)*h_3+64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*k**2+256))

omega5_keller(h_0,h_1,h_2,h_3,h_4)=sqrt((((((4*h_0*h_1**2+4*h_0**2*h_1)*h_2**2+4*h_0**2*h_1**2*h_2)*h_3**2+4*h_0**2*h_1**2*h_2**2*h_3)*h_4**2+4*h_0**2*h_1**2*h_2**2*h_3**2*h_4)*g*k**10+((((16*h_1+16*h_0)*h_2**2+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2+16*h_0*h_1**2+16*h_0**2*h_1)*h_3**2+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_4**2+(((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3**2+((64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_3+16*h_0**2*h_1**2*h_2**2)*h_4+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3**2+16*h_0**2*h_1**2*h_2**2*h_3)*g*k**8+(((64*h_2+64*h_1+64*h_0)*h_3**2+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3+(64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_4**2+((64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3**2+((256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_3+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_4+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3**2+((64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_3+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*g*k**6+((256*h_3+256*h_2+256*h_1+256*h_0)*h_4**2+(256*h_3**2+(1024*h_2+1024*h_1+1024*h_0)*h_3+256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_4+(256*h_2+256*h_1+256*h_0)*h_3**2+(256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_3+(256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*g*k**4+(1024*h_4+1024*h_3+1024*h_2+1024*h_1+1024*h_0)*g*k**2)/(h_0**2*h_1**2*h_2**2*h_3**2*h_4**2*k**10+((((4*h_1**2+16*h_0*h_1+4*h_0**2)*h_2**2+(16*h_0*h_1**2+16*h_0**2*h_1)*h_2+4*h_0**2*h_1**2)*h_3**2+((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3+4*h_0**2*h_1**2*h_2**2)*h_4**2+(((16*h_0*h_1**2+16*h_0**2*h_1)*h_2**2+16*h_0**2*h_1**2*h_2)*h_3**2+16*h_0**2*h_1**2*h_2**2*h_3)*h_4+4*h_0**2*h_1**2*h_2**2*h_3**2)*k**8+(((16*h_2**2+(64*h_1+64*h_0)*h_2+16*h_1**2+64*h_0*h_1+16*h_0**2)*h_3**2+((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3+(16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_4**2+(((64*h_1+64*h_0)*h_2**2+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2+64*h_0*h_1**2+64*h_0**2*h_1)*h_3**2+((64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*h_3+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_4+((16*h_1**2+64*h_0*h_1+16*h_0**2)*h_2**2+(64*h_0*h_1**2+64*h_0**2*h_1)*h_2+16*h_0**2*h_1**2)*h_3**2+((64*h_0*h_1**2+64*h_0**2*h_1)*h_2**2+64*h_0**2*h_1**2*h_2)*h_3+16*h_0**2*h_1**2*h_2**2)*k**6+((64*h_3**2+(256*h_2+256*h_1+256*h_0)*h_3+64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_4**2+((256*h_2+256*h_1+256*h_0)*h_3**2+(256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_3+(256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_4+(64*h_2**2+(256*h_1+256*h_0)*h_2+64*h_1**2+256*h_0*h_1+64*h_0**2)*h_3**2+((256*h_1+256*h_0)*h_2**2+(256*h_1**2+1024*h_0*h_1+256*h_0**2)*h_2+256*h_0*h_1**2+256*h_0**2*h_1)*h_3+(64*h_1**2+256*h_0*h_1+64*h_0**2)*h_2**2+(256*h_0*h_1**2+256*h_0**2*h_1)*h_2+64*h_0**2*h_1**2)*k**4+(256*h_4**2+(1024*h_3+1024*h_2+1024*h_1+1024*h_0)*h_4+256*h_3**2+(1024*h_2+1024*h_1+1024*h_0)*h_3+256*h_2**2+(1024*h_1+1024*h_0)*h_2+256*h_1**2+1024*h_0*h_1+256*h_0**2)*k**2+1024))

alpha=1.159
omega_gn(k,h)=sqrt(g*h*k**2*(1.+(alpha-1.)*(k*h)**2/3.)/(1.+alpha*(k*h)**2/3.))

set xlabel 'kH'
set ylabel 'c/c_e'
set key bottom left
set xr [0.1:100]
set yr [0.96:1.04]
set logscale x
set grid
plot  omega1_keller(x)/sqrt(tanh(x)) t '1 layer', \
      'log-1' u ($1):($3/sqrt(tanh($1))) t '' pt 5, \
      omega2_keller(x/2.,x/2.)/sqrt(tanh(x)) t '2 layers', \
      'log-2' u ($1):($3/sqrt(tanh($1))) t '' pt 6, \
      omega3_keller(x/3.,x/3.,x/3.)/sqrt(tanh(x)) t '3 layers', \
      'log-3' u ($1):($3/sqrt(tanh($1))) t '' pt 9, \
      omega4_keller(x/4.,x/4.,x/4.,x/4.)/sqrt(tanh(x)) t '4 layers', \
      'log-4' u ($1):($3/sqrt(tanh($1))) t '' pt 10, \
      omega5_keller(x/5.,x/5.,x/5.,x/5.,x/5.)/sqrt(tanh(x)) t '5 layers', \
      'log-5' u ($1):($3/sqrt(tanh($1))) t '' pt 13, \
      omega_gn(1,x)/sqrt(tanh(x)) t 'Green-Naghdi (1 layer)' lc 4, \
      '../../../../../../src/test/dispersion-gn.ref' u ($1):($3/sqrt(tanh($1))) pt 15 t ''
~~~

~~~gnuplot Relative energy evolution
set ylabel 'Energy variation per period (%)'
set yr [*:*]
set key bottom left
plot  'log-1' u ($1):($4*10.) t '1 layers' pt 6 lt 2, \
      'log-2' u ($1):($4*10.) t '2 layers' pt 6 lt 4, \
      'log-3' u ($1):($4*10.) t '3 layers' pt 9 lt 6, \
      'log-4' u ($1):($4*10.) t '4 layers' pt 10 lt 8, \
      'log-5' u ($1):($4*10.) t '5 layers' pt 13 lt 10, \
      '../../../../../../src/test/dispersion-gn.ref' u ($1):($4*10.) pt 15 lt 14 t 'Optimised Green-Naghdi'
~~~

## References

~~~bib
@article{nwogu1993,
  title={Alternative form of Boussinesq equations for nearshore 
         wave propagation},
  author={Nwogu, Okey},
  journal={Journal of waterway, port, coastal, and ocean engineering},
  volume={119},
  number={6},
  pages={618--638},
  year={1993},
  publisher={American Society of Civil Engineers}
}

@article{clamond2017,
  title={Conservative modified Serre--Green--Naghdi equations with 
         improved dispersion characteristics},
  author={Clamond, Didier and Dutykh, Denys and Mitsotakis, Dimitrios},
  journal={Communications in Nonlinear Science and Numerical Simulation},
  volume={45},
  pages={245--257},
  year={2017},
  publisher={Elsevier}
}

@article{chazel2011,
  title={Numerical simulation of strongly nonlinear and dispersive waves 
         using a Green--Naghdi model},
  author={Chazel, Florent and Lannes, David and Marche, Fabien},
  journal={Journal of Scientific Computing},
  volume={48},
  number={1-3},
  pages={105--116},
  year={2011},
  publisher={Springer}
}
~~~
*/