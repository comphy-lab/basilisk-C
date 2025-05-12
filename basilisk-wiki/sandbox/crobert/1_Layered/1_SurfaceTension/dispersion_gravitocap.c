/**
# Dispersion relations of gravito-capillary waves 

This test case measures the dependence of the oscillation frequency of
a gravito-capillary wave of wavenumber $k$ on the water depth $h$. The
corresponding phase velocity $c=\omega/k$ is then compared to the
exact linear dispersion relation
$$
c_e = \sqrt{\left(\frac{g}{k} + \gamma{\rho}k\right)\tanh(kh)}
$$ 

We use a 1D grid, the multi-layer solver and the surface tension additive. */

#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "../nh_PhiS.h"

#include "layered/remap.h"
#include "../tension.h"

double emax = 1.;
double h0 = 10.;
double k = 1.;
double L = 1;
double c0 = 1;
#define sig 1.

/**
We initialise a wave of small relative amplitude ($10^{-3}$) to be in
the linear regime. Surface tension is initialized to sig. 
Averaged depth h0 is 10.*/

event init (i = 0)
{
  foreach(){
    sigma[] = sig;
    foreach_layer()
      h[] =  beta[point.l]*h0*(1. + 0.001*cos(k*x));
  }
}

/**
The code below locates the time of zero "up-crossing" of the amplitude
in the middle of the wave and computes the corresponding averaged
wave period. */

double Tm = 0., Es = 1., Ee = 1.;
int nm = 0;

event logfile (i++) {
  double pe = 0., ke = 0.;
  
  foreach() {
    double H = 0;
    foreach_layer(){
      ke += h[]*Delta*(sq(u.x[]) + sq(w[]))/2.;
      H += h[];
    }
    pe += sigma[]*sq(eta[1]-eta[-1])/(8.*Delta);
    pe += Delta*G*sq(H - h0)/2.;
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
After ten (exact) wave periods, we stop and dump the solution. */

event dumps (t = 10.*L/c0) // ten wave periods
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

/**
We output the wave period and relative energy variation and compute the
relative error `emax`. */

event end (t = end)
{
  emax = (L/c0)/(Tm/nm);
  fprintf (stderr, "%g %g %g %g %g %g\n",
	   k, h0, c0, L/(Tm/nm), emax, (Ee - Es)/Es);
  fflush (stderr);
}

/**
We run for several numbers of layers and many water depths, stopping
when the error becomes too large. */

int main()
{
  nh_tension = true;
  periodic (right);
  N = 128;
  TOLERANCE = 1e-6;
  G = 1.;
  for (nl = 1; nl <= 8; nl*=2){
    char name[80];
    sprintf (name, "log-%d", nl);
    freopen (name, "w", stderr);
    emax = 1.;
    for (k = 0.01; k <= 10. && emax > 0.9; k *= 1.3) {
      c0 = sqrt(tanh(k*h0)*(sig*k + G/k));
      L = 2.*pi/k;
      size (L);
      run();
    }
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

set xlabel 'k'
set ylabel 'c/c_e'
set key bottom left
set xr [0.01:10]
set yr [0.96:1.04]
set logscale x
set grid
plot  omega1_keller(10*x)/sqrt(tanh(10*x)) t '1 layer', \
      'log-1' u 1:5 t '1 layers' pt 1 lt 2, \
      omega2_keller(10*x/2.,10*x/2.)/sqrt(tanh(10*x)) t '2 layers', \
      'log-2' u 1:5 t '2 layers' pt 2 lt 4, \
      omega4_keller(10*x/4.,10*x/4.,10*x/4.,10*x/4.)/sqrt(tanh(10*x)) t '4 layers', \
      'log-4' u 1:5 t '4 layers' pt 3 lt 6, \
      'log-8' u 1:5 t '8 layers' pt 4 lt 8, \
      'log-16' u 1:5 t '16 layers' pt 6 lt 10
~~~

~~~gnuplot Phase velocity
set ylabel 'c'
set xlabel 'k'
set logscale x
set xr [*:*]
set yr [*:*]
plot  'log-1' u 1:4 t '1 layers' pt 1 lt 2, \
      'log-2' u 1:4 t '2 layers' pt 2 lt 4, \
      'log-4' u 1:4 t '4 layers' pt 3 lt 6, \
      'log-8' u 1:4 t '8 layers' pt 4 lt 8, \
      sqrt(tanh(10.*x)*(x + 1./x)) t 'dispersion relation'
~~~

~~~gnuplot Energy loss
set ylabel 'energy loss'
set xlabel 'k'
set xr [*:*]
set yr [-0.1:0.1]
plot  'log-1' u ($1):($6*10.) t '1 layers' pt 1 lt 2, \
      'log-2' u ($1):($6*10.) t '2 layers' pt 2 lt 4, \
      'log-4' u ($1):($6*10.) t '4 layers' pt 3 lt 6, \
      'log-8' u ($1):($6*10.) t '8 layers' pt 4 lt 8
~~~
*/