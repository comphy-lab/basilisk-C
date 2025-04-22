/**
# Dispersion relations for various models

This test case measures the dependence of the oscillation frequency of
a gravity wave of wavenumber $k$ on the water depth $h$. The
corresponding phase velocity $c=\omega/k$ is then compared to the
exact linear dispersion relation
$$
c_e = \sqrt{gk\tanh(kh)}
$$

More details can be found in [Popinet (2020)](/Bibliography#popinet2020).

We use a 1D grid and either the Green-Naghdi solver or the multi-layer
solver. */

#include "grid/multigrid1D.h"
#if !ML
  #include "green-naghdi.h"
#else
  #include "layered/hydro.h"
  #include "layered/nh.h"
  // necessary for stability at high h0 of the box scheme
  #include "layered/remap.h"
#endif

double emax = 1.;
double h0 = 0.1;

/**
We initialise a wave of small relative amplitude ($10^{-3}$) to be in
the linear regime. For the multilayer model, we also run a case with
three "optimised" non-uniform layer thicknesses. */

event init (i = 0)
{
#if ML
  if (nl == 3 && emax == 1.5)
    beta[0] = 0.68, beta[1] = 0.265, beta[2] = 0.055;
  foreach()
    foreach_layer()
      h[] = beta[point.l]*h0*(1. + 0.001*cos(x));
#else
  foreach()
    h[] = h0*(1. + 0.001*cos(x));
#endif
}

/**
The code below locates the time of zero "up-crossing" of the amplitude
in the middle of the wave and computes the corresponding averaged
wave period. */

double Tm = 0., Es = 1., Ee = 1.;
int nm = 0;

event logfile (i++)
{
  double pe = 0., ke = 0.;
#if ML  
  foreach() {
    double H = 0.;
    foreach_layer() {
      ke += Delta*h[]*(sq(u.x[]) + sq(w[]))/2.;
      H += h[];
    }
    pe += Delta*G*sq(H - h0)/2.;
  }
#else
  foreach() {
    int l = 0;
    for (vector u in ul)
      ke += Delta*layer[l++]*h[]*sq(u.x[])/2.;
    pe += Delta*G*sq(h[] - h0)/2.;
  }
#endif

#if ML
  double H = 0.;
  foreach_point (L0/2.)
    foreach_layer()
      H += h[];
#else
  double H = interpolate (h, L0/2., linear = false);
#endif
  double dh = (H - h0)/h0;
#if 0
  static FILE * fp = fopen ("timeseries", "w");
  char name[80];
  sprintf (name, "w%d", nl - 1);
  scalar w = lookup_field (name);
  fprintf (fp, "%g %g %g %g %g\n", t/(2.*pi/sqrt(tanh(h0))),
	   H - h0, dh, ke, pe);
#endif
  
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

#if ML
event profile (t = 9.25*2.*pi/sqrt(tanh(h0)))
{
  if (nl == 5) {
    char name[80];
    sprintf (name, "profile-%g", h0);
    FILE * fp = fopen (name, "w");
    foreach_point (L0/2.) {
      double zl = zb[];
      foreach_layer() {
	fprintf (fp, "%g %g %g\n", zl + h[]/2., w[], phi[]);
	zl += h[];
      }
    }
    fclose (fp);
  }
}
#endif // ML

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
	   "set title 'nl = %d, h0 = %g, t = %.2f'\n"
	   "unset key\n"
	   "p []'-' u 1:4 w l, '' u 1:6 w l, '' u 1:8 w l, '' u 1:10 w l\n",
	   // "p []'-' u 1:2 w l\n",
	   nl, h0, t);
  foreach() {
    double etap = zb[];
    foreach_layer()
      etap += h[];
    fprintf (fp, "%g %g", x, eta[] - etap);
    //    eta[] = etap;
    foreach_layer()
      fprintf (fp, " %g %g", h[], u.x[]);
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
	   h0, sqrt(tanh(h0)), 2.*pi/(Tm/nm), (Ee - Es)/Es,
#if !ML
	   mgD.i, mgD.nrelax
#else	   
	   mgp.i, mgp.nrelax
#endif
	   , i);
  fflush (stderr);
  emax = 2.*pi/(Tm/nm)/sqrt(tanh(h0));
}

/**
We run for several numbers of layers and many water depths, stopping
when the error becomes too large. */

int main()
{
  periodic (right);
  size (2.*pi);
  DT = HUGE [0]; // dimensionless time
  N = 128;
#if ML
  CFL_H = 0.5;
  TOLERANCE = 1e-6;
  linearised = true;
#if 1
  for (nl = 1; nl <= 5; nl++)
#else
  nl = 1;
#endif
  {
    char name[80];
    sprintf (name, "log-%d", nl);
    freopen (name, "w", stderr);
    emax = 1.;
    DT = 2.*pi/sqrt(tanh(h0))/200.;
#if 1
    for (h0 = 0.1; h0 <= 100. && emax > 0.96; h0 *= 1.3)
#else
    h0 = 0.37;
#endif
    {
      run();
    }
  }
#if 1
  freopen ("log-opt", "w", stderr);
  nl = 3;
  emax = 1.;
  for (h0 = 0.1; h0 <= 100. && emax > 0.96; h0 *= 1.3) {
    emax = 1.5;
    run();
  }
#endif
#else
  alpha_d = 1.159;
  for (h0 = 0.1; h0 <= 100. && emax > 0.96 && emax < 1.04; h0 *= 1.2)
    run();
#endif
}

/**
Note that the standard Green-Naghdi model has a [2,2] Padé approximant
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
      omega3_keller(0.68*x,0.265*x,0.055*x)/sqrt(tanh(x)) \
          t 'Optimised 3 layers' lc 0, \
      'log-opt' u ($1):($3/sqrt(tanh($1))) t '' pt 14 lc 0, \
      omega_gn(1,x)/sqrt(tanh(x)) t 'Green-Naghdi (1 layer)' lc 4, \
      '../dispersion-gn/log' u ($1):($3/sqrt(tanh($1))) pt 15 t ''
~~~

~~~gnuplot Relative energy evolution
set ylabel 'Energy variation per period (%)'
set yr [*:*]
set key top left
plot  'log-1' u ($1):($4*10.) t '1 layer' pt 5 lt 2, \
      'log-2' u ($1):($4*10.) t '2 layers' pt 6 lt 4, \
      'log-3' u ($1):($4*10.) t '3 layers' pt 9 lt 6, \
      'log-4' u ($1):($4*10.) t '4 layers' pt 10 lt 8, \
      'log-5' u ($1):($4*10.) t '5 layers' pt 13 lt 10, \
      'log-opt' u ($1):($4*10.) t 'Optimised 3 layers' pt 14 lc 0, \
      '../dispersion-gn/log' u ($1):($4*10.) pt 15 lt 14 t 'Optimised Green-Naghdi'
~~~

~~~gnuplot Vertical profiles of vertical velocity at $t=37/4T$
reset
g = 1
k = 1
A = 1e-3
omega(H) = sqrt(g*H*k**2*tanh(k*H)/(k*H))
wmax(H) = A*H*g*k/omega(H)*sinh(k*H)/cosh(k*H)
ws(x,z,t,H) = - A*H*g*k/omega(H)*sinh(k*z)/cosh(k*H)*sin(omega(H)*t - k*x)
set xlabel 'z/H'
set ylabel 'w/w_{max}'
set key top left
range = "1.06045 2.32981 5.11859 24.7065"
plot [0:1][:1]							   \
  for [H in range]						   \
    ws(pi,x*H,9.25*2.*pi/omega(H),H)/wmax(H) w l t 'kH = '.H,	   \
  for [H in range]						   \
    'profile-'.H u ($1/H):($2/wmax(H)) w p t ''
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
