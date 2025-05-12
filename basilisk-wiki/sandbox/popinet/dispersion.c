/**
# Dispersion relations for various models
*/

#include "grid/multigrid1D.h"
#if SV
#  include "shallow1.h"
#elif GN
#  include "green-naghdi.h"
#else
#  include "shallow1-nh.h"
#endif

double h0 = 0.1;

event init (i = 0)
{
  foreach()
    h[] = h0*(1. + 0.01*cos(x));
}

double Tm = 0.;
int nm = 0;

event logfile (i++) {
  Point point = locate (L0/2.);
  double dh = (h[] - h0)/h0;
  static double told = 0., hold = 0., tsold = -1;
  if (i == 0) {
    Tm = 0., nm = 0;
    told = 0., hold = 0., tsold = -1;
  }
  else {
    if (i > 1 && hold*dh < 0.) {
      // this is a zero-crossing at time ts
      double ts = told - hold*(t - told)/(dh - hold);
      if (tsold > 0. && nm < 4) {
	// we average the periods
	Tm += 2.*(ts - tsold);
	nm++;
      }
      tsold = ts;
    }
    told = t;
    hold = dh;
  }
}

event dumps (t = 50)
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

event end (t = end)
{
  fprintf (stderr, "%g %g %g\n", h0, sqrt(tanh(h0)), 2.*pi/(Tm/nm));
}

int main()
{
  periodic (right);
  size (2.*pi);
#ifdef ALPHA
  alpha = ALPHA;
#endif
  for (h0 = 0.1; h0 <= 2.*pi; h0 += 0.1)
    run();
}

/**
~~~gnuplot Dispersion relation
set key top left
set xlabel 'h/lambda'
set ylabel 'c/sqrt(g h)'
g = 1.
k = 1.
omega(h,alpha) = sqrt(g*h*k**2/(1. + (k*h)**2/alpha))
set key top right
plot 'log' u ($1/(2.*pi)):($3/sqrt($1)) t 'alpha = 1 (numerical)', \
      omega(2.*pi*x,1.)/sqrt(2.*pi*x) t 'alpha = 1 (theory)',      \
      '../dispersion-alpha/log' u ($1/(2.*pi)):($3/sqrt($1))       \
           t 'alpha = 3.22 (numerical)',                           \
      omega(2.*pi*x,3.22)/sqrt(2.*pi*x) t 'alpha = 3.22 (theory)', \
      '../dispersion-gn/log' u ($1/(2.*pi)):($3/sqrt($1))          \
            t 'Green-Naghdi (numerical)',			   \
      sqrt(tanh(2.*pi*x)/(2.*pi*x)) t 'sqrt(k tanh(k h)/h)',       \
      '../dispersion-sv/log' u ($1/(2.*pi)):($3/sqrt($1))          \
            t 'Saint-Venant (numerical)'
~~~
*/
