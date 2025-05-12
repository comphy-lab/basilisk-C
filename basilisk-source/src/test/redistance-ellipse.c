/**
# Redistancing of a perturbed distance field

Test case originally written by [Alexandre Limare](/sandbox/alimare/README).

This case is extracted from [Russo et al., 2000](#russo2000). We
initialize a perturbed distance field, where the zero levelset is an
ellipse of the form
$$
\phi (x,y,0) = f(x,y) \times g(x,y)
$$
where the perturbation is
$$
f(x,y) = \epsilon  + (x - x_0)^2 +(y - y_0)^2
$$
and the ellipse is
$$
g(x,y) = \left( \sqrt{\frac{x^2}{A^2}+\frac{y^2}{B^2}} -R \right)
$$
with $A=2$, $B=1$, $R = 1$ , $x_0 = 3.5$, $y_0 = 2.$.

We want to recover a perfect distance field, *i.e.* remove the
initial perturbation. 

We show here the initial and final level-set for the same isovalues.

![Initial distance field](redistance-ellipse/start.png) 

![Final distance field](redistance-ellipse/final.png)

![Error, logscale between $10^{-4}$ and $10^{-1}$](redistance-ellipse/err.png)

~~~gnuplot Error analysis
f(x) = a + b*x
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b
fit f1(x) 'log' u (log($1)):(log($3)) via a1,b1
fit f2(x) 'log' u (log($1)):(log($4)) via a2,b2
ftitle(a,b) = sprintf("%.2f-%4.2f x", a, -b)
set logscale xy
set xrange [32:512]
set xtics 32,2,512
set format y "%.1e"
plot 'log' u 1:2 t 'avg', exp(f(log(x))) t ftitle(a,b), \
    'log' u 1:3 t 'rms', exp(f1(log(x))) t ftitle(a1,b1), \
    'log' u 1:4 t 'max', exp(f2(log(x))) t ftitle(a2,b2)
~~~

~~~gnuplot Error analysis 0-level-set
f(x) = a + b*x
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
unset logscale
unset xrange
fit f(x) 'log' u (log($1)):(log($5)) via a,b
fit f1(x) 'log' u (log($1)):(log($6)) via a1,b1
fit f2(x) 'log' u (log($1)):(log($7)) via a2,b2
ftitle(a,b) = sprintf("%.2f-%4.2f x", a, -b)
set logscale xy
set xrange [32:512]
set xtics 32,2,512
set format y "%.1e"
plot 'log' u 1:5 t 'avg', exp(f(log(x))) t ftitle(a,b), \
    'log' u 1:6 t 'rms', exp(f1(log(x))) t ftitle(a1,b1), \
    'log' u 1:7 t 'max', exp(f2(log(x))) t ftitle(a2,b2)
~~~

Here we study the value of the level-set function on a set of points where it is
theoretically 0, we show that we have also a $3^{rd}$ convergence.

## References

~~~bib
@article{russo2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}
~~~
*/

#include "distance_point_ellipse.h"
#include "redistance.h"
#include "view.h"

double perturb (double x, double y, double eps, coord center)
{
  return eps + sq(x - center.x) + sq(y - center.y);
}

void draw_isolines (scalar s, double smin, double smax, int niso, int w)
{
  vertex scalar vdist[];
  foreach_vertex()
    vdist[] = (s[] + s[-1] + s[0,-1] + s[-1,-1])/4.;
  for (double sval = smin ; sval <= smax; sval += (smax - smin)/niso)
    isoline ("vdist", sval, lw = w);
}

int main()
{
  origin (-5., -5.);
  L0 = 10;

  for (int MAXLEVEL = 6; MAXLEVEL < 9; MAXLEVEL++){
    init_grid (1 << MAXLEVEL);  
  
    scalar dist[];
    double A = 4., B = 2.;
    coord  center_perturb = {3.5, 2.};
    foreach(){
      double a, b;
      dist[] = DistancePointEllipse (A, B, x, y, &a, &b)*
	perturb (x, y, 0.1, center_perturb);
    }

    if (MAXLEVEL == 8) {
      squares ("dist", map = cool_warm, min = -2, max = 2);
      draw_isolines (dist, -2, 2, 20, 1);
      save ("start.png");
    }
      
    int nbit = redistance (dist, imax = 1 << (MAXLEVEL + 1));
    
    if (MAXLEVEL == 8) {
      squares ("dist", map = cool_warm, min = -2, max = 2);
      draw_isolines (dist, -2, 2, 20, 1);
      save ("final.png");
    }
    
    scalar err[], errt[], errt1[];
    foreach() {
      double a, b;
      err[] = dist[] - DistancePointEllipse (A, B, x, y, &a, &b);
      errt[] = dist[] > - 0.8 ? err[] : nodata;
      errt1[] = fabs (dist[]) < 1.2*L0/(1 << grid->maxdepth) ? err[] : nodata;
    }

    norm n = normf (errt), n2 = normf (errt1);
    fprintf (stderr, "%d %g %g %g %g %g %g %d\n",
	     1 << MAXLEVEL, n.avg, n.rms, n.max,
	     n2.avg, n2.rms, n2.max, nbit);

    if (MAXLEVEL == 8) {
      squares ("log(fabs(err) + 1e-16)", min = log(n.max) - 6, max = log(n.max));
      save ("err.png");
    }
  }
}
