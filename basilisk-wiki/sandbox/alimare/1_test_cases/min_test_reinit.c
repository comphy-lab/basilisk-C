/**
#LS_reinit() from Min2010

Circle perturbed. Test by [Min](#Min2010). Should show third order accuracy in
the vicinity of the 0-level-set and 2nd ordre accuracy elsewhere.

*/
#define BICUBIC 1
#define BGHOSTS 2
#include "../../popinet/distance_point_ellipse.h"
#include "../alex_functions.h"
#include "../LS_reinit.h"
#include "../basic_geom.h"
#include "view.h"

norm mynorm(scalar f, scalar f2, double threshold){
/**
For this test case, we use a special norm defined by the author, which ignores
the points below a specific threshold, $f > \text{threshold}$ 
*/
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
    reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0. && f2[] > threshold) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
    }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

norm mynorm2(scalar f, scalar f2, double threshold){
/**
This second norm selects points  such that $|f| < \text{threshold}$.
*/
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
    reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0. && fabs(f2[]) < threshold) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
    }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}


void draw_isolines(scalar s, double smin, double smax, int niso, int w){
  vertex scalar vdist[];
  cell2node(s,vdist);

  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso){
    isoline ("vdist", sval, lw = w);
  }
}

double distinit(double x, double y){
  return (sq(x-1) + sq(y-1) +0.1)*(sqrt(x*x+y*y)-1.);
}

double distcircle(double x, double y){
  return (sqrt(x*x+y*y)-1.);
}

scalar dist[];
scalar * level_set = {dist};

int main() {

  origin (-2., -2.);
  L0 = 4;

  int MAXLEVEL = 6;
  for(MAXLEVEL = 6; MAXLEVEL <= 8; MAXLEVEL++){
    init_grid (1 << MAXLEVEL);  
  
    foreach(){
      dist[] = distinit(x,y);
    }
    boundary({dist});
  
    squares ("dist", map = cool_warm, min = -1, max = 1);
    draw_isolines(dist, -0.5, 0.5, 4, 1);
    save("dist_init.png");

    // for (int k = 0; k<=2* (1 << MAXLEVEL); k++)
      LS_reinit(dist, it_max = 1 << (MAXLEVEL+1), RK2 = 1);

    squares ("dist", map = cool_warm, min = -1, max = 1);
    draw_isolines(dist, -0.5, 0.5, 4, 1);
    save("dist_reinit.png");


    scalar err[],LogErr[];
    foreach(){
      err[] = dist[] - distcircle(x,y);
      LogErr[] = log(fabs(err[])+1.e-16);
    }
    boundary({err,LogErr});
    norm n = mynorm(err,dist,-0.8);
    norm n2 = mynorm2(err,dist,1.2*L0/(1 << grid->maxdepth));
    fprintf(stderr, "%d %g %g %g %g %g %g\n",1<<MAXLEVEL, n.avg, n.rms,n.max,
      n2.avg, n2.rms,n2.max);

    squares ("LogErr", min = log(n.max)-6, max = log(n.max));
    save("LogErr.png");
    
  }
}


/**
We show here the initial and final level-set for the same isovalues.

![Initial level-set](min_test_reinit/dist_init.png) 

![Redistanced field. That's all folks !](min_test_reinit/dist_reinit.png)

![Error - logscale](min_test_reinit/LogErr.png)

~~~gnuplot error analysis
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

~~~gnuplot error analysis 0-level-set
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

## References

~~~bib

@article{russo_remark_2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}

@article{Min2010,
  author        = {Chohong Min},
  title         = {On reinitializing level set functions},
  year          = {2010},
  volume        = {229},
  pages         = {2764-2772},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2009.12.032},
}


~~~
*/