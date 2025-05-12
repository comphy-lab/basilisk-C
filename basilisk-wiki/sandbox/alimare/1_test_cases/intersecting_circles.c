/**
#Reinitialization of an interface with kinks

Taken from [Min](#Min2010).

## References

~~~bib

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
#define BICUBIC 1
#define BGHOSTS 2
#include "../../popinet/distance_point_ellipse.h"

#include "embed.h"
#include "../alex_functions.h"
#include "../LS_reinit.h"
#include "view.h"

double perturb (double x, double y, double eps, coord center){
  return eps + sq(x - center.x) + sq(y - center.y);
}

double init_circle(double x, double y, double a, double r){
  if((a-x)/sqrt(sq(a-x)+sq(y))>= a/r &&
     (a+x)/sqrt(sq(a+x)+sq(y))>= a/r){
    return -min(sqrt(sq(x) + sq(y + sqrt(sq(r) - sq(a))) ),
               sqrt(sq(x) + sq(y - sqrt(sq(r) - sq(a))) )) ;
  }
  else
  return min(sqrt(sq(x+a)+sq(y))-r,
   sqrt(sq(x-a)+sq(y))-r) ;
}

void draw_isolines(scalar s, double smin, double smax, int niso, 
  int w){
  scalar vdist[];
  cell2node(s,vdist);

  boundary ({vdist});
  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso){
    isoline ("vdist", sval, lw = w);
  }
}

norm mynorm(scalar f, scalar f2, double threshold){
/**
For this test case, we use a special norm defined by the author, which ignores
the points below a specific threshold, $f > \text{threshold}$ 
*/
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
    reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0. && fabs(f2[]) > threshold) {
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

scalar dist[];
scalar * level_set = {dist};

int main() {
  origin (-2., -2.);
  L0 = 4;

  for(int j = 0;j<=2;j++){  
    int MAXLEVEL = 6+j;  
    init_grid (1 << MAXLEVEL);

    coord  center_perturb = {1.,1.};
    double a=0.7, R=1;
    foreach(){
      // dist[] = init_circle(x,y,a,R);
      dist[] = init_circle(x,y,a,R)*
      perturb(x,y,0.1,center_perturb);
    }
    boundary({dist});

    view (fov = 20.);
    squares ("dist", map = cool_warm, min = -2, max = 1);
    draw_isolines(dist, -2., 1., 30, 1);
    save("dist_init.png");

    LS_reinit(dist, it_max = 1 << (MAXLEVEL+1));
    squares ("dist", map = cool_warm, min = -2, max = 1);
    draw_isolines(dist, -2., 1., 30, 1);
    save("dist_first_reinit.png");
    scalar err[],LogErr[];
    foreach(){
      err[] = dist[] - init_circle(x,y,a,R);
      LogErr[] = log(fabs(err[])+1.e-16);
    }
    boundary({err,LogErr});
    norm n = mynorm(err,dist,0.1);
    norm n2 = mynorm2(err,dist,1.2*L0/(1 << grid->maxdepth));
    fprintf(stderr, "%d %g %g %g %g %g %g\n",1<<MAXLEVEL, n.avg, n.rms,n.max,
      n2.avg, n2.rms,n2.max);

  }
}


/**
We show here the initial and final level-set for the same isovalues.

![Initial level-set](intersecting_circles/dist_init.png) 

![first reinit level-set](intersecting_circles/dist_first_reinit.png) 

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

Order of accuracy is 1 globally

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

Order of accuracy between 0.5 and 2.

*/