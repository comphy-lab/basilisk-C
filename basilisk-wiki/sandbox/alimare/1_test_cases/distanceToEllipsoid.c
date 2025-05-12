/**
#LS_reinit() test case in 3D.


This case is extracted from [Russo et al.,1999](#russo_remark_2000) we
initialize a perturbed distance field, where the zero level-set is an ellipsoid of
the form:
$$
\phi (x,y,z,t=0) = f(x,y,z) \times g(x,y,z)
$$
where g(x,y,z) is the actual distance to the ellipsoid and the perturbation is :
$$
f(x,y) = \epsilon  + (x - x_0)^2 +(y - y_0)^2 + (z - z_0)^2
$$
with $x_0 = 3.5$, $y_0 = 2.$, $z_0 = 1.$ and $\epsilon = 0.1$.

Here the considered ellispoid $\Gamma(\boldsymbol{x}) = \phi(\boldsymbol{x}) = 0$ is given by:
$$
\Gamma(\boldsymbol{x}) = \left( \sqrt{\frac{x^2}{A^2}+\frac{y^2}{B^2}+\frac{z^2}{C^2}} -R \right)
$$
with $A=4$, $B=2$, $C=1$, $R = 1$.

We want to recover a perfect distance field, \textit{i.e.} remove the initial
perturbation.

*/
#define QUADRATIC 1
#define BGHOSTS 2
#include "grid/octree.h"
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#include "../distance_point_ellipse3D.h"
#include "../alex_functions.h"
#include "../LS_reinit.h"
#include "../basic_geom.h"
#include "view.h"

double perturb (double x, double y, double z ,double eps, coord center){
  return eps + sq(x - center.x) + sq(y - center.y) + sq(z- center.z);
}




#define Pi 3.141592653589793

norm mynorm(scalar f, scalar f2, double threshold){
/**
For this test case, we use a special norm defined by the author, which ignores
values above a specific threshold, $f > \text{threshold}$ 
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


int main() {
  origin (-5., -5., -5.);
  L0 = 10;
 
  for(int MAXLEVEL = 5; MAXLEVEL <= 7; MAXLEVEL++){
    init_grid (1 << MAXLEVEL);  
  
    double A = 4., B=2., C=1.;

    // A > B > C, mandatory for the algo to work
    assert(A>B);
    assert(B>C);
    
    coord  center_perturb = {3.5,2.,1.};
    scalar dist[],distinit[];
    foreach(){
      double a,b,c;
      dist[] = DistancePointEllipsoid(A, B, C, x, y, z, &a, &b, &c)*
      perturb(x,y,z, 0.1, center_perturb);
      distinit[] = dist[];
    }
    boundary({dist,distinit});

 /**
And now we reinit
**/
    int nbit = LS_reinit(dist, it_max = 1 << (MAXLEVEL+1));
    scalar err[],LogErr[];
    foreach(){
      double a,b,c;
      err[] = dist[] - DistancePointEllipsoid(A, B, C, x, y, z, &a, &b, &c);
      LogErr[] = log(fabs(err[])+1.e-16);
    }
    boundary({err,LogErr});
    norm n = mynorm(err,dist,-0.8);
    norm n2 = mynorm2(err,dist,1.2*L0/(1 << grid->maxdepth));
    fprintf(stderr, "%d %g %g %g %g %g %g %d\n",1<<MAXLEVEL, n.avg, n.rms,n.max,
      n2.avg, n2.rms,n2.max, nbit);
/** 
We will display several isocontour using the level set ls1 variable and the
draw_vof function on dist2 and dist3 the vertex centered scalar before and
after the redistanciation process.
The iso-surfaces displayed will be for values with a step of 0.2, 
*/
    double ls1 = 0.;
    char name[80];
    int count = 0; // to alternate between planes
    if(MAXLEVEL == 6){
        view(quat = {0.,0.,0.,1.});
        sprintf (name, "error-plane%d.png", ++count);
        squares ("LogErr", min = log(n.max)-6, max = log(n.max));
        save(name);
        view(quat = {-0.707107,-0,-0,0.707107});
        sprintf (name, "error-plane%d.png", ++count);
        squares ("LogErr", min = log(n.max)-6, max = log(n.max), n = {0,1,0});
        save(name);
        view(quat = {-0,-0.707107,-0,0.707107});
        sprintf (name, "error-plane%d.png", ++count);
        squares ("LogErr", min = log(n.max)-6, max = log(n.max), n = {1,0,0});
        save(name);
    }

    for (int j = 0; j < 5; j++){

      count = 0;
      view(quat = {0.,0.,0.,1.});
      sprintf (name, "distinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("distinit", ls1);
      squares ("distinit", map = cool_warm, min = -2, max = 2);
      save(name);

      view(quat = {-0.707107,-0,-0,0.707107});
      sprintf (name, "distinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("distinit", ls1);
      squares ("distinit", map = cool_warm, min = -2, max = 2, n = {0,1,0});
      save(name);
      
      view(quat = {-0,-0.707107,-0,0.707107});
      sprintf (name, "distinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("distinit", ls1);
      squares ("distinit", map = cool_warm, min = -2, max = 2, n = {1,0,0});
      save(name);

      count = 0;
      view(quat = {0.,0.,0.,1.});
      sprintf (name, "distreinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("dist", ls1);
      squares ("dist", map = cool_warm, min = -2, max = 2);
      save(name);
      view(quat = {-0.707107,-0,-0,0.707107});
      sprintf (name, "distreinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("dist", ls1);
      squares ("dist", map = cool_warm, min = -2, max = 2, n = {0,1,0});
      save(name);
      view(quat = {-0,-0.707107,-0,0.707107});
      sprintf (name, "distreinit-plane%d-LS%g.png", ++count,ls1);
      isosurface ("dist", ls1);
      squares ("dist", map = cool_warm, min = -2, max = 2, n = {1,0,0});
      save(name);

      ls1+=0.2; // increment to the level-set value selected
    }
  }
}


/**
We show here the initial and final level-set for the same isovalues in front top and left planes.


|  Field   | Plane1   |      Plane2      |  Plane3 |
|:-------------:|:-------------:|:-------------:|:-------------:|
| Initial dist = 0   | ![](distanceToEllipsoid/distinit-plane1-LS0.png)     | ![](distanceToEllipsoid/distinit-plane2-LS0.png)     | ![](distanceToEllipsoid/distinit-plane3-LS0.png)
| Initial dist = 0.2 | ![](distanceToEllipsoid/distinit-plane1-LS0.2.png)   | ![](distanceToEllipsoid/distinit-plane2-LS0.2.png)   | ![](distanceToEllipsoid/distinit-plane3-LS0.2.png)
| Initial dist = 0.4 | ![](distanceToEllipsoid/distinit-plane1-LS0.4.png)   | ![](distanceToEllipsoid/distinit-plane2-LS0.4.png)   | ![](distanceToEllipsoid/distinit-plane3-LS0.4.png)
| Initial dist = 0.6 | ![](distanceToEllipsoid/distinit-plane1-LS0.6.png)   | ![](distanceToEllipsoid/distinit-plane2-LS0.6.png)   | ![](distanceToEllipsoid/distinit-plane3-LS0.6.png)
| Initial dist = 0.8 | ![](distanceToEllipsoid/distinit-plane1-LS0.8.png)   | ![](distanceToEllipsoid/distinit-plane2-LS0.8.png)   | ![](distanceToEllipsoid/distinit-plane3-LS0.8.png)
| Final dist = 0     | ![](distanceToEllipsoid/distreinit-plane1-LS0.png)   | ![](distanceToEllipsoid/distreinit-plane2-LS0.png)   | ![](distanceToEllipsoid/distreinit-plane3-LS0.png)
| Final dist = 0.2   | ![](distanceToEllipsoid/distreinit-plane1-LS0.2.png) | ![](distanceToEllipsoid/distreinit-plane2-LS0.2.png) | ![](distanceToEllipsoid/distreinit-plane3-LS0.2.png)
| Final dist = 0.4   | ![](distanceToEllipsoid/distreinit-plane1-LS0.4.png) | ![](distanceToEllipsoid/distreinit-plane2-LS0.4.png) | ![](distanceToEllipsoid/distreinit-plane3-LS0.4.png)
| Final dist = 0.6   | ![](distanceToEllipsoid/distreinit-plane1-LS0.6.png) | ![](distanceToEllipsoid/distreinit-plane2-LS0.6.png) | ![](distanceToEllipsoid/distreinit-plane3-LS0.6.png)
| Final dist = 0.8   | ![](distanceToEllipsoid/distreinit-plane1-LS0.8.png) | ![](distanceToEllipsoid/distreinit-plane2-LS0.6.png) | ![](distanceToEllipsoid/distreinit-plane3-LS0.8.png)
| Final error        | ![](distanceToEllipsoid/error-plane1.png)            | ![](distanceToEllipsoid/error-plane2.png)            | ![](distanceToEllipsoid/error-plane3.png)

~~~gnuplot error analysis
set pointsize 1.2
f(x) = a + b*x
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b
fit f1(x) 'log' u (log($1)):(log($3)) via a1,b1
fit f2(x) 'log' u (log($1)):(log($4)) via a2,b2
ftitle(a,b) = sprintf("%.2f-%4.2f x", a, -b)
set logscale xy
set xrange [16:256]
set xtics 16,2,256
set format y "%.1e"
plot 'log' u 1:2 t 'avg' pt 5 ps 1.2, exp(f(log(x)))  lt 2 lw 1.2 t ftitle(a,b)   , \
    'log' u 1:3 t 'rms'  pt 6 ps 1.2, exp(f1(log(x))) lt 3 lw 1.2 t ftitle(a1,b1), \
    'log' u 1:4 t 'max'  pt 7 ps 1.2, exp(f2(log(x))) lt 4 lw 1.2 t ftitle(a2,b2)
~~~
Between $1^{rst}$ and $2^{nd}$ order convergence.

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
set xrange [16:256]
set xtics 16,2,256
set format y "%.1e"
plot 'log' u 1:5 t 'avg' pt 5 ps 1.2, exp(f(log(x)))  lt 2 lw 1.2 t ftitle(a,b), \
    'log' u 1:6 t 'rms'  pt 6 ps 1.2, exp(f1(log(x))) lt 3 lw 1.2 t ftitle(a1,b1), \
    'log' u 1:7 t 'max'  pt 7 ps 1.2, exp(f2(log(x))) lt 4 lw 1.2 t ftitle(a2,b2)
~~~

Here we study the value of the level-set function on a set of points where it is
theoretically 0, we show that we have between a $2^{nd}$ and $3^{rd}$-order 
convergence.

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

~~~
*/