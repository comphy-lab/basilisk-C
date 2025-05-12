/**
#Calculation of the gradient on a sphere


~~~gnuplot Stats of the 2-norm of the gradient on the interface
set key left
set logscale x 2
set xrange [16:256]
set yrange [8:12]
f(x) = 10
plot 'log' u 1:2 lc rgb 'black' t 'plane area center min', '' u 1:3 lc rgb 'black' t 'max plane\_center', '' u 1:4 lc rgb 'black' t 'avg',\
  f(x) t 'theory'
~~~

![Initial temperature](sphere_gradient3D/T_init.png) 

![Initial gradient on the interface, 2-norm](sphere_gradient3D/norm_gtr_init.png) 

*/

#include "grid/octree.h"
#include "run.h"
#include "embed.h"
#include "view.h"
#define QUADRATIC 1
#include "../alex_functions.h"

scalar T[];

double T_init(double val,double V){
    // return -1.+exp(-V*val));
    return -val*V;
}

int main(){
	T[embed] = dirichlet(0.);
  T.third = true;
  for (int j = 0; j < 3; ++j){
    init_grid(1 << (5+j));
    origin(-L0/2., -L0/2., -L0/2.);
    fprintf(stderr, "##### %d\n", j);
    run();
  }
}

event init(t= 0){
  scalar dist[];
  coord center = {0.,0.,0.};
  double Radius =L0/(1<<3);
  foreach(){
    double r = sq(x-center.x)+ sq(y-center.y) + sq(z-center.z);
    dist[] = sqrt(r) - Radius;
  }

  boundary({dist});
  restriction({dist});
  vertex scalar distn[];
  cell2node(dist,distn);
  fractions (distn, cs, fs);
  foreach(){
   T[] = T_init(dist[],10);
  }
  boundary((scalar *){T});
  restriction((scalar *){T});

/**
Calculation of the gradient on the interface
*/
  vector gtr[],gtr2[];
  foreach(){
    foreach_dimension(){
      gtr.x[] = 0.;
      gtr2.x[] = 0.;
    }
    if(cs[]*(1-cs[]) != 0.){
      coord n, p;
      embed_geometry(point, &p, &n);
      double c    = 0.;
      double temp = 0.;
      double grad = dirichlet_gradient(point, T, cs , n, p, 
        temp, &c);

/**
For degenerate cases (stencil is too small), we add the term $tr[]*c$ to the 
gradient.
*/
      foreach_dimension(){
        gtr.x[] += grad*n.x+T[]*c;
      }
    }
  }
  boundary((scalar*){gtr});

  scalar norm_gtr[];
  foreach(){
    norm_gtr[] = 0.;
    foreach_dimension(){
      norm_gtr[] += sq(gtr.x[]);
    }
    norm_gtr[] = sqrt(norm_gtr[]);
    if(norm_gtr[] == 0.)norm_gtr[] = nodata;
  }
  boundary({norm_gtr});
  restriction({norm_gtr});

  stats s = statsf(norm_gtr);
  norm n = normf(norm_gtr);
  draw_vof("cs","fs");
  squares("T");
  save("T_init.png");
  squares("norm_gtr", min = 9.8, max = 10);
  save("norm_gtr_init.png");
  fprintf(stderr, "%d %g %g %g %g \n",1<<(grid->maxdepth), 
    s.min, s.max, n.avg, n.rms);
}


