/**
#Level set redistancing 1D test case

This is a test case for the LS_reinit() taken from [Russo et al.,1999](#russo_remark_2000) we
initialize a perturbed distance field.
*/
#define BICUBIC 1
#include "embed.h"
#include "curvature.h"
#include "alimare/alex_functions.h"
#include "alimare/LS_reinit.h"
#include "view.h"

/**
The initial distance is a parabola :
$$
y = (x-0.4\Delta)*(x+6)/2+1
$$

We want to study how the use of LS_reinit() modifies the position of the 0
level-set.
*/

double initial_dist(double Delta, double x){
  return (x-0.4*Delta)*(x+6)/2.+1. ;
}

scalar dist[];
scalar * level_set = {dist};

int main() {
/**
The domain size is 10. Our mesh is made of 32*32 cells.
*/
  origin (-5., -5.);
  L0 = 10;

  int MAXLEVEL = 5;  
  init_grid (1 << MAXLEVEL);

  foreach(){
    dist[] = initial_dist(Delta, y);
  }
  boundary({dist});

/**
The position of the 0 level-set will be given using the embedded boundaries and
the height function.
*/
  scalar dist_n[];
  cell2node(dist,dist_n);
  fractions (dist_n, cs, fs);
  double y_max = -5.;

  vector h[];
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs)){
      double yy = y+Delta*height(h.y[]);
      y_max = max(y_max,yy);
    }     
  }
  double LS0 = y_max;

    
  char filename [100];
  FILE * fp1;  
  snprintf(filename, 100,  "log0");
  fp1 = fopen (filename,"w");
  foreach(){
    if(x==Delta/2.){
      fprintf(fp1, "%g %g %g\n", y, dist[], 0.);
    }
  }
  fclose(fp1);

  LS_reinit(dist, it_max = 1 << (MAXLEVEL+1));

  snprintf(filename, 100,  "log1");
  fp1 = fopen (filename,"w");
  foreach(){
    if(x==Delta/2.){
      fprintf(fp1, "%g %g\n", y, dist[]);
    }
  }
  fclose(fp1);

  cell2node(dist,dist_n);
  fractions (dist_n, cs, fs);
  y_max = -5.;
  heights (cs, h);
  boundary((scalar *){h});
  foreach(reduction(max:y_max)){
    if(interfacial(point, cs)){
      double yy = y+Delta*height(h.y[]);
      y_max = max(y_max,yy);
    }     
  }
  double LS1 = y_max;
  fprintf(stderr, "%g %g\n", LS0 , LS1);
  squares("dist", min = -0.1, max = 0.1);
  draw_vof("dist_n");
  save("dist.png");
}


/**


~~~gnuplot Evolution of the dist function
plot 'log1' u 1:2 w lp lc 'blue' lw 2 t 'final distance', \
     'log0' u 1:2 w lp lc 'red' lw 2 t 'initial distance', \
     'log0' u 1:3 w lp lc 'black' t ''
~~~

As expected the redistance function gives us a linear distance function from the
y=0 value which seem to be correctly preserved.


~~~gnuplot Zoom. Evolution of the dist function
K = "`tail -1 log | awk '{print $1}'`"
K2 = "`tail -1 log | awk '{print $2}'`"
K = K + 0
K2 = K2 + 0
graph(n)  = sprintf("Initial, x = %.3e",n)
graph2(n) = sprintf("Final, x = %.3e",n)
set xrange[-0.7:0.5]
plot 'log1' u 1:2 w lp lc 'blue' lw 2 t graph2(K2), \
     'log0' u 1:2 w lp lc 'red' lw 2 t graph(K), \
     'log0' u 1:3 w lp lc 'black' t ''
~~~

After a little analysis, on can show that there is a 0.01 discrepency between
the position of the two 0-level-sets which is acceptable.

The two images can be compared with Figs. 5 and 6 of the paper by Russo
and Smereka. The main differences between our case and theirs is that we use a
cell-centered level-set and 32 grid cells and they used a vertex-centered
level-set function and 20 grid cells.


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