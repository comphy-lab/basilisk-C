/**
#Test of the my solvers

Save discretization in space... normally should be ENO3. To be tested.

A circle is advected diagonally and comes back to its original place (magic of
the periodic boundary conditions!).

We do plots of the error and show the final fields with isovalues in the
vicinity of the 0-level-set (what I'm interested in ...).


*/

/**
~~~gnuplot Error convergence
ftitle(a,b) = sprintf("%.3fx^\{%4.2f\}", exp(a), -b)

f(x) = a + b*x
fit f(x) 'log0' u (log($1)):(log($2)) via a,b

f2(x) = a2 + b2*x
fit f2(x) 'log1' u (log($1)):(log($2)) via a2,b2

f3(x) = a3 + b3*x
fit f3(x) 'log2' u (log($1)):(log($2)) via a3,b3

set ylabel 'Average error'
set xrange [16:256]
set yrange [*:*]
set xtics 16,2,256
set format y "%.1e"
set logscale
plot 'log0' u 1:2 pt 7 lc 'blue'  t 'FE ', exp(f (log(x))) t ftitle(a,b) lc 'blue',\
     'log1' u 1:2 pt 7 lc 'red'   t 'RK2', exp(f2(log(x))) t ftitle(a2,b2) lc 'red',\
     'log2' u 1:2 pt 7 lc 'green' t 'RK3', exp(f2(log(x))) t ftitle(a2,b2) lc 'green'

~~~

|  Time integration scheme   | FE   |      RK2      |  RK3 |
|:-------------:|:-------------:|:-------------:|:-------------:|
|  Final value   | ![final dist](testRK3/distFE.png)  | ![final dist](testRK3/distRK2.png)  | ![final dist](testRK3/distRK3.png)  |
|  Error plot   | ![error plot](testRK3/logeFE.png) |  ![error plot](testRK3/logeRK2.png) | ![error plot](testRK3/logeRK3.png)  |

![theoretical value of distance](testRK3/theoRK3.png)


*/

#define BGHOSTS 2
#define QUADRATIC 1
#include "utils.h"
#include "../basic_geom.h"
#include "../simple_discretization.h"
#include "../LS_reinit.h"
#include "view.h"

double perturb (double x, double y, double eps, coord center){
  return eps + sq(x - center.x) + sq(y - center.y);
}

void draw_isolines(scalar s, double smin, double smax, int niso, int w){
  vertex scalar vdist[];
  cell2node(s,vdist);
  
  boundary ({vdist});
  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso){
    isoline ("vdist", sval, lw = w);
  }
}

FILE * fp1;

int main() {
  origin (-L0/2.,-L0/2.);

  int j;
  double mytime;
  for(int jj=0;jj<=2;jj++){ // loop on method of integration
    char fname[100];
    snprintf(fname, 100,  "log%d", jj);
    fp1 = fopen(fname,"w");
    for (j=0;j<=2;j++){
      mytime = 0.;
      int MAXLEVEL;
      MAXLEVEL = 5+j;

      periodic(right);
      periodic(top);

      double nb_cell_NB  = 1<<3 ;
      double NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
      init_grid (1 << MAXLEVEL);

      double mydt = 4.e-3/pow(2.,j);
      coord center = {0.,0.};
      double R_init = 0.2;
      scalar dist[];
      foreach(){
        dist[] = circle(x,y,center,R_init);
      }
      boundary({dist});
      restriction({dist});
      vector vpc[];
      foreach(){
        vpc.x[] = 1./sqrt(2.);
        vpc.y[] = 1./sqrt(2.);
      }
      boundary((scalar *){vpc});
      restriction((scalar *){vpc});

/**
Time loop
*/
      
      for(int k = 0; mytime<sqrt(2.);k++){
        mytime += mydt;
        if(jj==0){
          scalar disti[];
          foreach(){
            disti[] = dist[];
          }
          boundary({disti});
          restriction({disti});
          FE(dist,disti,vpc,mydt, NB_width);
        }
        if(jj==1)
          RK2(dist,vpc,mydt, NB_width);
        if(jj==2)
          RK3(dist,vpc,mydt, NB_width);

        LS_reinit(dist, it_max = 10); // uses an RK3 time integration also.
      }
/**
error plot
*/
      center.x = 0.;
      center.y = 0.;
      scalar theo[], e[],loge[];
      foreach(){
        theo[] = circle(x,y,center,R_init);
        if(fabs(dist[])< NB_width){
          e[] = dist[] - theo[];
          loge[] = log(fabs(e[]));
        }
        else{
          e[] = nodata;
          loge[] = nodata;
        }
      }
      norm n = normf(e);
      fprintf(fp1, "%d %g %g\n", 1<<MAXLEVEL, n.avg, n.max);
      if(j==2){
        char filename [100],appendix [100];
        if(jj==0)
          strcpy(appendix,"FE.png");
        if(jj==1)
          strcpy(appendix,"RK2.png");
        if(jj==2)
          strcpy(appendix,"RK3.png");
        
        squares("loge", min = -10,max = -6);
        strcpy(filename, "loge");
        strcat(filename, appendix);
        save(filename);
        draw_isolines(theo, -0.2, 0.2, 20, 1);
        squares("theo");
        strcpy(filename, "theo");
        strcat(filename, appendix);
        save(filename);
        draw_isolines(dist, -0.2, 0.2, 20, 1);
        squares("dist");
        strcpy(filename, "dist");
        strcat(filename, appendix);
        save(filename);
      }
    }
    fclose(fp1);
  }
}