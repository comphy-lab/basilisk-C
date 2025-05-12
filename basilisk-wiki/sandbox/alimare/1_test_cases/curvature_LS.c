/**
#Curvature calculation

Simple test case for curvature calculation, we make a comparison between the
curvature calculated using a level-set field and a VOF field.
*/

#include "../alex_functions.h"
#include "fractions.h"
#include "curvature.h"
#include "view.h"

#define Pi 3.14159265358979323846

double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.;
  center.y = 0.;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = -( sqrt(R2)*(1+0.3*cos(6*theta)) - Radius);


  return s;
}


scalar cs[];
face vector fs[];
#include "../LS_curvature.h"

int main(){
  L0 = 2.;
  origin (-L0/2.,-L0/2.);
  for(int k = 5; k <=8; k++){
    int N = 1 << k;
    init_grid(N);
  
    scalar LS[];
  
    foreach_vertex(){
      LS[] = -geometry(x,y,L0/3.);
    }
  
    fractions (LS, cs, fs);

    scalar curve_LS[];
    curvature_LS(LS, curve_LS);
    stats s  = statsf(curve_LS);

    scalar curve[];
    curvature(cs,curve);
    boundary({curve});
    stats s2 = statsf(curve);

    scalar err[];
    foreach(){
      if(interfacial(point,cs)){
        err[] = fabs(curve[] - curve_LS[]);
      }
      else{
        err[] = nodata;
      }
    }
    boundary({err});
    stats s3 = statsf(err);  
    if(k == 7){
      draw_vof("cs");
      squares("err");
      save("err.png");
    }
    fprintf(stderr, "%d %g %g %g %g %g %g %g %g %g\n", 1<< k, s.min, s.max,
      s.stddev,s2.min,
     s2.max,s2.stddev,s3.min,s3.max,s3.stddev);
  }
}

/**
![Difference between the two curvatures](curvature_LS/err.png)


~~~gnuplot Curvature*
set key bottom right
set logscale
set xrange [16:512]
set xtics 16,2,512
set grid
plot 'log' u 1:4 w l t 'stddevLS', '' u 1:7 w l  t 'stddevVOF', \
  '' u 1:10 w l t 'difference-stddev' 
~~~

Both method show comparable results.
I should try to run discrimining results.
*/