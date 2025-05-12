/**
# solver 1D heat equation. 
$$\partial_{t}T = \partial_{xx}T$$

We tested several basilisk features in this code by tuning the settings. 

*/
#include "grid/cartesian1D.h"
#include "run.h"

// two initial condition
scalar T1[], dT1[]; 
scalar T2[], dT2[];
double dt;
#define EPS 0.1

int main(){
    L0 = 20;
    X0 = -L0/2.;
    N = 1 << 7;
    DT = sq(L0/N)/3;             // smaller than von neumann stability
    run();
}

event init(t = 0){
    foreach(){
        T1[] = 1./EPS*(fabs(x)<EPS)/2;
        T2[] = x>=0? 3: 0; 
    }
    boundary ({T1, T2});
}

// do note that the printed data are vary with different DT in this case

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata (t = 0; t <= 1; t += 0.2) {
  static FILE * fp = fopen ("output_H1DC","w");
  foreach()
    fprintf (fp, "%g %g %g %g %g\n", x, T1[], T2[], dT1[], t);
  fprintf (fp, "\n\n");
  fflush (fp);
}

event integration (i++) {
  // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext (dt);
  foreach(){
    dT1[] = (T1[1,0] - 2 *T1[0,0] + T1[-1,0])/sq(Delta); 
    dT2[] = (T2[1,0] - 2 *T2[0,0] + T2[-1,0])/sq(Delta); 
  }
  foreach(){
    T1[] += dt*dT1[];
    T2[] += dt*dT2[];
  }
  boundary ({T1, T2});
}


/**
~~~gnuplot
set terminal png size 800,600 enhanced font 'Times-Roman,16'
set key samplen 2 spacing 1.5 font 'Times-Roman,16'
file="H1DC"
set multiplot layout 1,2
set xlabel "x"
set ylabel "T"
p[][] 'output_'.file u ($1):($2) t 'IC1' w l
p[][-0.5:3.5] 'output_'.file u ($1):($3) t 'IC2' w l
unset multiplot
~~~
*/


