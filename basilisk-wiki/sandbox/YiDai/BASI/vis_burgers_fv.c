/**
# solving the viscid burgers equation using finite volume
$$\partial_t U + U \partial_x U = \nu \partial_{xx}U$$
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar u[], du[], flux[];
double dt;
double BT = 1;
double mu = 0.001;

int main(){
    L0 = 1;
    X0 = 0.;
    N = 1 << 8;
    DT = L0/N;            // smaller than von neumann stability
    run();
}

// flux[left] = neumann(0.);
// flux[right] = neumann(0.);

u[left] = dirichlet(0.);
u[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
        // u[] = cos(2 * pi * x);
        u[] = sin(2 * pi * x);
        // u[] = exp(-sq(2 * x-3));
    }
    boundary ({u});
}

// do note that the printed data are vary with different DT in this case

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata (t = 0) {
  static FILE * fp = fopen ("case_FV_T0","w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, u[], t);
  fprintf (fp, "\n\n");
  fflush (fp);
}

event printdata2 (t = BT) {
  static FILE * fp1 = fopen ("case_FV_T1","w");
  foreach()
    fprintf (fp1, "%g %g %g\n", x, u[], t);
  fprintf (fp1, "\n\n");
  fflush (fp1);
}

event printdata3 (t = 0; t <= BT; t += 0.01) {
  static FILE * fp3 = fopen ("case_FV_TS","w");
  foreach()
    fprintf (fp3, "%g %g %g\n", x, u[], t);
  fprintf (fp3, "\n\n");
  fflush (fp3);
}

event integration (i++) {
  // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext (dt);
  // Godunov scheme
  foreach(){
      flux[] = 0.5 * sq(u[]);
  }
//   foreach(){
//     if (u[-1]>=u[0]){
//         flux_inter[] = max(flux[-1], flux[0]);
//     }
//     else if(u[-1]<=0 && u[0]>=0){
//         flux_inter[] = 0;
//     }
//     else{
//         flux_inter[] = min(flux[-1], flux[0]);
//     }
//   }
    foreach(){
        // flux[] = u[]<0? 0.5 * sq(u[1]): 0.5 * sq(u[]);
        flux[] = (flux[1] + flux[0])/2;
        // flux[] = u[-1]>u[0]? max(flux[-1], flux[0]):min(flux[-1], flux[0]);
    }
    
    // flux[] = u[-1]>u[0]? max(flux[-1], flux[0]):min(flux[-1], flux[0]);
    // flux[] = u[]<0? 0.5 * sq(u[1]): 0.5 * sq(u[]);
    // flux[] = u[]<0? 0.5 * sq(u[1]): 0.5 * sq(u[]);

//   boundary ({flux_inter});
  foreach()
    du[] = (-flux[0] + flux[-1])/Delta + mu * (u[1] - 2 * u[] + u[-1])/sq(Delta);
  foreach()
    u[] += dt*du[];
  boundary ({u});
}

/**
~~~gnuplot
data1 = "case_FV_T0"
data2 = "case_FV_T1"
set xlabel "x"
set ylabel "u"
p[][] data1 u ($1):($2) t "T0" w l,\
 data2 u ($1):($2) t "T1" w l

~~~
*/
/**
I don't understand why the solution become unstable when $\nu$ = 0.01 
*/
