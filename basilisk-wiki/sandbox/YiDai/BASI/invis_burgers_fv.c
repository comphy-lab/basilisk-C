/**
# solving inviscid burgers equation using finite volume 
The equation 
$$\partial_t U + U \partial_x U = 0$$

can be written in conservative form
$$\partial_t U + \partial_x (\frac{1}{2}U^{2}) = 0$$

Finite volume method describes the variable in integral form and approximates the first term using slab average and second term interface flux $F(U) = \frac{1}{2}U^{2}$. The approximation of the flux depends on moving direction of the shock wave if there is one in the case of burgers equation. 

$$\frac{d}{d t} \int_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} U d x+\left.F(U)\right|_{x_{i+\frac{1}{2}}}-\left.F(U)\right|_{x_{i-\frac{1}{2}}}=0$$

In the presence of shock wave, the flux is always taken as the upwind scheme. 

*/


#include "grid/cartesian1D.h"
#include "run.h"

scalar u[], du[], flux[];
double dt;
double BT = 0.3;

int main(){
    L0 = 1;
    X0 = 0.;
    N = 1 << 8;
    DT = L0/N;            // smaller than von neumann stability
    run();
}

u[left] = dirichlet(0.);
u[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
        // u[] = cos(2 * pi * x);
        u[] = sin(2 * pi * x);
    }
    boundary ({u});
}

// do note that the printed data are vary with different DT in this case

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata (t = 0) {
  static FILE * fp = fopen ("case1_FV_T0","w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, u[], t);
  fprintf (fp, "\n\n");
  fflush (fp);
}

event printdata2 (t = BT) {
  static FILE * fp1 = fopen ("case1_FV_T1","w");
  foreach()
    fprintf (fp1, "%g %g %g\n", x, u[], t);
  fprintf (fp1, "\n\n");
  fflush (fp1);
}

event printdata3 (t = 0; t <= BT; t += 0.01) {
  static FILE * fp3 = fopen ("case1_FV_TS","w");
  foreach()
    fprintf (fp3, "%g %g %g\n", x, u[], t);
  fprintf (fp3, "\n\n");
  fflush (fp3);
}

event integration (i++) {
  // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext (dt);
  // #if upwind nonconservative
  foreach()
    flux[] = u[]<0? 0.5 * sq(u[1]): 0.5 * sq(u[]);
  boundary ({flux});
  foreach()
    du[] = (-flux[] + flux[-1])/Delta;
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
p[][-1.1:1.1] data1 u ($1):($2) t "T0" w l,\
 data2 u ($1):($2) t "T1" w l
~~~
*/
/**
with different initial condition
~~~gnuplot
data1 = "case1_FV_T0"
data2 = "case1_FV_T1"
set xlabel "x"
set ylabel "u"
p[][-1.1:1.1] data1 u ($1):($2) t "T0" w l,\
 data2 u ($1):($2) t "T1" w l
~~~
*/
