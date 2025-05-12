/**
# solving equation linear advection - diffusion equation
$$\partial_{t}u + c \partial_{x}u = \nu \partial_{xx}u$$
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar u[], du[];
double dt;
double BT = 1.;

int main(){
    L0 = 2;
    N = 1 << 8;
    DT = sq(L0/N)/3;            // smaller than von neumann stability
    run();
}

u[left] = dirichlet(0.);
u[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
        u[] = sin(2*pi*x);
    }
    boundary ({u});
}

// do note that the printed data are vary with different DT in this case

event printdata (t = 0; t <= 0.1; t += 0.01) {
  static FILE * fp = fopen ("cases0_c4","w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, u[], t);
  fprintf (fp, "\n\n");
  fflush (fp);
}


event integration (i++) {
  // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext (dt);
  // #if upwind nonconservative
  foreach()
    du[] = -4/(2*Delta)*(u[1]-u[-1]) + 0.5/sq(Delta)*(u[1]-2*u[]+u[-1]);
  foreach()
    u[] += dt*du[];

  boundary ({u});
}


/**
~~~gnuplot
reset
set xlabel "x"
set ylabel "u"
p[][] "cases0" u ($1):($2) t '' w l
~~~
*/
/**
~~~gnuplot
reset
set xlabel "x"
set ylabel "u"
set title "with BC"
p[][] "cases0_BC" u ($1):($2) t '' w l
~~~
*/
/**
~~~gnuplot
reset
set xlabel "x"
set ylabel "u"
set title "c = 4"
p[][] "cases0_c4" u ($1):($2) t '' w l
~~~
*/

