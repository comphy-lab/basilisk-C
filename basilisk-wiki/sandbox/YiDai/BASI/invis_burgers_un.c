/**
# solving inviscid burgers equation using different scheme. 
$$\partial_{t}U + U\partial_{x}U = 0$$
We solver this equation using upwind-non-conservative, Lax_Friedrichs, Lax_Wendroff scheme.
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar u[], du[];
double dt;
double BT = 1.;   // The breaking time  = 1/ min(f'(x))

int main(){
    L0 = 2*pi;
    X0 = -L0/2.;
    N = 1 << 8;
    DT = L0/N;            // smaller than von neumann stability
    run();
}

u[left] = dirichlet(0.);
u[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
      //  u[] = sin(2*pi*x);
      u[] = cos(x)>0? cos(x):0;
    }
    boundary ({u});
}

// do note that the printed data are vary with different DT in this case

// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata (t = 0) {
  static FILE * fp = fopen ("ini_UW_T1","w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, u[], t);
  fprintf (fp, "\n\n");
  fflush (fp);
}

event printdata2 (t = BT) {
  static FILE * fp1 = fopen ("ini_UW_T2","w");
  foreach()
    fprintf (fp1, "%g %g %g\n", x, u[], t);
  fprintf (fp1, "\n\n");
  fflush (fp1);
}

event integration (i++) {
  // DT = sq(L0/(1 << grid->maxdepth))/3.;    // smaller time step when the tree grid is refined
  double dt = DT;
  dt = dtnext (dt);
  // #if upwind nonconservative
   foreach()
     du[] = -0.5 * (sq(u[]) - sq(u[-1]))/Delta;   // center difference
   foreach()
     u[] += dt*du[];
  // #if Lax_friedrichs
//   foreach()
//     u[] = 0.5 * (u[-1]+u[1]) - dt/(4 * Delta) * (sq(u[1]) - sq(u[-1]));
  // #if lax_wendroff
// foreach(){
//      u[] = u[1] - dt/(4*Delta)*(sq(u[1])-sq(u[-1]))+sq(dt/Delta)/2*((0.25 * //(u[]+u[1])*(sq(u[1])-sq(u[]))) - (0.25*(u[]+u[-1])*(sq(u[])-sq(u[-1]))));
//  }
  boundary ({u});
}

/**
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "upwind non-conservative"
p[][] 'caseUW_T1' u ($1):($2) t 'T0' w l,\
'caseUW_T2' u ($1):($2) t 'T1' w l,
~~~
*/
/**
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "Lax-Friedrichs"
p[][] 'caseLF_T1' u ($1):($2) t 'T0' w l,\
'caseLF_T2' u ($1):($2) t 'T1' w l,
~~~
*/
/**
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "Lax-Wendroff"
p[][] 'caseLW_T1' u ($1):($2) t 'T0' w l,\
'caseLW_T2' u ($1):($2) t 'T1' w l,
~~~
*/

/**
With positive value of initial condition
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "upwind non-conservative"
p[][] 'ini_UW_T1' u ($1):($2) t 'T0' w l,\
'ini_UW_T2' u ($1):($2) t 'T1' w l,
~~~
*/
/**
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "Lax-Friedrichs"
p[][] 'ini_LF_T1' u ($1):($2) t 'T0' w l,\
'ini_LF_T2' u ($1):($2) t 'T1' w l,
~~~
*/
/**
~~~gnuplot
set xlabel "x"
set ylabel "u"
set title "Lax-Wendroff"
p[][] 'ini_LW_T1' u ($1):($2) t 'T0' w l,\
'ini_LW_T2' u ($1):($2) t 'T1' w l,
~~~
*/

