/**
# Couette flow

We test here the ability of the multilayer solver to reproduce a Couette flow.
*/ 

#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/implicit.h"

#define T_END 5.

/**
A flat film with homogeneous surface tension and viscosity is initialised.
A constant acceleration is imposed at the top.
*/ 
int main()
{
  periodic(right);
  
  nu = 1. [2,-1];
  CFL_H = 0.5;

  nl = 8;
  for (N = 1; N <= 64; N*=2)  
    run();

  N = 16;  
  for (nl = 1; nl <= 16; nl*=2)
    run();
}

event init (i = 0)
{
  foreach()
    foreach_layer()
      h[] = 1.[1]/nl;
  boundary({h});
}

const vector du[] = {1.,0.};
event viscous_term (i++)
{
  dut = du;
}

event output (t += T_END/10.; t <= T_END)
{
  char name[80];
  sprintf (name, "velocity-nl%d-N%d", nl, N);  
  static FILE * fp = fopen (name, "w");
  int J = 0;
  foreach() {
    double z = zb[];
    if (J==N/2)
      foreach_layer() {
        z += h[]/2;
        fprintf (fp, "%g %g %g\n", z, u.x[], x);
        z += h[]/2;
      }
    J += 1;
  }
  fprintf (fp, "\n");
}

event logfile (t = end) {
  fprintf(stderr, "%d %d %g\n", N, nl, fabs(statsf(u.x).sum*nl/0.5-1.));
}


/**
~~~gnuplot a
set terminal @SVG enhanced size 640,640 font ",14"
set output 'velocity.svg'
set xlabel "velocity"
set ylabel "height"
plot 'velocity-nl8-N64' u ($2):1 w l t 'nl = 16'
~~~

~~~gnuplot Final error with horizontal resolution
set output 'errorN.svg'
set xlabel "N"
set ylabel "relative error"
set logscale x
set logscale y
plot 'log' u 1:3 w p
~~~

~~~gnuplot Final error with horizontal resolution
set output 'errornl.svg'
set xlabel "nl"
plot 'log' u 2:3 w p
~~~

Accurancy is nearly independant of N and nl but stays surprisingly quite high.
This could be related to the fact that the BC is on acceleration and not velocity.
*/