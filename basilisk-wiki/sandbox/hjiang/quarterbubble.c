/**
# Axisymmetric quarter-drop benchmark for surfactant transport

This test case is used to validate the sharp-interface surfactant transport
algorithm in an axisymmetric configuration. A quarter circular drop is
initialized in the lower-left corner of the domain, and a uniform initial
surfactant concentration is prescribed on the interface.

The imposed velocity field corresponds to a purely expansive flow,

$\mathbf{u} = K(x-x_c,\, y-y_c)$,

so that the drop expands self-similarly in time. In this case, the
interfacial area increases continuously, and the surfactant concentration
decreases accordingly. The exact solution is

$\Gamma(t) = \exp(-2t)$

for the axisymmetric case.

This benchmark is used to assess the accuracy of the numerical method by
comparing the spatially averaged and extreme values of the surfactant
concentration against the exact solution at different grid resolutions.
*/

#define SURFACE 1
#define F_ERR 1.e-6
#define S_ERR 1.e-6

#include "axi.h"
#include "navier-stokes/centered.h"
#include "surfactant/vof-surfactant.h"
#include "tension.h"
#include "surfactant/adapt_wavelet_leave_interface.h"
#include "view.h"

/**
### Boundary Conditions

Open (outflow) boundary conditions are imposed on every boundary, since the
droplet is initialized in the center of the domain. The same boundary conditions
are copied in other velocity and pressure fields used by the phase change model.
*/

double density = 1.0;

scalar c[], *interfaces = {c};

double D0 = 1.e-3;
int maxlevel = 9;
int minlevel = 4;

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

/**
### Model Data

Additional data useful for the simulation are the maximum level of refinement
and the initial droplet radius.
*/

int main()
{
  L0 = 3. * D0;
  origin(0., 0.);
  c.sigma = 1.e-5;
  for (maxlevel = 5; maxlevel <= 7; maxlevel++) {
    init_grid(1 << maxlevel);
    // DT = 2.e-4 / (1 << (maxlevel - 5));
    run();
  }
}

event init(t = 0)
{
  fraction(c, sq(D0 / 2.) - pow(x, 2) - pow(y, 2));
  foreach () {
    surfactant[] = (c[] > F_ERR && c[] < 1. - F_ERR) ? density : 0.0;
  }
  boundary({surfactant});
  // foreach () {
  //   A[] = (c[] > F_ERR && c[] < 1. - F_ERR) ? density : 0.0;
  // }
  //boundary({A});
}

/**
## Prescribed expansive velocity field

A linear expansive velocity field is imposed so that the interface expands
self-similarly in time. 
*/

double K = 1.;
double xc = 0., yc = 0.;
event stability(i++)
{
  trash({u, uf});

  foreach_face(x) uf.x[] = K * (x - xc) * fm.x[];
  foreach_face(y) uf.y[] = K * (y - yc) * fm.y[];

  foreach () {
    u.x[] = K * (x - xc);
    u.y[] = K * (y - yc);
  }

  boundary((scalar *){u});
  boundary((scalar *){uf});
}

#if TREE
event adapt(i++) { adapt_wavelet_leave_interface({surfactant}, {c}, (double[]){1.e-3}, maxlevel, minlevel = 2, 1); }
// event adapt(i++) { adapt_wavelet_leave_interface({A}, {c}, (double[]){1.e-3}, maxlevel, minlevel = 2, 1); }
//
#endif

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/
#if AXI
double exact(double time) { return exp(-2 * time); }
#else
double exact(double time) { return exp(-time); }
#endif

event output(t += 0.1)
// event output(i++)
{
  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);

  double surfactantsum = 0.;
  double surfactantMin = HUGE;
  double surfactantMax = -HUGE;
  double Atot = 0.;

  foreach () {
    if (c[] > F_ERR && c[] < 1. - F_ERR) {
      coord n = interface_normal(point, c), p;
      double alpha = plane_alpha(c[], n);
      double A = plane_area_center(n, alpha, &p) * Delta;
#if AXI
      A *= y + p.y * Delta;
#endif
      if (surfactant[] > 0.) {
        surfactantsum += surfactant[] * A;
        Atot += A;
      }
      if (surfactant[] > 0. && surfactant[] < surfactantMin) surfactantMin = surfactant[];
      if (surfactant[] > 0. && surfactant[] > surfactantMax) surfactantMax = surfactant[];
    }
  }

  double surfactantAvg = surfactantsum / Atot;
  double L1 = fabs(surfactantAvg - exact(t)) / exact(t);
  double L8 = max(fabs(surfactantMax - exact(t)) / exact(t), fabs(surfactantMin - exact(t)) / exact(t));

  static FILE *fp = fopen(name, "w");
  fprintf(fp, "%.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", t, exact(t), surfactantAvg, L1, surfactantMax,
    surfactantMin, L8);
  fflush(fp);
}

// #if 0
event movie(t += 0.01)
{
  clear();
  draw_vof("c", lw = 1.5);
  squares("surfactant", min = 0., max = statsf(surfactant).max, linear = false, cbar = true, pos = {-0.75, -0.75},
    lc = {1.0, 1.0, 1.0});
  box();
  save("surfactant.mp4");
}
// #endif

event end(t = 0.7);

/** 

~~~gnuplot

reset

stats "OutputData-5" using 4 nooutput name "L1_LEVEL5"
stats "OutputData-6" using 4 nooutput name "L1_LEVEL6"
stats "OutputData-7" using 4 nooutput name "L1_LEVEL7"
##stats "OutputData-8" using 4 nooutput name "L1_LEVEL8"
##stats "OutputData-9" using 4 nooutput name "L1_LEVEL9"
##stats "OutputData-10" using 4 nooutput name "L1_LEVEL10"

stats "OutputData-5" using 7 nooutput name "L8_LEVEL5"
stats "OutputData-6" using 7 nooutput name "L8_LEVEL6"
stats "OutputData-7" using 7 nooutput name "L8_LEVEL7"
##stats "OutputData-8" using 7 nooutput name "L8_LEVEL8"
##stats "OutputData-9" using 7 nooutput name "L8_LEVEL9"
##stats "OutputData-10" using 7 nooutput name "L8_LEVEL10"

set print "Errors_L1_L8.csv"

print sprintf ("%d %.12f %.12f", 2**5, L1_LEVEL5_mean, L8_LEVEL5_mean)
print sprintf ("%d %.12f %.12f", 2**6, L1_LEVEL6_mean, L8_LEVEL6_mean)
print sprintf ("%d %.12f %.12f", 2**7, L1_LEVEL7_mean, L8_LEVEL7_mean)
##print sprintf ("%d %.12f %.12f", 2**8, L1_LEVEL8_mean, L8_LEVEL8_mean)
##print sprintf ("%d %.12f %.12f", 2**9, L1_LEVEL9_mean, L8_LEVEL9_mean)
##print sprintf ("%d %.12f %.12f", 2**10, L1_LEVEL10_mean, L8_LEVEL10_mean)

unset print

reset
set xlabel "N" font ", 28" offset 0,-2
set ylabel "Error" font ", 28" offset -5,0

set xtics font ",22"
set ytics font ",22"
set key font ",22"

set border lw 2
set tics lw 2
set grid lw 1.5

set logscale x 2
set logscale y
set format y "10^{%L}"
#set xtics ("2^{4}" 2**4, "2^{5}" 2**5, "2^{6}" 2**6, "2^{7}" 2**7, "2^{8}" 2**8)

set xr[2**4:2**8]
set yr[10**-6:10**-2]
set size ratio 0.75
set grid

set tmargin 3
set bmargin 8
set lmargin 8
set rmargin 3

f(x) = a*x**-b
fit f(x) "Errors_L1_L8.csv" u 1:2 via a,b

f1(x) = a1*x**-b1
fit f1(x) "Errors_L1_L8.csv" u 1:3 via a1,b1

ftitle(a,b) = sprintf("%.3f*x^{%4.2f}", a, -b)

plot "Errors_L1_L8.csv" u 1:2 w p ps 3 title "L^1", \
     "Errors_L1_L8.csv" u 1:3 w p ps 3 title "L^{inf}", \
     f(x) w l lw 3 title ftitle(a, b), \
     f1(x) w l lw 3 title ftitle(a1, b1)

~~~
*/
