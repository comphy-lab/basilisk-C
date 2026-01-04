/**
# Mass balance test with phase change in a porous medium
This test case verifies the mass conservation in a two-phase flow with phase change.
See [src/shrinking.h](src/shrinking.h) for details on why shrinking models are needed.
We check the error at different grid resolutions and for different porosity update strategies:

- Only shrinking (zeta_policy = 0)
- Constant-volume (zeta_policy = 1)
- Smooth function (zeta_policy = 3)

Analytical solution is possible in this simplified configuration, and used to compute the error.
*/

#define NO_ADVECTION_DIV 1
#define BALANCES_SPHERE // Use analytical expression for initial mass in sphere case in balances.h

#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "balances.h"
#include "view.h"

/**
## Boundary conditions
Top and right boundaries are open. Symmetry boundary conditions are applied at the left and bottom boundaries.
*/
u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
psi[top] = dirichlet(0.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
psi[right] = dirichlet(0.);

int maxlevel = 7, minlevel = 4;
double D0 = 1;
scalar omega[];

scalar fS[];
face vector fsS[];

int main() {
  eps0 = 0.4;
  rho1 = rho2 = 1.;
  mu1 = mu2 = 1.;
  muG = 1.e-3;
  L0 = 1.5*D0;
  DT = 1e-3;

  rhoS = 100.;
  rhoG = 1.;

  int zeta_to_run[3] = {0, 1, 3};
  for (int ii = 0; ii < 3; ii++) {
    zeta_policy = zeta_to_run[ii];
    fprintf(stderr, "Running zeta policy: %d\n", zeta_policy);
     for (maxlevel = 5; maxlevel <=7; maxlevel++) {
      init_grid(1 << maxlevel);
      run();
    }
  }
}

/**
## Initialization event
We initialize a circle of radius D0/2 at the bottom left corner of the domain.
The porosity is initialized based on the initial volume fraction field.
*/

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))
event init(i = 0) {
  fraction (f, circle(x, y, 0.5*D0));

  foreach()
    porosity[] = eps0*f[];

#ifdef BALANCES
  //specify the output file name
  sprintf(mb.name, "balances-%d-%d", maxlevel, zeta_policy);
#endif

  /**
  Output the initial fields for visualization
  */
  char name[80];
  sprintf(name, "OutputFields-init");
  static FILE *fs = fopen(name, "w");
  output_field({u.x, u.y, omega, zeta, porosity, f}, fs);
  fflush(fs);
  fclose(fs);
}

/**
## Chemistry event
We update the porosity field based on a constant reaction rate.
*/
event chemistry(i++) {
  foreach()
    omega[] = rhoS/10.;
  foreach() {
    if (f[] > F_ERR) {
      porosity[] = porosity[]/f[];
      porosity[] += (omega[]*(1.-porosity[])*(1.-zeta[])/rhoS)*dt;
      porosity[] *= f[];
    }
  }
}

// event adapt (i++) {
//     adapt_wavelet_leave_interface({porosity, u.x, u.y}, {f}, 
//       (double[]){1e-2, 1e-2, 1e-2}, maxlevel, minlevel);
// }

event log (t+=0.1) {
  fprintf(stderr, "%g\n", t);
}

/**
# Output event
[balances.h](src/balances.h) prints the mass balances to a file at each time step,
we ensure to call it at the end of the simulation.
At maxlevel 7, we also output the fields for visualization.
*/

event stop (t = 10) {
  if (maxlevel == 7) {
    // Output fields for visualization
    char name[80];
    sprintf(name, "OutputFields-%d", zeta_policy);
    FILE *fs = fopen(name, "w");

    output_field({u.x, u.y, omega, zeta, porosity, f}, fs);
    fflush(fs);
    fclose(fs);
  }
  write_balances();
}

/**
## Post-processing script
~~~gnuplot Error vs maxlevel
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "mass-conservation-error-maxlevel.svg"

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases
array diff_x3[3]  # For SMOOTH cases
analytical = 1.

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$4) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$4) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$4) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

set xlabel "N"
set ylabel "Error"
set format y "10^{%T}"
set key top right box opaque width 0.8
set logscale x 2
set logscale y
#set yr [0.001:1]
set yr [1e-5:0.5]
set xr [2**4:2**8]
set size square
set grid
set title "Total mass conservation"

# Define fitting functions
A0 = 10
B0 = -2
A1 = 10
B1 = -2
A3 = 10
B3 = -2

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)
f4(x) = 0.7*x**(-1)


set print $SHRINK_DATA
do for [i=1:3] { 
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print

set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print

set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print

# Perform curve fitting
fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

# Plot with equations in legend
plot  '$SHRINK_DATA' u 1:2 w p pt 2 ps 2 lc "dark-green" t "Only shrinking",\
      '$SWELLING_DATA' u 1:2 w p pt 8 ps 2 lc "black" t "Constant-volume",\
      '$SMOOTH_DATA' u 1:2 w p pt 4 ps 2 lc "blue" t "Smooth function", \
      f4(x) lw 2 lc "red" title "1^{st} order"
      #f0(x) lw 3 lc "dark-green" title sprintf("fit: %.2f x^{%.2f}", A0, B0), \
      #f1(x) lw 3 lc "black" title sprintf("fit: %.2f x^{%.2f}", A1, B1), \
      #f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3), \

~~~

~~~gnuplot Error in solid mass
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "solid-mass-error-maxlevel.svg"

analytical = exp(-1)

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$2) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$2) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$2) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

A0 = 0.01
B0 = -1
A1 = 0.13
B1 = -2
A3 = 0.01
B3 = -1

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)
f4(x) = 0.4*x**(-1)
f5(x) = 1*x**(-2)

set print $SHRINK_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print
set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print
set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print

fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

set key bottom right box width 1

set xlabel "N"
set ylabel "Error"
set format y "10^{%T}"
set key top right box opaque
set logscale x 2
set logscale y
set yr [1e-5:0.7]
set xr [2**4:2**8]
set size square
set grid
set title "Solid mass conservation"

plot  '$SHRINK_DATA'    u 1:2 w p pt 2  ps 2 lw 3 lc "dark-green" t "Only shrinking",\
      '$SWELLING_DATA'  u 1:2 w p pt 8  ps 2 lw 3 lc "black" t "Constant-volume",\
      '$SMOOTH_DATA'    u 1:2 w p pt 4  ps 2 lw 3 lc "blue" t "Smooth function", \
      f4(x) lw 2 lc "red" title "1^{st} order",\
      f5(x) lw 2 lc "orange" title "2^{nd} order"
      #f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3)
~~~

~~~gnuplot Error in gas mass
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "gas-mass-error-maxlevel.svg"

analytical = 1-exp(-1)

# Initialize arrays to store differences
array diff_x0[3]  # For SHRINK cases
array diff_x1[3]  # For SWELLING cases
array diff_x3[3]  # For SMOOTH cases

# Loop through levels 5 to 7
do for [i=0:2] {
    level = i + 5
    
    # Calculate difference for SHRINK case (x-0)
    stats sprintf("balances-%d-0", level) u (last_y=$3) nooutput
    diff_x0[i+1] = abs(last_y - analytical)
    
    # Calculate difference for SWELLING case (x-1)
    stats sprintf("balances-%d-1", level) u (last_y=$3) nooutput
    diff_x1[i+1] = abs(last_y - analytical)

    # Calculate difference for SMOOTH case (x-3)
    stats sprintf("balances-%d-3", level) u (last_y=$3) nooutput
    diff_x3[i+1] = abs(last_y - analytical)
}

array x_levels[3] = [2**5, 2**6, 2**7]

#set title "Gas mass conservation"
set key top right box opaque

set xlabel "N"
set ylabel "Error"
set format y "10^{%T}"
set logscale x 2
set logscale y
#set yr [0.001:1]
set yr [1e-5:0.7]
set xr [2**4:2**8]
set size square
set grid
set title "Gas mass conservation"

A0 = 0.01
B0 = -1
A1 = 0.13
B1 = -2
A3 = 0.01
B3 = -1

f0(x) = A0*x**(B0)
f1(x) = A1*x**(B1)
f3(x) = A3*x**(B3)
f4(x) = 0.6*x**(-1)

set print $SHRINK_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x0[i])
}
set print
set print $SWELLING_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x1[i])
}
set print
set print $SMOOTH_DATA
do for [i=1:3] {
    print sprintf("%g %g", x_levels[i], diff_x3[i])
}
set print
fit f0(x) '$SHRINK_DATA' using 1:2 via A0, B0
fit f1(x) '$SWELLING_DATA' using 1:2 via A1, B1
fit f3(x) '$SMOOTH_DATA' using 1:2 via A3, B3

plot  '$SHRINK_DATA'   u 1:2 w p pt 2 ps 2 lw 3 lc "dark-green" t "Only shrinking",\
      '$SWELLING_DATA' u 1:2 w p pt 8 ps 2 lw 3 lc "black" t "Constant-volume",\
      '$SMOOTH_DATA'   u 1:2 w p pt 4 ps 2 lw 3 lc "blue" t "Smooth function", \
      f4(x) lw 2 lc "red" title "1^{st} order"
      #f3(x) lw 3 lc "blue" title sprintf("fit: %.2f x^{%.2f}", A3, B3)
~~~

~~~gnuplot Maps for Shrink
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "massb-shrink-map.svg"
set pm3d map interpolate 3,3
load "/root/gnuplot-palettes/jet.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set border lw 2.5
set cbrange [0:0.8]
splot "OutputFields-0" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-0" u 1:2:8 lt 3 lc "black" lw 2, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "white" lw 2
unset multiplot
~~~

~~~gnuplot Maps for Swelling
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "massb-swelling-map.svg"
set pm3d map interpolate 0,0
load "/root/gnuplot-palettes/jet.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:0.8]
splot "OutputFields-1" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-1" u 1:2:8 lt 3 lc "black" lw 2, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "white" lw 2
unset multiplot
~~~

~~~gnuplot Maps for Smooth
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "massb-smooth-map.svg"
set pm3d map interpolate 3,3
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/jet.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
#unset colorbox
set cbrange [0:1]
set cblabel "Poroisty"
splot "OutputFields-3" u 1:2:($7 + (1-$8))

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-3" u 1:2:8 lt 3 lc "black" lw 2, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "white" lw 2
unset multiplot
~~~

~~~gnuplot Maps for initial condition
reset
set terminal svg size 450, 450
#set terminal epslatex size 3.6, 3.6 color colortext
set output "massb-init-map.svg"
set pm3d map interpolate 0,0
#load "/root/gnuplot-palettes/ylgn.pal"
load "/root/gnuplot-palettes/jet.pal"

set multiplot

set xr [0:1]
set yr [0:1]
set size square
unset key
unset xtics 
unset ytics
unset colorbox
set cbrange [0:0.8]
splot "OutputFields-init" u 1:2:7

set contour base
set cntrparam levels discrete 0.5
set cntrparam bspline
set cntrlabel onecolor

unset surface
splot "OutputFields-init" u 1:2:8 lt 3 lc "black" lw 2, "OutputFields-init" u 1:2:8 lt 3 dt 2  lc "white" lw 2
unset multiplot
~~~

~~~gnuplot Plot of porosity

reset
set terminal svg size 350, 350
#set terminal epslatex size 3.6, 3.6 color colortext
set output "massb-porosity.svg"

set xrange [0:1]
set yrange [0:1]
set xlabel "Radial coordinate"
set ylabel "Porosity [-]"
set grid 
set size square
set key bottom right box

#plot "OutputFields-init" u (($2==0 && $8==1 ? $1 : NaN)):(($2==0 && $8==1) ? $7 : NaN) w lp pt 4 t "Initial condition" 

!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-0.025; print last_x, last_y+0.025}}' OutputFields-init > init.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-0.025; print last_x, last_y+0.025}}' OutputFields-0 > shrink.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-0.025; print last_x, last_y+0.025}}' OutputFields-1 > swell.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-0.025; print last_x, last_y+0.025}}' OutputFields-3 > smooth.dat

 plot "swell.dat" u 1:2 w l lc "black" lw 2 t "Constant-volume", \
      "init.dat" u 1:2 w l lc "orange" lw 2 t "Initial condition",\
      "shrink.dat" u 1:2 w l lc "dark-green" lw 2 t "Only shrinking", \
      "smooth.dat" u 1:2 w l lc "blue" lw 2 t "Smooth function"
~~~

~~~gnuplot Plot of porosity
reset
set terminal svg enhanced size 450, 450
#set terminal epslatex size 3.6, 3.6
set output "massb-profiles.svg"


!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-1; print last_x, 1}}' OutputFields-init > init.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-1; print last_x, 1}}' OutputFields-0 > shrink.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-1; print last_x, 1}}' OutputFields-1 > swell.dat
!awk '$2==0 && $8==1 {print $1, $7; if ($7 != 0) {last_x=$1; last_y=$7}} END {if (last_x != "") {print last_x, last_y-1; print last_x, 1}}' OutputFields-3 > smooth.dat

unset key
set grid
set xrange [0:0.6]
set xtics 0.2
set yrange [0:1]
set ytics 0.5
unset key
set label "Porosity" at screen 0.02, 0.5 rotate by 90 center
set label "Radial coordinate" at screen 0.55, 0.1 center

set xtics format ""
set multiplot layout 3,1 margins 0.15, 0.95, 0.2, 0.95 spacing 0,0.05

set arrow 1 from 0.22, 0.23 to 0.22, 0.9
set label "t" at 0.23, 0.54
set label "$\\Gamma _{0} = \\Gamma _{f}$" at 0.35, 0.3
set label "$\\epsilon _{f}$" at 0.3, 0.91
set label "$\\epsilon _{0}$" at 0.3, 0.55
set title "Constant-volume" offset 0,-1
plot  "init.dat" u 1:2 w l lw 6 lc "black" notitle,\
      "swell.dat" u 1:2 w l lw 6 lc "orange" t "Constant-volume"

unset label
unset arrow 1
set arrow 2 from 0.42,0.5 to 0.22,0.5
set label "t" at 0.32, 0.6
set title "Only shrinking"
set label "$\\Gamma _{0}$" at 0.5, 0.5
set label "$\\Gamma _{f}$" at 0.24, 0.25
set label "$\\epsilon _{0} = \\epsilon _{f}$" at 0.35, 0.3
plot  "init.dat" u 1:2 w l lw 6 lc "black" notitle,\
      "shrink.dat" u 1:2 w l lw 6 lc "dark-green" t "Only shrinking"

unset label
unset arrow 2
set xtics format "%g"
set arrow 3 from 0.32,0.27 to 0.22, 0.9
set label "t" at 0.29, 0.56
set title "Smooth function"
set label "$\\Gamma _{0}$" at 0.5, 0.5
set label "$\\Gamma _{f}$" at 0.42, 0.7
set label "$\\epsilon _{0}$" at 0.35, 0.27
set label "$\\epsilon _{f}$" at 0.3, 0.91
plot  "init.dat" u 1:2 w l lw 6 lc "black" notitle,\
      "smooth.dat" u 1:2 w l lw 6 lc "blue" t "Smooth function"
unset multiplot
~~~
*/