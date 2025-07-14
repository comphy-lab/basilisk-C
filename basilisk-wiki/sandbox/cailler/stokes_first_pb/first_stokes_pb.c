/**
# Stokes First Problem

We want to solve the 
[Stokes first problem](https://en.wikipedia.org/wiki/Rayleigh_problem), 
with the 2D-incompressible centered *Navier--Stokes* solver. 
In this problem, a semi-infinite fluid domain of kinematic viscosity 
$\nu := \mu / \rho$, initially at rest,
is bounded at its bottom by a supposedly infinite flat plate, the latter 
being then moved at constant speed $U_0$ in the $x-$direction. 
It creates by diffusion a **rectilinear** flow due to the plate 
motion^[Hence, if $\mathbf{u} := (u, v)^T$, then $v = 0$.], 
without any pressure gradient in both directions.

This problem can be thus formalized according to the simple following set of 
equations:

  1. **Continuity Equation**

$$
\partial_x u = 0 \Rightarrow \mathbf{u} = u(y, t) \, \mathbf{e}_x
$$

  2. **Fluid domain (momentum balance)**

$$
\partial_t u = \nu \, \partial_{yy} u
$$

  3. **Initial and Boundary Conditions**

$$
u(y, 0) = 0, \quad u(0, t) = U_0, \quad u(+\infty, t) = 0
$$

![Diffusion of the velocity in the domain](first_stokes_pb/fstokes_phys_space_N7.mp4)(width="500" height="500")


## Self-similar solution

From the main equation of diffusion, one can find that we do have the scalings:
$$
\dfrac{U_0}{t} \sim \nu \, \dfrac{U_0}{y^2} \,\, \Rightarrow \,\, y \sim \sqrt{\nu t}
$$
We can then therefore introduce the **self-similar** variable $\eta$ and 
function $\overline{u}$:
$$
\eta(y,t) := \dfrac{y}{\sqrt{\nu t}}, \quad 
u(y, t) = U_0 \, \overline{u} \left( \eta(y, t) \right)
$$
which can be plugged into the diffusion equation to get a simple ODE:
$$
\dfrac{\mathrm{d}^2 \overline{u}}{\mathrm{d} \eta^2} 
+ \dfrac{\eta}{2} \, \dfrac{\mathrm{d} \overline{u}}{\mathrm{d} \eta}
=
0 
$$
Along with the BCs, a simple analytical solution is obtained:
$$
\boxed{
  u(y, t) = U_0 \, \mathrm{erfc} \left( \dfrac{y}{2 \, \sqrt{\nu t}} \right) 
          = U_0 \, \mathrm{erfc} \left( \dfrac{\eta}{2} \right) 
}
$$

We will show next how to recover the self-similar solution with `Basilisk`.


## Code

### General Parameters

*/

#define LEVEL 7
#define MAXLEVEL 8
#define SIZE 1.
#define T_END 1.e2
#define U_0 1.
#define MOVIE 1

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

scalar omega[]; // vorticity field
scalar u_L2[]; // for movie

/** 
### Boundary Conditions 
*/

// Plate moving at constant speed in the positive 'x' direction:
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(U_0);

// No flow far from the bottom plate:
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

// Outflow condition on the right:
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);

// Outflow condition on the left:
p[left] = dirichlet(0.);
pf[left] = dirichlet(0.);
u.n[left] = neumann(0.);

/** 
### Generic Events
*/

int main() {
  size(SIZE);
  init_grid(1 << LEVEL);
  // viscosity
  const face vector muc[] = {1.e-3, 1.e-3};
  mu = muc;

  run();
}

event init (t = 0){
  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

event vorti(i++){
  vorticity(u, omega);
  foreach()
    u_L2[] = sqrt(sq(u.x[]) + sq(u.y[]));
}

// ADAPTIVE REFINEMENT
#if TREE
event adapt (i++){
  adapt_wavelet ((scalar *){u}, (double []){5.e-4, 5.e-4}, 
                  maxlevel = MAXLEVEL, minlevel = LEVEL - 1);
}
#endif 

/** 
### Outputs 
*/

event monitoring (i++){
  printf("time = %g ; iteration #%i\n", t, i);
  fflush(stdout); // empty buffy
}

event end (t = T_END){}

/** 
We take some snapshots of the longitudinal velocity profile in the transverse 
direction at various times:
*/
event profiles (t = {0., 0.01, 0.05, 0.1, 0.5, 1., 5., 10., 50., 100.}){
  char filename[80];
  sprintf(filename, "ux_prof_%g.dat", t);
  FILE * fp = fopen(filename, "w");
  for (double y = 0. ; y <= SIZE ; y += 0.001)
    fprintf (fp, "%g %g\n", y, fabs(interpolate (u.x, SIZE/2., y)));
  fclose (fp);
}

/**
We can now plot the results, firstly in the physical space, and then 
by rescaling them according to the self-similar coordinates to show the 
scale invariance of this problem:

~~~gnuplot Self-similarity of the axial velocity
reset session

set style line 1  lc rgb '#0025ad' lt 1 lw 2.5 # --- blue
set style line 2  lc rgb '#0042ad' lt 1 lw 2.5 #      .
set style line 3  lc rgb '#0060ad' lt 1 lw 2.5 #      .
set style line 4  lc rgb '#007cad' lt 1 lw 2.5 #      .
set style line 5  lc rgb '#0099ad' lt 1 lw 2.5 #      .
set style line 6  lc rgb '#00ada4' lt 1 lw 2.5 #      .
set style line 7  lc rgb '#00ad88' lt 1 lw 2.5 #      .
set style line 8  lc rgb '#00ad6b' lt 1 lw 2.5 #      .
set style line 9  lc rgb '#00ad4e' lt 1 lw 2.5 #      .
set style line 10 lc rgb '#00ad31' lt 1 lw 2.5 #      .
set style line 11 lc rgb '#00ad14' lt 1 lw 2.5 #      .
set style line 12 lc rgb '#09ad00' lt 1 lw 2.5 # --- green

set term pngcairo enhanced size 1000,500
set output 'first_stokes_pb.png'

set size ratio 1
set multiplot layout 1,2

# Plot various velocity profiles in the physical space
set title "Time evolution in the physical space"
set xlabel "y"  
set ylabel "u_x"
set xrange [0:1]
set yrange [0:1]

p "ux_prof_0.01.dat" u 1:2 w l ls 1 t "t = 0.01", \
  "ux_prof_0.05.dat" u 1:2 w l ls 2 t "t = 0.05", \
  "ux_prof_0.1.dat" u 1:2 w l ls 3 t "t = 0.1", \
  "ux_prof_0.5.dat" u 1:2 w l ls 4 t "t = 0.5", \
  "ux_prof_1.dat" u 1:2 w l ls 6 t "t = 1", \
  "ux_prof_10.dat" u 1:2 w l ls 8 t "t = 10", \
  "ux_prof_50.dat" u 1:2 w l ls 10 t "t = 50", \
  "ux_prof_100.dat" u 1:2 w l ls 12 t "t = 100"


# Plot various velocity profiles in the self-similar space
set title "Time evolution in the self-similar space"
set xlabel "η = y/√νt"  
set ylabel "u_x/U"
set xrange [0:5]
set yrange [0:1]

f(x) = erfc(x/2.)

p "ux_prof_0.01.dat" u ($1/sqrt(0.001*0.01)):2 w l ls 1 t "t = 0.01", \
  "ux_prof_0.05.dat" u ($1/sqrt(0.001*0.05)):2 w l ls 2 t "t = 0.05", \
  "ux_prof_0.1.dat" u ($1/sqrt(0.001*0.1)):2 w l ls 3 t "t = 0.1", \
  "ux_prof_0.5.dat" u ($1/sqrt(0.001*0.5)):2 w l ls 4 t "t = 0.5", \
  "ux_prof_1.dat" u ($1/sqrt(0.001*1)):2 w l ls 6 t "t = 1", \
  "ux_prof_10.dat" u ($1/sqrt(0.001*10)):2 w l ls 8 t "t = 10", \
  "ux_prof_50.dat" u ($1/sqrt(0.001*50)):2 w l ls 10 t "t = 50", \
  "ux_prof_100.dat" u ($1/sqrt(0.001*100)):2 w l ls 12 t "t = 100", \
  f(x) dt 2 lc rgb "red" lw 2 t "Analytical solution"

unset multiplot
~~~
*/

/**
### Movie

We can make a movie:
*/

#if MOVIE 

#include "view.h"
#include "cmaps.h"

event mov_p_u (t+=1e-1) {
  char legend[100];
  sprintf(legend, "t = %0.2g", t);
  view (tx = -0.5, ty = -0.4, width = 1000., height = 1000.); 
  box();
  squares("u_L2", spread=-1, linear=true, map=viridis);
  vectors("u", scale=0.001, lc = {1, 1, 1});
  draw_string(legend, 2, size = 20., lc = {.851, 0., 0.043}, lw = 2.1);
  save ("fstokes_phys_space_N7.mp4");
}

#endif