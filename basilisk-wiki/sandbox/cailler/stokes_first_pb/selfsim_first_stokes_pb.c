/**
# Stokes First Problem with the dedicated self-similar solver

We want to converge spontaneously towards the self-similar solution of the
[Stokes first problem](https://en.wikipedia.org/wiki/Rayleigh_problem), 
with a dedicated 2D-incompressible centered *Navier--Stokes* solver, 
to compare with [the results obtained in the physical space](/sandbox/cailler/stokes_first_pb/first_stokes_pb.c) 
(see this link for a complete presentation of the problem).

<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
**The self-similar solver has been explained extensively on other problems:**

  + ***[here](/sandbox/cailler/self_sim_DNS/README)*** (theory);
  + ***[or here](/sandbox/cailler/lamb_oseen/selfsim_lamb_conv.c)*** 
  (summary).
</div>


![Diffusion of the velocity in the self-similar domain](selfsim_first_stokes_pb/fstokes_selfsim_space_N7.mp4)(width="500" height="500")

We recall the self-similar (analytical) solution:
$$
  u(y, t) = U_0 \, \overline{u}(\eta)  
          = U_0 \, \mathrm{erfc} \left( \dfrac{\eta}{2} \right), 
  \quad 
  \eta = \dfrac{y}{\sqrt{\nu t}} 
$$

In practice, to solve in the self-similar space with `Basilisk`, we introduce 
a *logarithmic time* $\tau = \ln t$. Therefore, the diffusion equation in the 
physical space 
$$
\partial_t u = \nu \, \partial_{yy} u
$$
now becomes in the self-similar coordinates:
$$
\partial_\tau \overline{u} 
- \dfrac{1}{2} \eta \, \partial_\eta \overline{u} 
=
\partial_{\eta \eta} \overline{u} 
$$
with $u(y, t) = U_0 \, \overline{u}(\eta, \tau)$. 
This last equation is the one solved by the self-similar solver.

## Code

### General Parameters

The simulation is set in a squared box of size $[0 \, ; 10] \times [0 \, ; 10]$, 
as the analytical solution vanishes for $\eta \geqslant 5$.
*/


#define LEVEL 7
#define MAXLEVEL 8
#define SIZE 10.
#define T_END 1.e2
#define U_0 1.
#define MOVIE 1

#include "grid/multigrid.h"
#include "selfsim_centered_stokes.h"

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
  const face vector muc[] = {1.e-0, 1.e-0};
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
event profiles (t = {0., 0.001, 0.01, 0.1, 0.5, 1., 10., 20., 50., 100.}){
  char filename[80];
  sprintf(filename, "selfsim_ux_prof_%g.dat", t);
  FILE * fp = fopen(filename, "w");
  for (double y = 0. ; y <= SIZE ; y += 0.001)
    fprintf (fp, "%g %g\n", y, fabs(interpolate (u.x, .8*SIZE, y)));
  fclose (fp);
}

/**
We can now plot the results, showing the spontaneous convergence towards 
the steady solution with the self-similar solver:

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

set term pngcairo enhanced size 500,500
set output 'selfsim_first_stokes_pb.png'

set size ratio 1

# Plot various velocity profiles in the self-similar space
set title "Time evolution in the self-similar space [SELF-SIM. SOLVER]"
set xlabel "η = y/√νt"  
set ylabel "u_x/U"
set xrange [0:5]
set yrange [0:1]

f(x) = erfc(x/2.)

p "selfsim_ux_prof_0.001.dat" u 1:2 w l ls 1 t "τ = 0.001", \
  "selfsim_ux_prof_0.01.dat" u 1:2 w l ls 2 t "τ = 0.01", \
  "selfsim_ux_prof_0.1.dat" u 1:2 w l ls 3 t "τ = 0.1", \
  "selfsim_ux_prof_0.5.dat" u 1:2 w l ls 4 t "τ = 0.5", \
  "selfsim_ux_prof_1.dat" u 1:2 w l ls 5 t "τ = 1", \
  "selfsim_ux_prof_10.dat" u 1:2 w l ls 6 t "τ = 10", \
  "selfsim_ux_prof_20.dat" u 1:2 w l ls 7 t "τ = 20", \
  "selfsim_ux_prof_50.dat" u 1:2 w l ls 8 t "τ = 50", \
  "selfsim_ux_prof_100.dat" u 1:2 w l ls 9 t "τ = 100", \
  f(x) dt 2 lc rgb "red" lw 2 t "Analytical solution"
~~~

Some discrepancies can be seen between the exact solution and the one computed 
with the self-similar solver. This is *a priori* due to a wrong behaviour of the 
*left* outflow boundary condition (see the movie at the start of this page), 
modifying slightly the obtained steady solution. 
A countermeasure would be to increase the simulation box and take transverse 
profiles of the longitudinal velocity far away from the left boundary, to not 
be impacted/reduce the impact on the results.


### Movie
*/

#if MOVIE 

#include "view.h"
#include "cmaps.h"

event mov_p_u (t+=1e-1) {
  char legend[100];
  sprintf(legend, "tau = %0.2g", t);
  view (tx = -0.5, ty = -0.4, width = 1000., height = 1000.); // for t_end=1e-2
  box();
  squares("u_L2", spread=-1, linear=true, map=viridis);
  vectors("u", scale=0.01, lc = {1, 1, 1});
  draw_string(legend, 2, size = 20., lc = {.851, 0., 0.043}, lw = 2.1);
  save ("fstokes_selfsim_space_N7.mp4");
}

#endif 
