/**
![Mixing of viscous fluids may be a challange. Image via [hub pages](https://hubpages.com/games-hobbies/Mixing-Colors-To-Painting-Miniatures-And-Models-A-How-to-Guide)](https://usercontent2.hubstatic.com/4605727.jpg)

# Laminar chaotic advection

On this page we demonstrate how laminar advection can be chaotic. As
such it may be employed to achieve mixing with a Stokes flow. The
setup mimics an `industrial' setup where two rotors are used to mix
passive tracer particles.
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#define BVIEW 1
#include "particles.h"

face vector muc[];
int maxlevel = 8;

int main() {
  init_grid (1 << maxlevel);
  L0 = 20. + pi;
  stokes = true;
  mu = muc;
  X0 = Y0 = -L0/2;
  run();
}
/**
The rotation alternatingly switches between the left and right
cylinder with a time interval $\Delta T = \pi$. 
 */
#define LEFT_CYL (U*(cos(t) > 0)*cos(t))
#define RIGHT_CYL (U*(cos(t) < 0)*cos(t))

double U = -7.5;
double a1 = 2;
#define RAD(a) (pow(sq(x - a) + sq(y), 0.5))
u.n[embed] = dirichlet (x > 0? -y*RIGHT_CYL : -y*LEFT_CYL);
u.t[embed] = dirichlet (x > 0? (x - a1)*RIGHT_CYL : (x + a1)*LEFT_CYL);
/**
Chaos is introduced to the system by randomizing the rotation velocity
(`U`) of the cylinders between each turn.
 */
event randomize_U(t = pi/2.; t += pi)
  U = (5. + 2.*noise())*sign(noise());

/**
We set the viscosity for this Stokes-flow problem, taking into
account the cylinders' boundaries.
 */
event properties (i++){
  foreach_face()
    muc.x[] = 5.*fm.x[];
  boundary ((scalar*){muc});
}

event init (t = 0){
  init_particles_2D_square_grid (20, 0, 0, 1);
  srand (time(NULL));
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (RAD(a1) - 1) * (x > 0) + (x < 0)*(RAD(-a1) - 1);
  fractions (phi, cs, fs);
  NITERMIN = 2;
  TOLERANCE = 1E-4;
  DT = 0.015;
}

event adapt (i++)
  adapt_wavelet ((scalar*){cs, u}, (double[]){0.0005, 0.01, 0.01}, maxlevel);
/**
## Results

We can view how the particles are redistributed over space as time progresses:

![Mixing seems to be achieved to some degree](laminarmixing/mixing.mp4)

*/

event movie (t += 0.05; t <= 100){
  view (fov = 12);
  draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.5,0.5});
  scatter (loc, s = 10);
  scalar fl[], fr[];
  foreach(){
    fl[] = (1 - cs[])*(x < 0);
    fr[] = (1 - cs[])*(x > 0);
  }
  float dc[3] = {U/7*(U > 0), 0 , -U/7*(U < 0)};
  if (fmod(t + pi/2, 2*pi) < pi){
    draw_vof ("fl", lw = 5, lc = {dc[0], dc[1], dc[2]});
    draw_vof ("fr", lw = 5);
  
  }else{
    draw_vof ("fr", lw = 5, lc = {dc[0], dc[1], dc[2]});
    draw_vof ("fl", lw = 5);
  }
  save ("mixing.mp4");
}
/**
Furthermore, we are interested to see the *chaotic* motion of the
tracers. Therefore, we follow the paths of two particles that are
initially placed close to each other at the bottom left corner of the
initial square (see movie).
 */

event trace_paths(i+=5){
  static FILE * fpt = fopen ("paths", "w");
  fprintf (fpt, "%g\t%g\t%g\t%g\t%g\t%g\n",
           loc[0].x, loc[0].y, loc[1].x, loc[1].y, t,
           sqrt(sq(loc[0].x - loc[1].x) + sq(loc[0].y - loc[1].y)));
}  

/**
The resulting paths can be plotted:

~~~gnuplot The paths diverge 
set xr [-5 : 5]
set yr [-3 : 3]
set xlabel 'x'
set ylabel 'y'
set size ratio -1
plot 'paths' u 1:2 w l , 'paths' u 3:4 w l
~~~

Eventough the initial conditions for the particles are *almost*
identical, their paths in the $x-y$ parameter space are very
different. Well done chaotic laminar advection!

~~~gnuplot Their distance over time
set xr [0 : 100]
set yr [0 : 5]
set xlabel 'time'
set ylabel 'distance'
set size ratio 0.5
set key off
plot 'paths' u 5:6 w l
~~~
 */
