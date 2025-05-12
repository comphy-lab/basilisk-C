/**
# Vof tracer and particle tracers.

The case is based on a [test case](/src/test/revered/c) for the `vof`
advection. Results are visual:

Particle Legend:  

Green: Inside drop  
Magenta: Outside drop  
Red: On interface  
Blue: Unbounded  

![Swirling tracers](reversedp/mov.mp4)

![The Magenta circle is the exact solution](reversedp/result.png)
 */

#include "advection.h" //diffusive tracer field
#include "vof.h"       //Interface tracer
#include "view.h"
#include "vof-tracer-particles.h"
#include "scatter2.h"

scalar f[], * interfaces = {f}, *tracers = NULL;
Particles Pin, Pout, Pon, Pun;
double T = 15;

int main () {
  L0 = 1;
  X0 = Y0 = -L0/2;
  DT = 0.1;
  N = 64;
  run();
}

double R = 0.2, xo = -0.2, yo = -.236338;
#define circle(x,y) (sq(R) - (sq(x - xo) + sq(y - yo)))

event init (t = 0) {
  fraction (f, circle(x,y)); 
  Pin  = new_vof_tracer_particles (10, 1);  // f[] == 1 phase
  Pout = new_vof_tracer_particles (10, 0);  // f[] == 0 phase
  Pon  = new_vof_tracer_particles (10, 3);  // On the interface
  Pun  = new_vof_tracer_particles (10, -1); // Unbounded
  foreach_particle() {
    p().x = xo + R*sin(2.5*(double)_j_particle);
    p().y = yo + R*cos(2.5*(double)_j_particle);
  }
}

event velocity (i++) {
  vertex scalar psi[];
  foreach_vertex()
    psi[] = - 1.5*sin(2.*pi*t/T)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;
  trash ({u});
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    uf.x[] = f.x*(psi[0,1] - psi[])/Delta;
  boundary ((scalar *){uf});
}

event movie (t += 0.1) {
  draw_vof("f");
  foreach_P_in_list (tracer_particles) { 
    scatter (P, pc = {sin(P), cos(P), sin(2.4*P)});
  }
  save ("mov.mp4");
}

event stop (t = T) {
  scalar b[];
  fraction (b, circle(x,y));
  draw_vof("b", lc = {1,0,1}, lw = 2);
  draw_vof("f");
  translate (z = 0.01) {
    foreach_P_in_list (tracer_particles) { 
      scatter (P, pc = {sin(P), cos(P), sin(2.4*P)});
    }
  }
  save ("result.png");
  system ("ffmpeg -i mov.mp4 mov_1.mp4");
}
/**
<div class="figure">
<video controls="" preload="metadata" width="600px">
<source src="http://www.basilisk.fr/sandbox/Antoonvh/reversedp/mov.mp4" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Use HTML
</p>
</div>

<div class="figure">
<video controls="" preload="metadata" width="600px">
<source src="http://www.basilisk.fr/sandbox/Antoonvh/reversedp/mov_1.mp4" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Use HTML and converted movie
</p>
</div>
*/