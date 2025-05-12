/**
![Confetti is fun. Photo by [Ben
 Aveling](https://commons.wikimedia.org/wiki/User_talk:BenAveling) via
 [wikimedia](https://commons.wikimedia.org/wiki/File:IMG_2566_Confetti.jpg)
 ](https://upload.wikimedia.org/wikipedia/commons/thumb/5/50/IMG_2566_Confetti.jpg/1280px-IMG_2566_Confetti.jpg){width="400px"}

# A confetti cannon

![New conteffi particles are added at fixed intervals](confetti/movp.mp4)

~~~gnuplot Averaged location of first (pink) confetti drop
set size ratio -1
set key right outside
plot 'out' u 3:4:2 w l palette lw 2 t 'stddev'
~~~

## Set-up

Note that we do not (have to) have names for the added particles.
 */
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

#define POSX (R1*sin(t/6))
#define POSY (R1*cos(t/6))

double R1 = 1, R2 = 0.2;

Particles initial;

int main() {
  const face vector muc[] = {1./500., 1./500.};
  mu = muc;
  L0 = 10.;
  X0 = Y0 = -L0/2.;
  DT = 0.1;
  run();
}

event init (t = 0)
  initial = init_tp_circle();

event forcing (i++) {
  double angle = t/2.;
  double U = 1;
  foreach()
    if (sq(x - POSX) + sq(y - POSY) < sq(R2)) {
      u.x[] = U*sin(angle);
      u.y[] = U*cos(angle);
    }
}

event add (t += 2)
  init_tp_square (xm = POSX, ym = POSY, l = R2);

event adapt (i++)
  adapt_wavelet ({u.x, u.y}, (double[]){0.1, 0.1}, 8, 6);

event mov (t += 0.3) {
  foreach_P_in_list(tracer_particles) {
    scatter (P, pc = {sin(P), cos(P), sin(P*2.4)});
  }
  box();
  save ("movp.mp4");
  
  pstats ps = statsp (tracer_particles[1]);
  if (pid() == 0) 
    printf ("%g %g %g %g\n",t,
	    sqrt(sq(ps.stddev.x) + sq(ps.stddev.y)), //length of stddev vector
	    ps.avg.x, ps.avg.y);
}

event stop (t = 200);
