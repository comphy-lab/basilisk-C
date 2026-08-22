/**
[<img src="https://img.youtube.com/vi/sYdXa_yp0fc/sddefault.jpg">](https://www.youtube.com/watch?v=sYdXa_yp0fc)  


## A string in an instrument

Inspired by [this movie](https://www.youtube.com/watch?v=sYdXa_yp0fc)
by *Minute physics*, we model the deformation of string in various instruments 
($h(x,t)$) using their equation of motion, $$\frac{\partial^2 h}{\partial
t^2} = k \frac{\partial h^2}{\partial x^2},$$ With $k$ some stiffness
and mass related parameter.

### Plucked guitar string

![As seen on TV!](string/mov0.mp4)(loop)

### Hammer-struck piano string

![As seen on TV!](string/mov1.mp4)(loop)

### Violin string

![Not quite as seen on TV! Compared to the other two, we can see the characteristic "attack" of a violin, i.e. a smooth and growing bluid-up of the note](string/mov2.mp4)(loop)
*/
#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"

scalar h[], h2[], h3[], * hl = {h, h2, h3},
  dhdt[], dh2dt[], dh3dt[], * dhdtl = {dhdt, dh2dt, dh3dt};
h[left] = dirichlet (0);
h[right] = dirichlet (0);
h2[left] = dirichlet (0);
h2[right] = dirichlet (0);

double saw_tooth (double t) {
  double ta = fmod (t, 2.);
  return (ta < 0.5 ? 2*ta : ((1 + 0.5/1.5) - (1./1.5)*ta))/5;
}

h3[left] = dirichlet (saw_tooth(t));
h3[right] = dirichlet (0.1);

FILE * fps[3];

double k = 1;

int main() {
  init_grid (512);
  DT = 1e-3;
  run();
}

event init (t = 0) {
  /**
     Following the aforementioned video, we initialize a plucked
     guitar string from a smoothed sawtooth, using the Helmholtz
     filter.
  */
  double w = 0.05, A = 10;
  foreach() {
    h[] = x < 0.25 ? -3*x : x - 1;
    dh2dt[] = A*exp(-sq((x - 0.25)/w));
    h3[] = x/10.;
  }
  diffusion(h, 2e-4);
  /**
     Set timestep and movie stuff
  */
  dt = dtnext(DT);
 
  for (int i = 0; i < 3; i++) {
    fps[i] = popen ("gnuplot", "w");
    fprintf(fps[i],
	    "set term pngcairo\n"
	    "set xr [0: 1]\n"
	    "set yr [-1: 1]\n"
	    "set key off\n"
	    "set grid\n"
	    "set title 'String'\n"
	    "set xlabel 'x'\n"
	    "set ylabel 'h'\n");
  }
}

event stability (i++) {
  dt = dtnext(DT);
}
/**
We use a leap-frog scheme for time advancement. Noting that the `dhdt`
field lags the solution by half a time step.
 */
event advance (i++, last) {
  foreach() {
    scalar hi, dhdti;
    for (hi, dhdti in hl, dhdtl) 
      dhdti[] += dt*k*(hi[-1] - 2*hi[] + hi[1])/sq(Delta);
  }
  foreach() {
    scalar hi, dhdti;
    for (hi, dhdti in hl, dhdtl) 
      hi[] += dt*dhdti[];
  }
}

int frame = 0;
event movie(t += 0.01){
  for (int i = 0; i < 3; i++) {
    fprintf(fps[i], "set output 'plo%dt%d.png'\n",i, frame);
    fprintf(fps[i], "plot '-' w l lw 5\n");
    scalar hi = hl[i];
    foreach()
      fprintf(fps[i], "%g %g\n",x, hi[]);
    fprintf(fps[i], "e\n");
  }
  frame++;
}

event stop (t = 6){
  system("rm mov*.mp4");
  system("ffmpeg -r 60 -f image2 -i plo0t%d.png -c:v libx264 -vf format=yuv420p -y mov0.mp4");
  system("ffmpeg -r 60 -f image2 -i plo1t%d.png -c:v libx264 -vf format=yuv420p -y mov1.mp4");
  system("ffmpeg -r 60 -f image2 -i plo2t%d.png -c:v libx264 -vf format=yuv420p -y mov2.mp4");
  system("rm plo*");
  return 1;
}
