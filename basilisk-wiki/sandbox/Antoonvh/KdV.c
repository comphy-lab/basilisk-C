/**
![Solitary waves are a source of happiness to a broad audiance. Image
 via [ESA
 science](http://sci.esa.int/cluster/42436-soliton-in-water/).](https://sci.esa.int/documents/34817/35562/1567218782764-CTS-070308-image1-410b.jpg)

# Soliton Solutions to the Korteweg-De Vries Equation

So-called soliton solutions may exist as a special case to special
equations. This entails that a *solution* shape, while it moves
through space, does not deform over time. Meaning that observations of
such phenomenon may be expected to be either very rare (special
solutions to special equations) or very common (since they persist
over time for many to observe). In this example we follow the soliton
solution as listed by
[wikipedia](https://en.wikipedia.org/wiki/Korteweg%E2%80%93de_Vries_equation)
page on the Korteweg-De Vries equation. Fortunately, there is a solver
to this equation in my sandbox.

*/
#include "grid/multigrid1D.h"
#include "KdV.h"
#include "run.h"

/**
We study two soliton solutions, both described by: 

$$c(x,t) = -\frac{v}{2 \mathrm{cosh}^2\left(\frac{\sqrt{v}}{2}(x - vt - b)\right)}$$

for $v = 1$ and $v = 2$. To omit the most prominent issues at the
boundaries, we place them far away from the region of interest.
 */
scalar fast[], slow[], both[];

int main(){
  L0 = 50;
  X0 = -L0/2;
  init_grid(512);
  run();
}
/**
## Initialization

We initialize the solution fields and we also open a pipeline for output
plots with `gnuplot`.
 */
FILE * gnuplotPipe;
event init(t = 0){
  double v = 2;
  foreach() {
    fast[] = -v/(2.*sq(cosh(pow(v, 0.5)*0.5*(x + 12))));
    slow[] = -1./(2.*sq(cosh(0.5*(x + 3))));
    both[] = fast[] + slow[];
  }
  DT = 0.01;
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
        "set term pngcairo\n"
        "set xr [-15: 25]\n"
        "set yr [-1.8: 0.1]\n"
        "set key bottom left\n"
        "set grid\n"
        "set title 'KdV Solitons'\n"
        "set xlabel 'x'\n"
        "set ylabel 'c'\n");
}

/**
## Time integration

User convienience is provided via the `KdV()`-fuction interface.
 */
event advance(i++) { 
  dt = dtnext (DT);
  KdV({fast, slow, both}, dt);
}
/**
## Output many plots

Each two time steps we output a plot that is produced via the pipe. 
 */
int frame = 0;
event movie(i += 2){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot \
          '-' w l lw 5 t 'fast', \
          '-' w l lw 5 t 'slow', \
          '-' w l lw 3 t 'both'\n");
  for (scalar s in {fast, slow, both}) {
    foreach()
      fprintf(gnuplotPipe, "%g %g\n",x, s[]);
    fprintf(gnuplotPipe, "e\n");
  }
  frame++;
}
/**
## Creation of a movie from plots

We use Gnu/Linux commands and `ffmpeg` to generate a movie
from all our plots.
 */
event stop (t = 15){
  system("rm mov.mp4");
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  return 1;
}

/**
## Result
We can check if the solutions are indeed solitons. 

![Well done `KdV()` function](KdV/mov.mp4)

## Addional analysis

We check the mean position and the variance of the solitons.

~~~gnuplot The moments of overtaking coincide
set xlabel 'time'
set ylabel 'position'
set grid
plot 'out' u 1:2 t 'fast' , '' u 1:4 t 'slow' , '' u 1:6 t 'both'
~~~

~~~gnuplot The evolution of the variance
set xlabel 'time'
set ylabel 'Variance'
set grid
plot 'out' u 1:3 t 'fast' , '' u 1:5 t 'slow' , '' u 1:7 t 'both'
~~~


 */
event positions (t += 0.1) {
  printf ("%g ", t);
  for (scalar s in {fast, slow, both}) {
    double xp = 0, sm = 0, mean;
    foreach() {
      xp += s[]*x*Delta;
      sm += s[]*Delta;
    }
    printf ("%g ", mean = xp/sm);
    xp = 0;
    foreach() 
      xp += s[]*sq(x - mean)*Delta;
    printf ("%g ", xp/sm);
  }
  putchar ('\n');
}

/**
## See also

* [A Vortex Soliton Colliding with a Wall](lamb-dipole.c)
*/
