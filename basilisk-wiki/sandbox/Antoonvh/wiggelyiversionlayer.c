/**
# A wiggely inversion layer

We simulate a sharp temperature/density/buoyancy inversion that
wiggles a bit.  This example aims to be a minimalistic example of the
generation of errouneous vortex structures as was observed in a much
more complex simulation of the atmospheric diurnal cycle.A possible 
solution is also presented. 
*/ 

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "fractions.h"

int maxlevel  = 9;
scalar b[];
scalar * tracers = {b};
face vector av[];
int j = 0;

/**
   We run three cases,

1. Using the default settings (which shows strange behaviour)
2. Using a slope limiter for the advection of the buoyancy parameter
3. Using some more non-default settings that seem important for atmospheric simulations (for other reasons). 

Case 1, 2 and 3 are labelled with $j = 0, 1$ and $2$ respectively. 

*/

int main(){
  periodic(left);
  L0 = M_PI;
  a = av;
  run();
  j++;
  b.gradient = minmod2;
  run();
  j++;
  p.refine = p.prolongation = refine_linear;
  foreach_dimension()
    u.x.refine = refine_linear;
  run();
}

event init(t = 0){
  init_grid(1 << maxlevel);
  DT = 0.05; 
  fraction(b, (y - (1. + 0.02*sin(2.*x))));
}

event acceleration(i++){
  coord gr = {0,1,0};
  foreach_face(y)
    av.y[]+=(gr.y*(b[]+b[0,-1])/2.);
}

event adapt (i++)
  adapt_wavelet((scalar*){b,u}, (double[]){0.0125, 0.0125, 0.0125}, maxlevel);

/**
## Output

We only analyze movies
 */
event movies (t += 0.2 ; t <= 50){
  char fname[99];
  scalar lev[], omega[];
  vorticity(u, omega);
  foreach()
    lev[] = level;
  if (j == 0){
    sprintf (fname, "ppm2mp4 b%d.mp4", j);
    static FILE * fp = popen (fname, "w");
    output_ppm(b, fp, n = 512, min = 1.49, max = -1.51);
    sprintf (fname, "ppm2mp4 lev%d.mp4", j);
    static FILE * fp3 = popen (fname, "w");
    output_ppm(lev, fp3, n = 512, min = 1, max = maxlevel);
  }
  sprintf (fname, "ppm2mp4 omg%d.mp4", j);
  static FILE * fp2 = popen (fname, "w");
  output_ppm(omega, fp2, n = 512, min = -1, max = 1);
}

/**
## Results:

For all cases, the evolutions of the buoyancy field looks a lot
like the result that is obtained for $j = 0$:

![The buoyancy field for a wiggely inversion layer ($j = 0$)](wiggelyiversionlayer/b0.mp4)

Nothing seems wrong. However, the evolution of the grid structure shows some
unexpected behaviour;

![The used grid in (secret?) color code ($j = 0$)](wiggelyiversionlayer/lev0.mp4)

I happen to know from plotting the [$\lambda_2$ iso-surfaces](http://www.basilisk.fr/src/lambda2.h) for the 3D
simulation (at some more expense) that there may be errouneous vortex
structures. Therefore, we plot the vorticity field:

![Vorticity field ($j = 0$)](wiggelyiversionlayer/omg0.mp4)

A Oei!, something is clearly wrong. Maybe the usage of a slope limiter for `b` helps.

![Vorticity field ($j = 1$)](wiggelyiversionlayer/omg1.mp4)

The most prominent issues seem mostly gone now. Great succes! Before we redirect our
attention back to the 3D case, we check if the conservative attributes for
the pressure inteprolation and momentum refinement do not ruin the
results.

![Vorticity filed ($j = 2$)](wiggelyiversionlayer/omg2.mp4)

Results are still OK-ish. Any remaining issues will be "diffused" away in a more representative set-up for the atmospheric boundary layer.

## Note

If you are interested in the effects of gravity and sharp density interfaces, you may want to
consider using the [reduced gravity approach](http://www.basilisk.fr/src/reduced.h).
*/