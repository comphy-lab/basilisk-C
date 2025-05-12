/**
![Water droplets may rest on a hydrophobic surface. Image via [the university of Wisconsin](https://www.chem.wisc.edu/content/hydrophobic-interactions).](https://www.chem.wisc.edu/deptfiles/Water.png)

# Equilibrium shapes of a droplet on a hydrophobic surface

An atmosphere with a density $\rho_2$ hosts a denser, watery droplet with a
density $\rho_1 = \rho_2 + \Delta\rho$. The droplet is affected by the
acceleration due to gravity ($g$) and rests on a hydrophobic
surface. A characteristic lengthscale of the droplet ($R$) may be
introduced by considering it's volume; $V=4/3 \pi R^3$. We assume that
there exists a single equilibrium state of the system that is governed
by the so-called Bond number:

$$Bo = \frac{\Delta \rho g R^2 }{\sigma}$$

On this page we aim to find the shape of the droplet for a range of
$Bo$ values. Because of symmetry reasons, we choose to use the
axi-symmetric solver to save some computational effort.
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "contact.h"
#include "tension.h"
#include "axi.h"
#include "reduced.h"
#include "view.h"

int maxlevel = 8;
double R = 1.;
vector h[];
h.t[left] = contact_angle (165*pi/180.);
u.t[left] = dirichlet (0.);
int sw = 1;
char fname[99];
int it = 0;
double AT = 1e-3;
/**
## set-up

We choose to use normalized values for all parameters that define the
Bond number except for the value of the surface tension
($\sigma$). This parameter is used to control the aforementioned
dimensionless ratio $Bo$.
 */
int main(){
  G.x = -1; 
  mu1 = 1.; 
  mu2 = 1.; 
  rho1 = 2.; 
  rho2 = 1.;
  L0 = 4*R;
  f.height = h;
  init_grid (64);
  run();
}
/**
For $Bo \rightarrow 0$, a spherical droplet shape is expected. We
therefore start with an small value for the Bond number and initialize
the droplet as a sphere.
 */
event init (t=0){
  TOLERANCE = 10E-4;
  f.sigma = 32.;
  refine (sq(x-R) + sq(y) < 1.1*R && level < maxlevel);
  fraction (f, R - sq(x - R*0.95) - sq(y));
}

/**
## Finding equilibrium solutions

Due to the effect of gravity, the droplet will deform. Deformation is
characterized by velocities and we assume that equilibrium solutions
have no kinetic energy in the system. We apply an arbritrary threshold
(`AT`) that sets a lower limit for the kinetic energy. When the
kinetic energy is lower than this value, the solution is considered to
be in equilibrium. At this point we output the shape of the bubble to
a file. Next, the value of the Bond number is increased by a factor of
two and a new equilibrium is calculated, starting from the
previous solution. Since this is not really a physically consistent set-up we also monitor the volume of the droplet to diagnose possible changes in the $R$ parameter.
 */

double e = 0.;
event check_equi (i+=25; t<100){
  double en=0.;
  foreach(reduction(+:en)){
    foreach_dimension()
      en+=sq(u.x[])*sq(Delta);
  }
  if (sw == 0){
    it++;
    if (it > 10)
      sw=1;
  }
  if (en < AT && en < e  && sw == 1){ 
    printf ("steady %g %d %g %g\n", f.sigma, i, t, en);
    sprintf(fname, "%g.fac", f.sigma);
    FILE * fp = fopen (fname, "w"); 
    output_facets (f, fp);
    fclose (fp);
    f.sigma /= 2.;
    sw = 0;
    it = 0;
  }
  e = en;
}
/**
## Visual reference

It is important to show what happends. Therefore, we plot the evolution
of the droplet as different equilibria are found.
 */

event movie (i += 100){
  char str[99];
  sprintf (str, "1/Bo = %g", f.sigma);
  clear();
  view(fov = 15, ty=-0.25, width = 800, height = 512);
  glRotatef (90, 0, 0, 1); // Is this a utility function?
  draw_vof ("f", lw=6);
  squares ("f");
  draw_string (str, pos=1, size = 25);
  mirror ({0,-1}){
    draw_vof ("f", lw=6);
    squares ("f");
  }
  save ("mov.mp4");
}
/**
The movie looks like this:
![The evolution of the droplet](rest/mov.mp4)

## Adaptation

It seems attractive to focus refinement on the bubble interface. Special care is taken to prevent errouneous energy injection thay may occur naively, as is illustrated [here](bootstrapbubble.c).
 */
event adapt (i++){
  scalar ff[];
  foreach()
    ff[] = f[];
  boundary ({ff});
  adapt_wavelet ((scalar*){ff}, (double[]){0.0001}, maxlevel);
}

/**
The run is stopped once we have found the solution for $Bo=2$.
 */
event stop (i += 100){
  if (f.sigma < 0.5)
    return 1;
}

/**
## Results

Here are some droplet shapes:

~~~gnuplot Facets of the droplets
set xr [0:1.5]
set key box top right
set size ratio -1
plot '32.fac' u 2:1  w l lw 2 t '1/Bo = 32'  ,\
     '8.fac' u 2:1 w l lw 2 t '1/Bo = 8'   ,\
     '2.fac'  u 2:1 w l lw 2 t '1/Bo = 2'   ,\
     '0.5.fac' u 2:1 w l lw 2 t '1/Bo = 0.5'
~~~
 */
