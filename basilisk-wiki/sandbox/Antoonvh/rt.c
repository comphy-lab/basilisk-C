/**
# The Rayleigh-Taylor instability and the 'Black body' colourbar

Inspired by the reference for the so-called 'cool-warm' colourbar, I
visited the website of [Kenneth
Moreland](http://www.kennethmoreland.com/color-advice/) and found that
the 'black body' colour bar is simply stunning. We examplify it with
the Rayleigh-Taylor instability:

![A nice colourbar not only helps with the interpretation of data, it can also help
 expand the artistic freedom for visualizations](rt/rt.png)

This looks an awfull lot like [this $\mu$-HH based
video](https://vimeo.com/84675235), which puts the artistic 'freedom'
argument into its perspective.
*/
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

scalar b[];
scalar * tracers = {b};
face vector av[];
int maxlevel = 11;
long unsigned int n;
/**
## Implementation

We follow Stephane's implementation of the 'cool-warm' diverging
colour map and use the colour codes from Kenneth's website.
 */
void black_body (double cmap[NCMAP][3])
{
  /* black body color map from:
   * http://www.kennethmoreland.com/color-advice/
   */
  static double basemap[33][3] = {
    {0.0,0.0,0.0},
    {0.0857913205762,0.0309874526184,0.0173328711915},
    {0.133174636606,0.0588688899571,0.0346802666087},
    {0.180001956037,0.0730689545154,0.0515393237212},
    {0.22981556179,0.0840603593119,0.0647813713857},
    {0.281397607223,0.093912584278,0.075408501413},
    {0.334521638801,0.102639499627,0.0842454688083},
    {0.388957802186,0.110254429637,0.0927990674821},
    {0.444611925648,0.116732501721,0.101402659637},
    {0.501422312285,0.122025816585,0.110058408122},
    {0.559331322331,0.126067584009,0.118767796491},
    {0.618285970576,0.128767919785,0.127531801155},
    {0.678237857955,0.130007052818,0.136351016263},
    {0.712849583079,0.181721849923,0.13081678256},
    {0.743632057947,0.232649759358,0.120991817028},
    {0.774324938583,0.279315911516,0.108089917959},
    {0.804936242903,0.323627020047,0.0907961686083},
    {0.835473266757,0.366524681419,0.0662363460741},
    {0.865942668698,0.408541395043,0.026029485466},
    {0.876634426153,0.46401951695,0.0173065426095},
    {0.883455346031,0.518983528803,0.0149628730405},
    {0.88905246237,0.572164381169,0.013499801006},
    {0.893375939063,0.624108797455,0.0130334871745},
    {0.89637036663,0.675180034619,0.013680092215},
    {0.897973818846,0.725630730259,0.015555776796},
    {0.898116710502,0.775642817733,0.0187767015864},
    {0.896720396485,0.825350944866,0.023459027255},
    {0.927670131094,0.859991226192,0.319086199143},
    {0.956158602738,0.893933112845,0.503316730316},
    {0.97827065392,0.92856476667,0.671307024002},
    {0.993196411712,0.963913323002,0.83560909192},
    {1.0,1.0,1.0},
  };
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(31 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

/**
## Simulating the Rayleigh-Taylor instability

All steps are rather straight forward.
 */

int main(){
  periodic(left);
  a = av; //Link gravity acceleration
  init_grid (32);
  b.gradient = minmod2; //Sharp interface advection
  const face vector muc[] = {1E-6, 1E-6};
  mu = muc;
  run();
}

event init (t = 0){
  DT = 0.01;
  do{ // Initialize an unstable stratification with a sharp interface.
    foreach(){
      b[] = (y < Y0 + L0/2.);
      u.y[] = 0.001*noise(); // Add noise to kick-off the growth of the instability
    }
  }while (adapt_wavelet({b}, (double[]){0.01}, maxlevel).nf);
}

event acceleration (i++){ //The effect of Gravity
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++)
  diffusion (b, dt, mu);

event adapt (i++)
  adapt_wavelet ((scalar*){b,u}, (double[]){0.01, 0.01, 0.01}, maxlevel);

event output (t += 0.01; t <= 1.001){ // Output a .mp4 movie and a .png image, then stop the simulation
  boundary ({b});
  scalar db[];
  foreach(){
    db[] = 0.;
    foreach_dimension()
      db[] += sq((b[1] - b[-1])/(2*Delta));
    if (db[] > 0.)
      db[] = log (sqrt (db[]) + 1.);
  }
  boundary({db});
  output_ppm(db, file = "rt.mp4", n = (1 << (maxlevel)), map = black_body,
	     linear = true, box = {{0., 0.45}, {1., 0.55}}, min = 0., max = 6.);
  if (t > 0.999){
    output_ppm(db, file = "rt.png", n = (1 << (maxlevel)), map = black_body,
	       linear = true, box = {{0., 0.45}, {1., 0.55}}, min = 0., max = 6.);
    foreach()
      db[] *= - 1;
    boundary({db});
    output_ppm(db, file = "rt_banner.png", n = (1 << (maxlevel)), map = black_body,
	       linear = true, box = {{0., 0.47}, {1., 0.53}}, min = -6.5, max = -1);
  }
}
/**
The banner with inverted colors looks like this:

![](rt/rt_banner.png)
*/

event stop (t = 1.001){
  return 1;
}
  
/**
The movie reveals how the solution got into the state depicted above. The original rendering is 2048 pixels wide, we display it with half (or 1/4-th) of the pixels.

<video width="1024" height="154" controls>
<source src="rt/rt.mp4" type="video/mp4">
</video> 

 */