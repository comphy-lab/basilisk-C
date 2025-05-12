/**
# Weird behaviour of the event inheritance

To illustrate the problem, we aim to simulate the time evolution of an
elliptical bubble positioned at a certain distance from a flat
interface. Initially, the fluid is at rest. The surface tension will
act to minimize the bubble's surface area from the very beginning,
setting the fluid into motion. Therefore, the libraries to be used
are:
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#define L0 10
#define L1 0.2
#define R1 1.

uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;

int LEVEL;

int main() {
  size (10 [0]);
  DT = HUGE [0];
  LEVEL = 9;
  origin (-L0/2., 0.);
  init_grid (1 << (5));

  /**
  Set the densities and viscosities of the phases... */
  
  rho2 = 1e-3, mu2 = 0.01/100;
  rho1 = 1., mu1 = 0.01;
  f.sigma = 1.;
  run();
}

double geometry(double x, double y) {

  double C1 = sq((x + R1 + L1)/R1) + sq(y/(1.5*R1)) - 1.0;
  double D1 = -x - 1e-8;

  return min(C1, D1);
} 

/**
Initialize the geometry and grid... */

event init (t = 0) {  
  double eps = 0.05;
  refine ( sq((x + R1 + L1)/(R1-eps)) + sq(y/(1.5*(R1-eps))) > 1 && sq((x + R1 + L1)/(R1+eps)) + sq(y/(1.5*(R1+eps))) < 1. && level < LEVEL);
  refine (x < eps && x > -(2*0.1+eps) && level < LEVEL);
  fraction(f, geometry(x,y));
}

/**
Use adaption also... */

event adapt (i++) {
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1}, (double[]){1.e-3}, minlevel = 7, maxlevel = LEVEL);
}

/**
But if you wish to monitorize for example the timestep set by the
stability events at a certain instant (t=0.2 in this case) using the
event inheritance (vof for instance) It seems to ruin the control of
the iterations of the event erasing the event in the first iterations . */

#if 0
event vof (t = 0.2; i += 10; t <= 0.65) {
printf("vof: %g %g %d\n",t, dt, i);
}
#endif

event abrupt (i = 1) {
  return 1;
}

/**
# Results

With the vof event above activated the vof events in the .h files are not longer
set in each iteration step

~~~ bash
 events (i = 0, t = 0)    
  ...                                                                                        
  stability                 src/tension.h:36                                      
  stability                 src/navier-stokes/conserving.h:72           
  stability                 src/vof.h:140
  stability                 src/navier-stokes/centered.h:224
  set_dtmax                 src/navier-stokes/centered.h:222
  stability                 src/tension.h:36
  stability                 src/navier-stokes/conserving.h:72
  stability                 src/vof.h:140
  stability                 src/navier-stokes/centered.h:224
~~~
<div style="background-color: #add8e6; padding: 10px; border-radius: 5px;">
<pre><code class="language-bash">
  vof                       src/navier-stokes/conserving.h:117
  vof                       src/vof.h:380
  vof                       src/navier-stokes/centered.h:234
</code></pre>
</div>
~~~bash
  tracer_advection          src/navier-stokes/conserving.h:192
  tracer_advection          src/two-phase-generic.h:50
  tracer_advection          src/navier-stokes/centered.h:235
  tracer_diffusion          src/navier-stokes/centered.h:236
  ...
~~~
*/
