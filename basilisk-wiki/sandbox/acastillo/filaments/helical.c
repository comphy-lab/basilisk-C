/**
# Simulated Navier-Stokes helical vortex

The evolution of a helical vortex is simulated as in
[Antoon's sandbox](http://basilisk.fr/sandbox/Antoonvh/helical.c). The goal
of the physical space initialisation is to map an analytically defined trefoil
vortex onto an Eulerian (static) numerical mesh. The helical trajectory
discussed in here is defined by:
$$
x = R \cos{(t)}
$$
$$
y = R \sin{(t)}
$$
$$
z = \frac{H}{2\pi} t
$$
where $t$ ranges between 0 and $2\pi$, $R$ is the radius and $H$ the pitch.

For this example, we use the (compressible) Navier-Stokes equations inside a
periodic box and the same parameters as in Antoon's example but fixing the
circulation $\Gamma=1$.

![Space-curve of a helical vortex relative to the domain size](helical/helical0.png)

Which results in something like this
![The movie shows a $\lambda_2=0$ isosurface](helical/helical1.mp4)

![The movie shows a $\lambda_2=0$ isosurface on top of a slice of $u_z$](helical/helical2.mp4)

*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"
#include "view.h"
#include "lambda2.h"

#define MINLEVEL 4
#define RAD (sqrt(sq(x) + sq(y)))

int n_seg;
double as = 0.1, Hs = 0.5, Rs = 1.0, n_turns;

scalar l2[];
vector omega[];

int main() {
  L0 = 8;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.1;
  N = 1<<MINLEVEL;

  n_turns = L0/Hs;
  n_seg = 64 * n_turns;

  double reynolds= 8000;
  const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds};
  mu = muc;

  periodic(back);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){1e-3, 1e-3, 1e-3}, MAXLEVEL, MINLEVEL);
}


/**
# Initial conditions
*/

#ifdef _INITIAL

double delta_init = 0.0, k_init = 0.5;

#include "helical_init.h"

event stop (i = 2)
  dump();

#else
#ifdef _RELAX
#ifdef _FRAME

event init (t = 0){

  double u_frame=0.;
  if (Hs == 4.0){
    u_frame=0.22894555168948946;
  }
  else if (Hs == 2.0){
    u_frame=0.3357585169629363;
  }
  else if (Hs == 1.0){
    u_frame=0.534953599887698;
  }

  if (pid()==0) printf( "Restoring from previous run, moving frame \n" );
  restore("./profiles/level10/0_init/dump");
  foreach(){
    u.z[] -= u_frame;
  }
  boundary ((scalar *){u});
  if (pid()==0) printf( "... Done restoring! \n" );
}

#else

event init (t = 0){
  if (pid()==0) printf( "Restoring from previous run \n" );
  restore("./profiles/level10/0_init/dump");
  if (pid()==0) printf( "... Done restoring! \n" );
}

#endif

/**
# Outputs
*/

double tsample1=0.01, tsample2=0.1, tsample3=1.0, tsample4=5.0;

#include "3d/ellipticity.h"
event end_timestep (t += tsample1, last){
  lambda2 (u, l2);
  vorticity3d(u, omega);
}

#include "helical_moments.h"
#include "helical_export.h"
#include "helical_tracers.h"

event stop (t = 1.0)
  dump();

#else
event stop (i = 5)
  dump();
#endif
#endif
