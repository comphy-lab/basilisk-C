/**
## Compare of same radius cylinder

Here we compare the womersley flow in a cylinder with same BC as the descending aorte (same mean radius), we can find the differences caused by the little geometric deformation in the real models.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h" 
#define MIN_LEVEL 6
#define LEVEL 8
#define MAX_LEVEL 9
#define tmax   5*2.*M_PI
#define radius 1.11471 
#define alpha_w 10

FILE * fp1;
FILE * fp2;
FILE * fp3;

int main(){
  size (10);
  N=64;
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("u.csv", "w");
  fp3 = fopen ("position.csv", "w");
  TOLERANCE = 1e-4;
  run();
}

/**
we set the no slip Boundary conditions for lateral wall*/
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);
u.n[right] = neumann(0.);



/** 
## initial event 
here we calculate the $\mu$ based on womersley number $\alpha$*/
event init (t = 0) {
  mask( y > radius? top : none );
  refine (y < radius && level < MAX_LEVEL);
  double viscosity = 1. / sq(alpha_w);
  const face vector muc[] = {viscosity , viscosity};
  mu = muc;
  fprintf(stderr,"Corresponding viscosity: %g\n", viscosity); 
}


/**
##  Mesh adaptation
We adapt the mesh according to the error on the volume fraction field
and the velocity. */
event adapt (i++) {
  adapt_wavelet((scalar *){u}, (double[]){5e-4,1e-3}, MAX_LEVEL,LEVEL) ;
}

/**
## Time integration 

We change the value of the external force, and we record the pressure value at one point in each iteration to check if the result converges(periodic cosus form) after the simulation.*/

event midpressure(t <= tmax; i += 1) {
  const face vector  g[] = {cos(t) , 0.};
  a = g;
  double px = L0/2.;
  double py = 0.;
  fprintf(fp1,"%d %g %g %g \n" , i, t, dt, interpolate(u.x, px, py));
  if (t == tmax)
    fprintf(fp1,"\n \n"); //to avoid the line between each draw 
}


/**
We output the values of velocity and pressure in the middle position at the last period with time step $\delta t = 0.1 T$, so we can see the change of velocity profil, which we can compare with analytic results to check The accurecy of simulation.*/
event tracer (t <= tmax; t += (0.1* 2.*M_PI)) {
  if(t>= tmax - 2.*M_PI){
    double vxf = L0/2.;
    for (double y = 0.; y < 1.5; y += L0/100. ){
      fprintf(fp2,"%g %g %g %g \n", t, vxf, y / radius, interpolate(u.x , vxf, y));
    }
  }
}

/** we output the mesh grid and velocity field at tmax*/
event final(t = tmax){	
  foreach()
    fprintf(fp3,"%d %g %g %g \n",alpha_w, x, y, u.x[] );
}


/**
# Results N fixe
~~~gnuplot convergence testpoint 
plot 'testpoint' us 2:4 w l 
~~~

Compare the convergence 

~~~gnuplot compare of velocity
plot [0:1.05][] 'u.csv' u ($4 < 10 ? $3:NaN ):($4*9.) w l t'64',\
'../2dwo/thwo10' us 1:2 t'theory' w l ,\
'../3dendotube/test-1' us 2:3 t't1' w lp,\
'../3dendotube/test-2' us 2:3 t't2' w lp,\
'../3dendotube/test-3' us 2:3 t't3' w lp,\
'../3dendotube/test-4' us 2:3 t't4' w lp,\
'../3dendotube/test-5' us 2:3 t't5' w lp,\
'../3dendotube/test-6' us 2:3 t't6' w lp,\
'../3dendotube/test-7' us 2:3 t't7' w lp,\
'../3dendotube/test-8' us 2:3 t't8' w lp,\
'../3dendotube/test-9' us 2:3 t't9' w lp,\
'../3dendotube/test-10' us 2:3 t't10' w lp
~~~

## Bibliography

* Womersley 1955
* Ghigo 

*/

