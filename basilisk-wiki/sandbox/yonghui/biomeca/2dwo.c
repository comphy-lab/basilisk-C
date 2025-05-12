/**
# An axi-symmetric Womersley flow in a tube 

It is almost identical to the code [2dwoflow.c](http://basilisk.fr/sandbox/yonghui/biomeca/2dwoflow.c).
Here, we study the effects of womersley numbers and compare the results with theoretical results. Mesh adaptation is used by default.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
//difine if we want to use mesh adaption
#define ADAPT 1
// some mesh level and maxtime control
#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 7
#define tmax   10*2.*M_PI
FILE * fp1;
FILE * fp2;
FILE * fp3;
int alpha_w;

/**
## main fonction
We change the value of womersley number at each simulation, which is set in the array womun[]*/
int main(){
  periodic (right);
  N = 2 << LEVEL;
  int wonum[] = {3,10,20};
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("u.csv", "w");
  fp3 = fopen ("position.csv", "w");
  for(int K= 0 ; K <= 2 ; K++){
    alpha_w = wonum[K];
    fprintf(stderr,"# now womersley number = %d \n", alpha_w );
    run();
  }
}

/**
we set the no slip boundary conditions for lateral wall*/
u.n[top]  = dirichlet(0);
u.t[top]  = dirichlet(0);

/** 
##initial event 
here we calculate the $\mu$ based on womersley number $\alpha$, we apply a cosinus force in x direction to manipulate the system periodic in time*/
event init (t = 0) {
  double viscosity = 1. / sq(alpha_w);
  const face vector muc[] = {viscosity , viscosity};
  mu = muc;
  fprintf(stderr,"Corresponding viscosity: %g\n", viscosity);  
}

/**
## Mesh adaptation
We adapt the mesh according to the error on the volume fraction field
and the velocity. */
#if ADAPT
event adapt (i++) {
  adapt_wavelet((scalar *){u}, (double[]){5e-4,1e-3}, MAX_LEVEL, MIN_LEVEL) ;
}
#endif

/**
## Time integration 

We record the pressure value of a choosen point in order to check if the result is converged after the simulation (periodic form)*/
event midpressure(t <= tmax; i += 1) {
  const face vector  g[] = {cos(t) , 0.};
  a = g;
  double px = L0/2.;
  double py = 0.;
  fprintf(fp1,"%d %g %g %g \n" , alpha_w, t, cos(t), interpolate(u.x, px, py));
  if (t == tmax)
    fprintf(fp1,"\n \n"); //to avoid the line between each draw 
}


/**
We output the values of velocity and pressure in the middle position at the last period with time step $\delta t = 0.1 T$, so we can see the change of velocity profil, which we can compare with analytic results to check The accurecy of simulation.*/
event tracer (t <= tmax; t += (0.1* 2.*M_PI)) {
  if(t>= tmax - 2.*M_PI){
    for (double y = 0.; y < 1.; y += L0/50. ){
      double vxf = L0/2.;
      fprintf(fp2,"%d %g %g %g \n", alpha_w, t, y, interpolate(u.x , vxf, y));
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
plot 'testpoint' u 2:4 w l t'testpointvelocity'
~~~

~~~gnuplot Wo 20 compare
plot 'u.csv' u ($1==20?$3:NaN):4 t'Wo-20',\
'thwo20' us 1:2 t'theory20' w l lc 1
~~~

~~~gnuplot Wo 20 mesh grid
plot 'position.csv' u ($1==20?$2:NaN):3 t'Mesh_Wo20'
~~~


~~~gnuplot Wo 10 compare
plot 'u.csv' u ($1==10?$3:NaN):4 t'Wo-10',\
'thwo10' us 1:2 t'theory10' w l lc 1
~~~

~~~gnuplot Wo 10 mesh grid
plot 'position.csv' u ($1==10?$2:NaN):3 t'Mesh_Wo10'
~~~


~~~gnuplot Wo 3 compare
plot 'u.csv' u ($1==3?$3:NaN):4 t'Wo-3',\
'thwo3' us 1:2 t'theory3' w l lc 1
~~~

~~~gnuplot Wo 3 mesh grid
plot 'position.csv' u ($1==3?$2:NaN):3 t'Mesh_Wo3'
~~~

## Bibliography

* Womersley 1955
* R. Ghigo Reduced-Order Models for Blood Flow in Networks of Large Arteries

*/