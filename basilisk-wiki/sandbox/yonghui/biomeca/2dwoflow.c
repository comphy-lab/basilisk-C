/**
# An axi-symmetric Womersley flow in a tube

We want to simulate the simplest model of the blood flow in an artery:
an axi-symmetric rigid periodic in time and space flow.
We started by nondimensionaliss the Navier-Stokes equation with 
$$
x = R_0 \bar x, r= R_0 \bar r, t = \tau \bar t, u_x= U_x \bar u_x, u_r= U_r \bar u_r , \rho= \rho_0$$
The pressure is decomposed in two, the effective pressure which will be 0 and the oscillating pressure 
$$P = P_0 + \hat P \bar P $$
with $\bar P = cos( \bar t)$
so the scale of pressure and velocity are related by:
 $\hat P = \frac{ U_x \rho_0 \tau }{ R_0 }$
  we also define the womersley number
$\alpha^2 = \frac{ R_0^2 }{ \tau \mu }$.

the nondimensionnal N-S forced equation becomes:
$$
\frac{\partial \bar{u}}{\partial \bar{t}} = 
-\cos\bar{t}+ \frac{1}{\alpha ^ 2} \frac{1}{\bar r} \frac{\partial}{\partial \bar r}  \bar r  \frac{\partial \bar u }{\partial \bar r}
$$
here in our case $\rho= 1$, $R_0 = 1$, $\tau =2 \pi$. Since the sytem is forced, we run the simulation in such a long time to make sure the result is converge.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "view.h"
/** if we set {ADAPT, MSAS, PBC} as {0/1, 0/1, 0/1} : 

- {0, 1, 0} default setting: We can study Mesh Size Accuracy when using acceleration BC \n
- {0, 1, 1} : We can study Mesh Size Accuracy when using pressure BC
- {0, 0, 0} : We use a fixed mesh when using acceleration BC
- {0, 0, 1} : We use a fixed mesh when using pressure BC
- {1, 0, 0} : We use a adaptive mesh when using acceleration BC
- {1, 0, 1} : We use a adaptive mesh when using pressure BC

*/ 
#define ADAPT 0 //adaptive mesh
#define MSAS 1 // Mesh Size Accuracy Study
#define PBC 0 //periodic pressure B.C
//check the logfile to make sure using the right define
#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 7

#define tmax   10*2.*M_PI
#define alpha_w  10.

FILE * fp1;
FILE * fp2;
FILE * fp3;

/**
## main fonction
all the code for "it's in XXX" is useless for runnign the code, it's just for check you have set the right situation*/
int main(){
  fp1 = fopen ("testpoint.csv", "w");
  fp2 = fopen ("u.csv", "w");
  fp3 = fopen ("position.csv", "w");
  L0 = 1.;

  #if PBC
  // in this case we use periodic pressure as the source term so no espace periodic B.C 
  fprintf(stderr,"it's in PBC \n" );
  #else
  // In this case periodic acceleration is the source term, so we apply the  espace periodic boundary condition.
  fprintf(stderr,"it's in acc BC \n" );
  periodic (right);
  #endif

  #if MSAS
  /**
  we change the grid cell to study the convergence(MASA 1).*/
  fprintf(stderr,"it's in MASA N= %d \n",N );
  for(N= 4 ; N <= 64 ; N *=2){
    fprintf(stderr,"it's in MASA N= %d \n",N );
    run();
  }
  #else //case normal
  N = 2 << LEVEL;
  fprintf(stderr,"it's in normal N= %d \n" ,N);
  run();
  #endif
}

/**
we set the no slip Boundary conditions for lateral surface*/
u.n[top]  = dirichlet(0);
u.t[top]  = dirichlet(0);

#if PBC  
/**
we set the B.C for inlet & outlet surface*/
p[left] = dirichlet(sin(t));
pf[left]  = neumann(0.);
p[right] = dirichlet(0);
pf[right]  = dirichlet(0.);
u.n[right]  = neumann(0);
u.n[left]  = neumann(0);
#endif

/** 
##initial event 

here we calculate the $\mu$ based on womersley number $\alpha$, we apply a cosinus force in x direction to manipulate the system periodic in time*/
event init (t = 0) {
  double viscosity = 1. / sq(alpha_w);
  const face vector muc[] = {viscosity , viscosity};
  mu = muc;
  #if PBC //we set the intial velocity field ,not sure it's useful
  fprintf(stderr,"it's in PBC init velo \n" );  
  foreach(){
    u.x[] = 0.;
    u.y[] = 0.;
  }
  #endif
}
/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */
#if ADAPT
event adapt (i++) {
   adapt_wavelet((scalar *){u}, (double[]){1e-3,5e-4}, MAX_LEVEL, MIN_LEVEL) ;
}
#endif

/**
## Time integration 

We change the value of the external force, and we record the pressure value at one point in each iteration to check if the result converges(periodic sinus form) after the simulation.*/
event midpressure(t <= tmax; i += 1) {
  #if PBC
  if (i == 10)
    fprintf(stderr,"it's in PBC i++ \n" );
  #else
  const face vector	g[] = {cos(t) , 0.};
  a = g;
  if (i == 10)
    fprintf(stderr,"it's in acc i++ \n" );  
  #endif 
  //to check our results is "converged" or stational, we output the value in one point to see if it's periodic 
  double px = L0/2.;
  double py = 0.;
  fprintf(fp1,"%d %g %g %g \n" ,N, t, cos(t), interpolate(u.x, px, py));
  if (t == tmax)
    fprintf(fp1,"\n \n"); //to avoid the line between each draw 
}


/**
We output the values of velocity and pressure in the middle position at the last period with time step $\delta t = 0.1 T$, so we can see the change of velocity profil, which we can compare with analytic results to check The accurecy of simulation.*/
event tracer (t <= tmax; t += (0.1* 2.*M_PI)) {
  if(t>= tmax - 2.*M_PI){
    for (double y = 0.; y < 1.; y += L0/50. ){
      double vxf = L0/2.;
      fprintf(fp2,"%d %g %g %g \n", N, t, y, interpolate(u.x , vxf, y));
    }
  }
}


event final(t = tmax){	
  foreach()
    fprintf(fp3,"%d %g %g %g \n",N, x, y, u.x[] );
}

/**
# Results N fixe
//While we use the pressure BC, le theory value is meaningless since there can be some phase differente.
//In the case of adapting mesh size (or not),  plot 1-3 is useful
~~~gnuplot convergence testpoint 
set xlabel "t"
set ylabel "u(1/2,0,t)"
plot 'testpoint.csv' us 2:4 w l 
~~~

Compare the convergence 

~~~gnuplot compare of velocity
plot 'u.csv' u 3:4 t'64',\
'utheory.csv' us 1:2 t'theory' w l lc 1
~~~

~~~gnuplot mesh grid
plot 'position.csv' u 2:3
~~~

# Results change N
//In the case of studing the convergence with mesh size,  plot 4-6 is useful
~~~gnuplot testpoint N
set format y "%g"
plot 'testpoint.csv' us ($1==4?$2:NaN):4 t'4',\
'' us ($1==16?$2:NaN):4 t'16',\
'' us ($1==64?$2:NaN):4 t'64'
~~~

~~~gnuplot compare of velocity N
plot 'u.csv' u ($1==4?$3:NaN):4 t'4' lc 6 ,\
''u ($1==16?$3:NaN):4 t'16' lc 8 ,\
''u ($1==64?$3:NaN):4 t'64' lc 4 ,\
'utheory.csv' us 1:2 t'theory' w l lc 1
~~~

~~~gnuplot mesh grid N
plot 'position.csv' u ($1==4?$2:NaN):3 t'4' lt 5 lc 1 ,\
''u ($1==16?$2:NaN):3 t'16' lt 2 lc 2 ,\
''u ($1==32?$2:NaN):3 t'32' lt 1 lc 4
~~~

## Bibliography

* Womersley 1955
* R. Ghigo Reduced-Order Models for Blood Flow in Networks of Large Arteries

*/

