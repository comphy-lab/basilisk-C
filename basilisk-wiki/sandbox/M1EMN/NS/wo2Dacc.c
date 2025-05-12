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
//apapt=0 for msas=1, else it's nonsense
#define ADAPT 0
#define MSAS 1
#define MIN_LEVEL 3
#define LEVEL 5
#define MAX_LEVEL 6

#define tmax   10*2.*M_PI
#define alpha_w  10.

FILE * fp1;
FILE * fp2;
FILE * fp3;

/**
## main fonction*/
int main(){
  fp1 = fopen ("testpoint.csv", "w");
  fp2 = fopen ("u.csv", "w");
  fp3 = fopen ("position.csv", "w");
  L0 = 1.;
//we apply the  espace periodic boundary condition.
  periodic (right);

//we change the grid cell to study the convergence(MASA 1).
#if MSAS
	for(N= 4 ; N <= 64 ; N *=2){
  	run();
	}
#else //case normal
	N = 2 >>LEVEL;
  run();
#endif
}

//we set the no slip Boundary conditions for lateral surface
	u.n[top]  = dirichlet(0);
	u.t[top]  = dirichlet(0);

/** 
##initial event 

here we calculate the $\mu$ based on womersley number $\alpha$, we apply a cosinus force in x direction to manipulate the system periodic in time*/
event init (t = 0) {
	double viscosity = 1. / sq(alpha_w);
//  printf ("%g",viscosity);
	const face vector muc[] = {viscosity , viscosity};
	mu = muc;
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
	const face vector	g[] = {cos(t) , 0.};
	a = g;
//to check our results is "converged", we output the value in one point to see if it's periodic 
	double px = L0/2.;
	double py = 0.;
		fprintf(fp1,"%d %g %g %g \n" ,N, t, cos(t), interpolate(u.x, px, py));
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

#if MASA
/**
# Results change N
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
*/

#else
/**
# Results adapt mesh

The velocity in the middle of the domain as function of time 

~~~gnuplot testpoint
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


## Bibliography

* Womersley 1955

* Ghigho 

*/
#endif


